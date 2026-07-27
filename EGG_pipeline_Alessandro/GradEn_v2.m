function H = GradEn_v2(img, params)
% GRADEN  Gradient Entropy of a 2D image (grayscale matrix).
%
%   H = GRADEN(img) computes the normalized Gradient Entropy (GradEn) of the
%   2D array IMG, the two-dimensional extension of Slope Entropy proposed by
%   Jiang & Shang, "Gradient entropy (GradEn): the two dimensional version of
%   slope entropy for image analysis", arXiv:2502.18516.
%
%   H = GRADEN(img, params) overrides the default symbolization thresholds.
%
%   INPUTS
%     img    : 2D grayscale image (nRows x nCols). RGB input is converted to
%              grayscale. The image is treated as a real-valued matrix.
%     params : (optional) struct with the symbolization thresholds:
%                - delta : lower threshold (small-gradient / "flat" cutoff)
%                - gamma : upper threshold (large-gradient cutoff)
%              If omitted, the fixed standard-normal quantile thresholds
%              recommended in the reference are used (see below).
%
%   OUTPUT
%     H      : Normalized Gradient Entropy in [0, 1]. H -> 1 for maximally
%              irregular textures (uniform pattern distribution), H -> 0 for
%              highly regular ones.
%
%   METHOD (following the reference)
%     1. For each 2x2 pixel block, three gradients are computed: horizontal
%        (Gh = x(i,j+1)-x(i,j)), vertical (Gv = x(i+1,j)-x(i,j)) and diagonal
%        (Gd = x(i+1,j+1)-x(i,j)), forming a 3-D gradient vector per block.
%     2. Each gradient direction is z-score standardized across all blocks.
%     3. Every standardized value is mapped to a symbol in {-2,-1,0,1,2}
%        using the thresholds delta and gamma:
%            a <= -gamma        -> -2
%           -gamma < a <= -delta-> -1
%           -delta < a <=  delta->  0
%            delta < a <=  gamma->  1
%            a >  gamma         ->  2
%     4. Each 3-symbol block is encoded as a base-5 pattern (5^3 = 125
%        possible patterns), and the pattern probability distribution p(pi_k)
%        is estimated by relative frequency (Eq. 1 of the reference).
%     5. GradEn is the Shannon entropy of that distribution normalized by
%        log(5^3), its theoretical maximum (Eq. 2 of the reference).
%
%   THRESHOLDS
%     After z-score standardization the gradient values are approximately
%     standard normal, so the reference fixes delta and gamma as quantiles of
%     N(0,1): delta = Phi^{-1}(a), gamma = Phi^{-1}(b), with a = 0.55 and
%     b = 0.80 (delta ~ 0.1257, gamma ~ 0.8416). These are constants, i.e.
%     the same reference levels for every image, which keeps GradEn values
%     comparable across images/subjects.
%
%   Reference:
%     R. Jiang, P. Shang, "Gradient entropy (GradEn): the two dimensional
%     version of slope entropy for image analysis", arXiv:2502.18516.
%     Based on Slope Entropy: D. Cuesta-Frau, "Slope entropy: a new time
%     series complexity estimator...", Entropy, 2019.

    % --- 1. Preprocessing -------------------------------------------------
    if size(img, 3) == 3
        img = rgb2gray(img);
    end
    img = double(img);

    % --- 2. Directional gradients over each 2x2 block ---------------------
    Gh = img(1:end-1, 2:end)   - img(1:end-1, 1:end-1);   % horizontal
    Gv = img(2:end,   1:end-1) - img(1:end-1, 1:end-1);   % vertical
    Gd = img(2:end,   2:end)   - img(1:end-1, 1:end-1);   % diagonal

    % --- 3. Per-direction z-score standardization -------------------------
    G     = [Gh(:), Gv(:), Gd(:)];
    mu    = mean(G, 1);
    sigma = std(G, 0, 1);
    sigma(sigma == 0) = eps;               % avoid division by zero (flat image)
    GS = (G - mu) ./ sigma;

    % --- 4. Symbolization thresholds delta and gamma ----------------------
    % Default: fixed standard-normal quantiles as in the reference (a = 0.55,
    % b = 0.80). norminv(p) = sqrt(2)*erfinv(2p-1), so no Statistics Toolbox
    % is required.
    if nargin < 2 || ~isfield(params, 'delta')
        a = 0.55;   b = 0.80;
        delta = sqrt(2) * erfinv(2*a - 1);   % ~ 0.1257
        gamma = sqrt(2) * erfinv(2*b - 1);   % ~ 0.8416
    else
        delta = params.delta;
        gamma = params.gamma;
    end
    if delta >= gamma                        % enforce monotonicity
        delta = gamma - eps;
    end

    % --- 5. Map to symbols {-2,-1,0,1,2} ----------------------------------
    edges   = [-Inf, -gamma, -delta, delta, gamma, Inf];
    symbols = [-2, -1, 0, 1, 2];
    symb_map = zeros(size(GS));
    for k = 1:3
        symb_map(:, k) = discretize(GS(:, k), edges, symbols);
    end

    % --- 6. Encode each 3-symbol block as a base-5 pattern (0..124) -------
    symb_shift = symb_map + 2;               % shift {-2..2} -> {0..4}
    patterns = symb_shift(:, 1)*25 + symb_shift(:, 2)*5 + symb_shift(:, 3);

    % --- 7. Pattern probability distribution and normalized entropy -------
    counts = histcounts(patterns, 0:125);    % 125 bins for the 5^3 patterns
    prob = counts / sum(counts);
    prob = prob(prob > 0);
    % Normalize by log(5^3) = the maximum entropy over 125 patterns -> [0,1].
    H = -sum(prob .* log(prob)) / log(5^3);
end