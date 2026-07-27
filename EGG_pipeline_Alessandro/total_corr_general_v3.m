function TC = total_corr_general_v3(signals, alpha)
%TOTAL_CORRELATION  Total correlation via matrix-based Renyi alpha-entropy.
%
%   TC = TOTAL_CORRELATION(SIGNALS, ALPHA) estimates the total correlation
%   (also called multi-information) among two or more scalar channels:
%
%       TC(X1,...,Xn) = sum_i S_alpha(Xi)  -  S_alpha(X1,...,Xn)
%
%   where S_alpha is the matrix-based Renyi alpha-order entropy of
%   Sanchez Giraldo et al. [2], and the joint entropy follows the
%   multivariate extension of Yu et al. [3], based on the Hadamard
%   (element-wise) product of the normalized kernel Gram matrices.
%   No probability density is ever estimated: all quantities are computed
%   directly from Gram matrices built on the samples.
%
%   INPUTS
%     signals : [N x nCh] matrix whose columns are the channels, or a
%               cell array {x1, x2, ...} of vectors of equal length N.
%     alpha   : Renyi entropy order (alpha > 0, alpha ~= 1).
%
%   OUTPUT
%     TC      : scalar total correlation in bits.
%               TC >= 0 up to finite-sample estimation error, and
%               TC = 0 if and only if the channels are jointly
%               independent (in the infinite-sample limit).
%
%   IMPLEMENTATION NOTES
%     ALPHA = 2 (recommended default) uses the EXACT identity
%         sum_i lambda_i(A)^2 = trace(A^2) = ||A||_F^2,
%     which requires only sums of squared kernel values. These sums are
%     accumulated in row blocks, so the full N x N Gram matrices are
%     NEVER stored: memory is O(blockRows * N) regardless of N, and the
%     result is bit-identical to the full-matrix computation. This makes
%     N of the order of 10^5 feasible (time still scales as O(nCh*N^2)).
%
%     Generic ALPHA requires the full eigenspectrum, i.e. storing an
%     N x N matrix and an O(N^3) symmetric eigendecomposition. For large
%     N this is intractable by nature (not by implementation); consider
%     ALPHA = 2, or subsample/average over blocks of samples.
%
%   COMPLEXITY
%     alpha = 2     : O(nCh * N^2) time,  O(blockRows * N) memory.
%     generic alpha : O(nCh * N^3) time,  O(N^2) memory.
%
%   EXAMPLE
%     x = randn(1000,1);  y = x + 0.3*randn(1000,1);  z = randn(1000,1);
%     TC = total_correlation([x y z], 2);   % driven by the (x,y) coupling
%
%   REFERENCES
%     [1] S. Watanabe, "Information theoretical analysis of multivariate
%         correlation," IBM Journal of Research and Development, 4(1),
%         66-82, 1960.                                (total correlation)
%     [2] L.G. Sanchez Giraldo, M. Rao, J.C. Principe, "Measures of
%         entropy from data using infinitely divisible kernels," IEEE
%         Trans. Information Theory, 61(1), 535-548, 2015.
%     [3] S. Yu, L.G. Sanchez Giraldo, R. Jenssen, J.C. Principe,
%         "Multivariate extension of matrix-based Renyi's alpha-order
%         entropy functional," IEEE TPAMI, 42(11), 2960-2966, 2020.
%
%   See also EIG, NORM.

% -------------------------------------------------------------------------
% 1) Input handling and validation
% -------------------------------------------------------------------------
if ~iscell(signals)
    if ~ismatrix(signals) || ~isnumeric(signals)
        error('total_correlation:input', ...
              'SIGNALS must be an [N x nCh] numeric matrix or a cell array of vectors.');
    end
    signals = num2cell(signals, 1);          % columns -> 1 x nCh cell array
else
    signals = cellfun(@(x) x(:), signals, 'UniformOutput', false);
end

nCh = numel(signals);
if nCh < 2
    error('total_correlation:input', 'At least 2 channels are required.');
end

nSamp = numel(signals{1});
if any(cellfun(@numel, signals) ~= nSamp)
    error('total_correlation:input', 'All channels must have the same length.');
end

if ~isscalar(alpha) || alpha <= 0 || alpha == 1
    error('total_correlation:input', ...
          'ALPHA must be a positive scalar different from 1 (Shannon limit not implemented).');
end

% -------------------------------------------------------------------------
% 2) Kernel bandwidths (Silverman rule of thumb, 1-D case)
% -------------------------------------------------------------------------
%   sigma_i = 1.06 * std(x_i) * N^(-1/5)
sigma = zeros(1, nCh);
for i = 1:nCh
    sigma(i) = 1.06 * std(signals{i}) * nSamp^(-1/5);
    if sigma(i) == 0                         % constant channel
        sigma(i) = eps;                      % avoids 0/0 = NaN in the exponent
    end
end

if alpha == 2
    TC = tc_alpha2_blocked(signals, sigma, nSamp, nCh);
else
    TC = tc_generic_alpha(signals, sigma, nSamp, nCh, alpha);
end
end

% =========================================================================
function TC = tc_alpha2_blocked(signals, sigma, N, nCh)
%TC_ALPHA2_BLOCKED  Exact alpha = 2 total correlation without storing Gram
% matrices.
%
% With the Gaussian kernel k_i(j,l) = exp(-(x_j-x_l)^2 / (2*sigma_i^2)) and
% the unit-trace normalization A_i = G_i / N:
%
%   marginal:  tr(A_i^2) = (1/N^2) * sum_jl k_i(j,l)^2
%   joint   :  J = A_1 o ... o A_n (Hadamard), trace(J) = N^(1-n) exactly
%              (unit-diagonal kernels), and after renormalization
%              tr(Jn^2) = sum_jl [prod_i k_i(j,l)^2] / N^2
%
%   S_2(A)  = -log2(tr(A^2))
%   TC      = sum_i S_2(A_i) - S_2(Jn)
%
% All quantities are plain sums over sample pairs (j,l), so they are
% accumulated over blocks of rows: peak memory is 2 * blockRows * N
% doubles, independent of the total number of samples. The result is
% identical (to machine precision) to the full-matrix computation.

% Block size: cap each temporary [blockRows x N] matrix at ~512 MB.
blockRows = max(1, min(N, floor(512e6 / 8 / N)));

sumMarg  = zeros(1, nCh);                    % sum_jl k_i(j,l)^2 per channel
sumJoint = 0;                                % sum_jl prod_i k_i(j,l)^2

for r0 = 1:blockRows:N
    rows = r0:min(r0 + blockRows - 1, N);
    P = ones(numel(rows), N);                % running Hadamard product
    for i = 1:nCh
        x  = signals{i};
        % Squared kernel: k^2 = exp(-d^2 / sigma^2)  (note: sigma^2, not
        % 2*sigma^2, because the kernel value is squared).
        K2 = exp(-(x(rows) - x.').^2 ./ sigma(i)^2);
        sumMarg(i) = sumMarg(i) + sum(K2(:));
        P = P .* K2;
    end
    sumJoint = sumJoint + sum(P(:));
end

TC = sum(-log2(sumMarg ./ N^2)) + log2(sumJoint / N^2);
end

% =========================================================================
function TC = tc_generic_alpha(signals, sigma, N, nCh, alpha)
%TC_GENERIC_ALPHA  Generic-order path: requires the full eigenspectrum.

% Fail early with a helpful message instead of exhausting memory.
bytesNeeded = 2 * 8 * N^2;                   % A and the Hadamard accumulator
if bytesNeeded > 16e9
    error('total_correlation:memory', ...
         ['N = %d requires ~%.1f GB for the Gram matrices and an O(N^3) ', ...
          'eigendecomposition. Use ALPHA = 2 (memory-safe exact path), ', ...
          'or subsample / average TC over blocks of samples.'], ...
          N, bytesNeeded / 1e9);
end

sumMarginals = 0;
J = [];                                      % joint (Hadamard) accumulator
for i = 1:nCh
    x = signals{i};
    G = exp(-(x - x.').^2 ./ (2 * sigma(i)^2));
    A = G ./ N;                              % unit trace (G_jj = 1)
    sumMarginals = sumMarginals + renyi_entropy_gram(A, alpha);
    if isempty(J), J = A; else, J = J .* A; end
end

J = J ./ trace(J);                           % renormalize the joint [3]
TC = sumMarginals - renyi_entropy_gram(J, alpha);
end

% =========================================================================
function S = renyi_entropy_gram(A, alpha)
%RENYI_ENTROPY_GRAM  Matrix-based Renyi entropy of a unit-trace Gram matrix.
%
%   S_alpha(A) = 1/(1-alpha) * log2( sum_i lambda_i(A)^alpha )
%
% (A+A')/2 removes round-off asymmetry so EIG takes the symmetric
% (real-spectrum) code path; tiny negative eigenvalues produced by
% floating-point noise are clipped and the spectrum renormalized.
lambda = eig((A + A.') / 2, 'vector');
lambda = max(lambda, 0);
lambda = lambda ./ sum(lambda);
S = (1 / (1 - alpha)) * log2(sum(lambda.^alpha));
end