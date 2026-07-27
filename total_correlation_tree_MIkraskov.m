function TC = total_correlation_tree_MIkraskov(data_cells)
% TOTAL_CORRELATION_TREE_MIKRASKOV 
% Estimates Total Correlation using a tree-based decomposition and
% Kraskov mutual information (k-nearest-neighbor estimator). 
% Based on Bai et al., 2023
%
% This function handles time series that may contain NaN values by:
%   - extracting the longest valid (non-NaN) contiguous segment for each variable
%   - aligning all variables on the common valid temporal window
% This approach preserves the dynamical structure of the data and avoids
% artificial discontinuities in the reconstructed phase space.
%
% INPUT:
%   data_cells : cell array {1×D}, each element is 1xN or Nx1 vector
%
% OUTPUT:
%   TC : scalar estimate of the Total Correlation
%
% REFERENCES 
% Bai, K., Cheng, P., Hao, W., Henao, R., & Carin, L. (2023, April). 
% Estimating total correlation with mutual information estimators. 
% In International Conference on Artificial Intelligence and Statistics (pp. 2147-2164). PMLR.
%

%% --------------------------
% 0. Basic checks
%% --------------------------
if ~iscell(data_cells)
    error('Input must be a cell array {x1,x2,...}.');
end

D = numel(data_cells);
for i = 1:D
    data_cells{i} = data_cells{i}(:)';  % force row vector
end

%% --------------------------
% 1. Extract the longest Nan-free segment for each variable
%% --------------------------
valid_segments = cell(1,D);
seg_lengths    = zeros(1,D);

for i = 1:D
    x = data_cells{i};

    % find contiguous non-NaN segments
    isValid = ~isnan(x);
    d = diff([false, isValid, false]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;

    % find the longest segment
    [seg_lengths(i), idx] = max(ends - starts + 1);
    if seg_lengths(i) < 50
        error('Variable %d has insufficient valid data for embedding.', i);
    end

    valid_segments{i} = x(starts(idx):ends(idx));
end

%% --------------------------
% 2. Align all variables to their common valid temporal window
%% --------------------------
min_len = min(seg_lengths);
for i = 1:D
    valid_segments{i} = valid_segments{i}(1:min_len);
end

%% replace data_cells with cleaned aligned data
data_cells = valid_segments;
N = min_len;

%% --------------------------
% 3. Estimate embedding dimension & delay
%% --------------------------
dims = zeros(1,D);
taus = zeros(1,D);

for i = 1:D
    [~, tau_i, m_i] = phaseSpaceReconstruction(data_cells{i});
    dims(i) = m_i;
    taus(i) = tau_i;
end

%% --------------------------
% 4. Configure MI estimator
%% --------------------------
cfg.metric  = 'euclidean';
cfg.verbose = 0;
cfg.mass    = min(max(round(sqrt(N)/10),3),10); % k for Kraskov
cfg.dims    = [2 2];  % placeholder, updated per MI computation
cfg.taus    = [1 1];

%% --------------------------
% 5. Build tree decomposition scheme
%% --------------------------
scheme = build_tree_scheme(1, D); % rows = [left midpoint right]
TC = 0;

%% --------------------------
% 6. Compute MI across tree splits
%% --------------------------
for k = 1:size(scheme,1)
    L = scheme(k,1);
    M = scheme(k,2);
    R = scheme(k,3);

    % build left block
    Xleft = [];
    for i = L:(M-1)
        Xleft = [Xleft; data_cells{i}(:)'];
    end

    % build right block
    Xright = [];
    for i = M:(R-1)
        Xright = [Xright; data_cells{i}(:)'];
    end

    if ~isempty(Xleft) && ~isempty(Xright)
        % determine embedding parameters for each block
        dimsL = max(dims(L:M-1));
        dimsR = max(dims(M:R-1));
        tausL = max(taus(L:M-1));
        tausR = max(taus(M:R-1));

        cfg.dims = [dimsL dimsR];
        cfg.taus = [tausL tausR];

        % Kraskov MI
        res = nta_MIkraskov(Xleft, Xright, cfg);
        TC = TC + res.MI;
    end
end

end

%% -----------------------------------------------------
function scheme = build_tree_scheme(left_idx, right_idx)
% Recursively build a balanced binary-tree partition scheme
scheme = [];
if right_idx - left_idx > 1
    mid = floor((left_idx + right_idx)/2);
    scheme = [scheme; left_idx mid right_idx];
    scheme = [scheme; build_tree_scheme(left_idx, mid)];
    scheme = [scheme; build_tree_scheme(mid, right_idx)];
end
end

