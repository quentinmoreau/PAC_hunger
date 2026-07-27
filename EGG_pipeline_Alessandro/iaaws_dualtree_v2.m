function [sgates, errors] = iaaws_dualtree_v2(data, num_surr, max_iter, channel_id)
% IAAWS_DUALTREE_V2 Iterative Amplitude Adjusted Wavelet Transform (IAAWT)
%
% Generates surrogate data preserving both the amplitude distribution and 
% the dual-tree complex wavelet transform (DTCWT) amplitudes of the 
% original time series. This is a robust, native-MATLAB implementation 
% based on the algorithm logic presented in Keylock (2017).
%
% INPUTS:
%   data     - 1D array, original time series.
%   num_surr - Integer, number of surrogates to generate (default: 1).
%   max_iter - Integer, maximum number of iterations per surrogate (default: 1000).
%
% OUTPUTS:
%   sgates   - Matrix of size [length(data), num_surr] containing the surrogates.
%   errors   - Vector of size [num_surr, 1] containing the final total error 
%              for each surrogate upon termination.
%
% REQUIREMENTS:
%   Wavelet Toolbox (MATLAB R2021a or later).

    %% --- Input Validation & Defaults ---
    if nargin < 3, max_iter = 1000; end

    %% --- Step 0: Set Random Seed for Reproducibility ---
    rng(channel_id, 'twister');
    
    if ~exist('dualtree', 'file')
        error('This function requires the MATLAB Wavelet Toolbox (R2021a+).');
    end

    % Convergence thresholds (Adjustable parameters)
    accerror = 0.01;       % Absolute acceptable error threshold
    error_change = 10;     % Relative error change threshold to prevent stalling

    data = double(data(:));
    n = length(data);

    % --- 1. Dyadic Zero-Padding ---
    % DTCWT requires the signal length to be a power of 2.
    % The signal is right-aligned in the padded array to match R logic.
    numlevels = floor(log2(n)); 
    len_pad = 2^(numlevels + 1);
    count = len_pad - n + 1;
    
    x_pad = zeros(len_pad, 1);
    x_pad(count:end) = data;
    
    % --- 2. Original Signal Decomposition ---
    % A_orig: Approximation coefficients (Low-pass)
    % D_orig: Detail coefficients (High-pass, complex cell array)
    [A_orig, D_orig] = dualtree(x_pad, 'Level', numlevels);
    
    % Extract target wavelet amplitudes
    Amp_D_orig = cellfun(@abs, D_orig, 'UniformOutput', false);
    
    % Prepare target time-domain amplitudes
    sorted_data = sort(data);
    std_val = std(sorted_data); 
    if std_val == 0, std_val = 1; end % Prevent division by zero
    
    % Pre-allocate outputs
    sgates = zeros(n, num_surr);
    errors = zeros(num_surr, 1);
    
    fprintf('Generating %d IAAWFT surrogates\n', num_surr);
    
    % --- 3. Surrogate Generation Loop ---
    for k = 1:num_surr
        
        % Initialize with a random permutation of the original data
        z = zeros(len_pad, 1);
        z(count:end) = data(randperm(n)); 
        
        % Initialize error trackers
        amperror = 100;
        waverror = 100;
        
        for iter = 1:max_iter
            old_z = z; % Store for wavelet error calculation
            
            % A. Surrogate Decomposition
            [~, D_curr] = dualtree(z, 'Level', numlevels);
            
            % B. Wavelet Domain Projection (Phase Preservation)
            D_new = cell(size(D_curr));
            for j = 1:numlevels
                % Retain surrogate's phase, impose original's amplitude
                Phase_curr = angle(D_curr{j});
                D_new{j} = Amp_D_orig{j} .* exp(1i * Phase_curr);
            end
            
            % C. Signal Reconstruction
            % Impose original approximation trend (A_orig)
            z_temp = real(idualtree(A_orig, D_new));
            
            % --- Calculate Wavelet Error ---
            wavdiff = mean(abs(z_temp - old_z));
            curr_waverror = wavdiff / std_val;
            
            z = z_temp;
            old_z_amp = z; % Store for amplitude error calculation
            
            % D. Time Domain Projection (Rank Ordering)
            valid_part = z(count:end);
            [~, idx] = sort(valid_part); 
            
            % Replace values with original amplitudes based on sorted ranks
            valid_part(idx) = sorted_data;
            
            % Rebuild padded signal
            z = zeros(len_pad, 1);
            z(count:end) = valid_part;
            
            % --- Calculate Amplitude Error ---
            ampdiff = mean(abs(z - old_z_amp));
            curr_amperror = ampdiff / std_val;
            
            % E. Convergence Criteria
            toterror = curr_amperror + curr_waverror;
            oldtoterr = amperror + waverror;
            
            % Stop if either error domain reaches the absolute tolerance
            if (curr_amperror <= accerror) || (curr_waverror <= accerror)
                break;
            end
            
            % Stop if the relative improvement stalls
            if abs((oldtoterr - toterror) / toterror) < (accerror / error_change)
                break;
            end
            
            % Update error states for the next iteration
            amperror = curr_amperror;
            waverror = curr_waverror;
        end
        
        % Store final surrogate (strip padding) and final error
        sgates(:, k) = z(count:end);
        errors(k) = toterror;
        
        % Progress tracker
        if mod(k, 10) == 0, fprintf('.'); end
    end
    fprintf('\nGeneration complete.\n');
end