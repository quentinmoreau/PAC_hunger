function [index, x_new] = find_nyquist_spike_v5(x, t, fs, delay)
% FIND_NYQUIST_SPIKE_V4
% Detects and corrects spikes or jumps in a signal.
% Handles:
%   - Persistent jumps: signal shifts and stays at a new level -> corrected via offset shift
%   - Short transient spikes: flattened locally
%
% INPUTS:
%   x     - input signal
%   t     - time vector
%   fs    - sampling frequency
%   delay - window radius for local mean/std estimation
%
% OUTPUTS:
%   index - indices of detected spikes/jumps
%   x_new - corrected signal

%% Step 1: Amplitude threshold based on Nyquist power contribution
lim = pi * (fs)/8;   % empirical factor

%% Step 2: First-order difference
d = diff(x);

%% Step 3: Detect large changes
index = find(abs(d) >= lim);

%% Step 4: Initialize corrected signal
x_new = x;

%% Step 5: Correct spikes/jumps
if ~isempty(index)

    % Plot original signal with detected spike locations
    figure;
    subplot(2,1,1)
    plot(t, x)
    title('Original signal')
    hold on
    for i = 1:length(index)
        xline(t(index(i)), '--r');
    end
    ylabel('\muV')

    % Parameters for robust statistics
    k = 2;                        % tolerance in number of standard deviations
    min_duration_for_jump = 10;   % min number of samples for persistent jump
    min_pts = max(5, floor(delay/4));  % minimum samples required to estimate stats

    for i = 1:length(index)

        %% ----------------------------------------------------------
        %  (A) ROBUST ESTIMATION OF LOCAL MEAN AND STD (BEFORE the spike)
        % ----------------------------------------------------------
        start_ind = max(1, index(i) - delay);
        left_window = x(start_ind : index(i)-1);

        if length(left_window) >= min_pts
            % Case 1: left window has enough samples
            local_mean = mean(left_window);
            local_std  = std(left_window);

        else
            % Case 2: use the right window instead
            end_ind = min(length(x), index(i) + delay);
            right_window = x(index(i)+1 : end_ind);

            if length(right_window) >= min_pts
                % Right window has enough samples
                local_mean = mean(right_window);
                local_std  = std(right_window);

            else
                % Case 3: fallback to robust baseline estimation
                N = min(30, length(x));
                baseline = x(1:N);

                local_mean = median(baseline);
                local_std  = 1.4826 * median(abs(baseline - local_mean)); % MAD-based std estimate

                % Additional fallback if variance is zero
                if local_std == 0
                    local_std = std(baseline);
                    if local_std == 0
                        local_std = eps;   % avoid division-by-zero
                    end
                end
            end
        end

        %% ----------------------------------------------------------
        %  (B) DETERMINE DURATION OF THE DEVIATION
        % ----------------------------------------------------------
        j = index(i) + 1;
        while j <= length(x) && abs(x(j) - local_mean) > k * local_std
            j = j + 1;
        end
        duration = j - index(i);

        %% ----------------------------------------------------------
        %  (C) APPLY CORRECTION BASED ON DURATION
        % ----------------------------------------------------------
        if duration >= min_duration_for_jump
            % Persistent jump: estimate level AFTER the jump
            end_ind_after = min(length(x), index(i) + delay);
            after_window = x(index(i)+1 : end_ind_after);

            if ~isempty(after_window)
                level_after = median(after_window);
            else
                level_after = x(index(i)+1);
            end

            % Offset introduced by the jump
            offset = level_after - local_mean;

            % Shift the rest of the signal to remove the offset,
            % preserving its original shape/variability
            x_new(index(i)+1:end) = x(index(i)+1:end) - offset;

            % Set the jump sample itself to match the corrected level
            if index(i)+1 <= length(x_new)
                x_new(index(i)) = x_new(index(i)+1);
            else
                x_new(index(i)) = local_mean;
            end

        else
            % Short spike: flatten only the affected samples
            x_new(index(i):j-1) = local_mean;
        end

    end

    % Plot corrected signal
    subplot(2,1,2)
    plot(t, x_new)
    title('Adjusted signal')
    xlabel('Time (s)')
    ylabel('\muV')
end

end