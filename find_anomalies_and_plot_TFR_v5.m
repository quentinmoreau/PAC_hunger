function [filled_x, an_edges_filt_opt, faslt, time_pow, norm_pow] = find_anomalies_and_plot_TFR_v5(x, t, wind, fs, w_ext, channel_id)
% FIND_ANOMALIES_AND_PLOT_TFR
% Detects anomalies in a filtered signal, fills them, and plots the time-frequency representation.
%
% Inputs:
%   x      - original signal
%   t      - time vector
%   wind   - sliding window size
%   fs     - sampling frequency
%   w_ext  - number of samples to extend anomaly boundaries
%   channel_id - EGG channel ID number
%
% Outputs:
%   filled_x         - signal with anomalies filled
%   an_edges_filt_opt- start and end indices of extended anomalies
%   faslt            - frequency axis for time-frequency representation
%   time_pow         - power over time
%   norm_pow         - normalized power


    % === Step 1: Setting default channel number ===
    % If the function is called without a channel ID, assign a default value (1).
    % This allows standalone use without explicitly specifying the channel.
    if nargin < 6
        channel_id = 1;
    end

    % === Step 2: Frequency range for TFR ===
    freq_faslt = 0.016:0.0001:0.16;

    % === Step 3: Bandpass filter the signal ===
    filtered_signal = cheb_EGG_filt(x, fs);

    % === Step 4: Detect anomalies ===
    [~, ~, an_edges_filt_opt] = find_anomalies_Tuk_COF_RRCF_v4(x, wind, fs);

    % === Step 5: Extend anomaly edges by w_ext samples ===
    extended_anom = [];
    for i = 1:2:length(an_edges_filt_opt)-1
        extended_anom = [extended_anom, ...
                         max(1, an_edges_filt_opt(i) - w_ext), ...
                         min(an_edges_filt_opt(i+1) + w_ext, length(x))];
    end
    an_edges_filt_opt = extended_anom;

    % === Step 6: Clean short anomalies and merge if needed ===
    an_edges_filt_opt = anomaly_intervals_mod_v2(an_edges_filt_opt, t, 3, fs);  % at least 3 sec

    % === Step 7: Replace anomalies with NaN ===
    x_nan = filtered_signal;
    for i = 1:2:length(an_edges_filt_opt)-1
        x_nan(an_edges_filt_opt(i):an_edges_filt_opt(i+1)) = NaN;
    end

    % === Step 8: Fill gaps (e.g., via interpolation) ===
    filled_x = fillgaps(x_nan, 900);

    % === Step 9: Plot filtered signal with highlighted anomalies ===
    figure;
    margin_factor = 0.05;  % Extend Y-axis by 5%
    plot(t, filtered_signal, 'k'); hold on;

    y_limits = ylim;
    range_y = y_limits(2) - y_limits(1);
    extended_y_limits = [y_limits(1) - margin_factor * range_y, ...
                         y_limits(2) + margin_factor * range_y];

    % Highlight anomaly regions
    if ~isempty(an_edges_filt_opt)
        for i = 1:2:length(an_edges_filt_opt)-1
            fill([t(an_edges_filt_opt(i)), t(an_edges_filt_opt(i)), ...
                  t(an_edges_filt_opt(i+1)), t(an_edges_filt_opt(i+1))], ...
                 [extended_y_limits(1), extended_y_limits(2), ...
                  extended_y_limits(2), extended_y_limits(1)], ...
                 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
    end

    hold off;
    xlim([t(1), t(end)]);
    ylim(extended_y_limits);
    title(sprintf('Filtered channel %d signal with anomalies', channel_id));
    xlabel('Time (seconds)');
    ylabel('\muV');

    % === Step 10: Compute Time-Frequency Representation ===
    [faslt, time_pow, norm_pow] = clean_faslt_plot_no_filt_v3( ...
        filled_x, t, fs, freq_faslt, an_edges_filt_opt, 'off', 0, 'on');

end