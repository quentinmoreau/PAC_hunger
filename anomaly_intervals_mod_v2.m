function [anomaly_intervals_new] = anomaly_intervals_mod_v2(anomaly_intervals, t, factor, fs)
%ANOMALY_INTERVALS_MOD Adjusts anomaly intervals by removing short gaps between them
%   anomaly_intervals: vector of start and end indices of anomalies (e.g., [s1 e1 s2 e2 ...])
%   t: time vector
%   factor: multiplier for time threshold (in minutes)
%   fs: sampling frequency
%   Returns: modified anomaly_intervals with short gaps merged and ends adjusted

    if isempty(anomaly_intervals)
        anomaly_intervals_new = anomaly_intervals;
        return;
    end

    % Initialize index list for intervals to remove (short gaps)
    indices_to_remove = [];

    % Threshold duration in samples
    threshold = factor * fs;

    % Loop through each pair of anomaly intervals (between end of one and start of next)
    for i = 2:2:length(anomaly_intervals)-1
        gap_length = anomaly_intervals(i+1) - anomaly_intervals(i);
        if gap_length < threshold
            indices_to_remove = [indices_to_remove i i+1]; % Mark gap for removal
        end
    end

    % If the final anomaly ends very close to the end of the signal, extend it
    if anomaly_intervals(end) ~= length(t) && (length(t) - anomaly_intervals(end) < threshold)
        anomaly_intervals(end) = length(t);
    end

    % Remove intervals corresponding to short gaps
    anomaly_intervals(indices_to_remove) = [];

    % If the first anomaly starts very close to the beginning, shift it to the start
    if anomaly_intervals(1) ~= 1 && (anomaly_intervals(1) - 1 < threshold)
        anomaly_intervals(1) = 1;
    end

    anomaly_intervals_new = anomaly_intervals;
end
