function [filt_buff, filtered_signal, an_edges_filt_opt] = find_anomalies_Tuk_COF_RRCF_v4(x, wind, fs)
    % FIND_ANOMALIES_TUK_COF_RRCF
    % Detects anomalies in a signal using Tukey, COF, and RRCF methods.
    % 
    % Inputs:
    %   x   - original signal
    %   wind - window size
    %   fs  - sampling frequency
    %
    % Outputs:
    %   filt_buff          - buffered filtered signal
    %   filtered_signal    - Chebyshev-filtered signal
    %   an_edges_filt_opt  - start and end indices of final detected anomalies

    % === Step 1: Filter the signal with a Chebyshev filter ===
    filtered_signal = cheb_EGG_filt(x, fs);
    if isrow(filtered_signal)
        filtered_signal = filtered_signal';
    end

    % === Step 2: Buffer the filtered signal into overlapping windows ===
    filt_buff = buffer(filtered_signal, wind, wind-1, 'nodelay');

    % === Step 3: Anomaly detection using liberal Tukey method ===
    %---windowed (local) Tukey
    r_window = isanomaly_FE(filtered_signal, 'liberal tukey', wind);
    %---global (local) Tukey
    r_tot = isanomaly_FE(filtered_signal, 'liberal tukey', Inf);
    r_window = ensure_row(r_window);
    r_tot = ensure_row(r_tot);

    % === Step 4: COF (Connectivity Outlier Factor) anomaly detection ===
    [~, cof, ~] = COF(decimate(filtered_signal, 10), 60);
    threshold = median(cof) + 4 * mad(cof);
    r_cof = double(cof > threshold);

    % === Step 5: Merge short gaps in COF and global Tukey results ===
    r1_cof = merge_short_gaps(r_cof, 2);
    r_tot_d_m = merge_short_gaps(downsample(r_tot, 10), 2);

    % === Step 6: Combine COF and global Tukey anomalies ===
    r1_updatedCOFTuk = anomaly_intersection(r_tot_d_m', r1_cof');

    % === Step 7: Interpolate back to original signal length ===
    r_tot = interp1( ...
        linspace(1, length(filtered_signal), length(r1_updatedCOFTuk)), ...
        double(r1_updatedCOFTuk), ...
        1:length(filtered_signal), ...
        'nearest', 0);
    r_tot = logical(r_tot);

    % === Step 8: Combine windowed and global Tukey anomalies ===
    r1 = (r_tot == 1) | (r_window == 1);
    r1 = ensure_row(r1);

    % === Step 9: Extract anomaly segment boundaries ===
    an_edges_filt_tot = find(diff([0, r1(:)', 0]) ~= 0);
    an_edges_filt_tot(an_edges_filt_tot > length(x)) = length(x);

    % === Step 10: Random Cut Forest anomaly detection ===
    
    % Fissiamo il seed per garantire riproducibilità totale al 100%
    %old_rng = rng('twister'); 
    
    rng('default'); 
    [~, ~, scores] = rrcforest(filtered_signal, ...
        ContaminationFraction = sum(r1) / length(r1), ...
        CollusiveDisplacement = "maximal", ...  % Ottimizzato per artefatti contigui
        NumLearners = 1000, ...                 % Stabilizza gli score ed elimina sfarfallii
        NumObservationsPerLearner = 256 );
        
    % Ripristiniamo il generatore casuale
    %rng(old_rng);

    % === Step 11: Tukey test on RCF scores ===
    r_window_sc = isanomaly_FE(scores, 'conservative tukey', Inf);
    r_window_sc = ensure_row(r_window_sc);
    r_window_sc = merge_short_gaps(r_window_sc, 20);

    an_edges_filt_sc = find(diff([0, r_window_sc(:)', 0]) ~= 0);
    an_edges_filt_sc(an_edges_filt_sc > length(x)) = length(x);

    % === Step 12: Final anomaly mask (Tukey ∩ RCF) ===
    r1_updated = anomaly_intersection(r_window_sc, r1);
    %an_edges_filt_opt = find(diff([0, r1_updated, 0]) ~= 0);
    an_edges_filt_opt = find(diff([0, r1_updated(:)', 0]) ~= 0);
    an_edges_filt_opt(an_edges_filt_opt > length(x)) = length(x);

end

% === Helper function to ensure row vector ===
function row = ensure_row(vec)
    if size(vec,1) > size(vec,2)
        row = vec';
    else
        row = vec;
    end
end
