function [parameters, tp, ist_ef1, if_curve_ridge] = par_est_TFR_v5(x, t, plots, anomalies, fs, faslt, filled_sig, numPart)
% PAR_EST_TFR_V5
% This function estimates multiple spectral and nonlinear parameters from a time-frequency representation (TFR).
% It uses a dual-track approach:
% 1. Spectral features: Uses 'faslt' (with anomalies excluded) for unbiased statistics.
% 2. Dynamic features: Recalculates TFR on 'filled_sig' for continuous phase space embedding.

    if nargin < 8
        numPart = 1;
    end

    %% === Valid mask creation for statistics ===
    valid_mask = true(1, length(t));
    if ~isempty(anomalies)
        for i = 1:2:length(anomalies)-1
            start_idx = max(anomalies(i), 1);
            end_idx = min(anomalies(i+1), length(t));
            valid_mask(start_idx:end_idx) = false;
        end
    end
    norm_time = sum(valid_mask);
    column_mask = valid_mask;

    %% === Power spectrum and instantaneous energy ===
    freq_faslt = 0.016:0.0001:0.16;
    faslt_abs = abs(faslt)/2;
    norm_pow = sum(faslt_abs.^2, 2); % Normalize over time
    freq_faslt = freq_faslt(:);      % Ensure column vector
    
    % Instantaneous energy per unit time (normalized by frequency)
    time_pow = sum(faslt_abs.^2 ./ freq_faslt, 1);
    tp = time_pow;

    %% === Instantaneous spectral entropy ===
    ist_f1 = abs(faslt).^2 ./ sum(abs(faslt).^2, 1); % Instantaneous power distribution
    ist_ef1 = zeros(1, size(faslt, 2));             % Preallocate
    valid_columns = valid_mask;                     % Valid columns only
    for i = find(valid_columns)
        ist_ef1(i) = -sum(ist_f1(:, i) .* log2(ist_f1(:, i)), 'omitnan');
        ist_ef1(i) = ist_ef1(i) / log2(size(faslt, 1)); % Normalize entropy
    end
    if ~any(valid_columns)
        ist_ef1(:) = NaN;
    end

    %% === Dominant frequency and power ===
    peak_f = freq_faslt(norm_pow == max(norm_pow));
    peak_pow = max(norm_pow);

    %% === Spectral centroid and spread ===
    spectral_centroid = sum(freq_faslt .* norm_pow) / sum(norm_pow);
    spectral_spread = sqrt(sum(norm_pow .* (freq_faslt - spectral_centroid).^2) / sum(norm_pow));

    %% === Mean and std of instantaneous energy ===
    mean_En = mean(time_pow(valid_mask));
    std_En = std(time_pow(valid_mask));

    %% === Ridge-based normogastria fraction ===
    % Ensure row vector for inst_freq_frag_v2 to avoid dimension mismatch
    if size(freq_faslt, 2) == 1
        freq_faslt=freq_faslt';
    end
    [if_curve_ridge, perc] = inst_freq_frag_v2(faslt, freq_faslt, anomalies, t, 5e4, plots, 'off', numPart);
    
    freq_faslt = freq_faslt(:); % Re-ensure column vector for subsequent steps

    %% === Power distribution into brady, normo, tachy ranges ===
    bradi = find_closest(freq_faslt, 0.016);
    low   = find_closest(freq_faslt, 0.033);
    high  = find_closest(freq_faslt, 0.067);
    tachi = find_closest(freq_faslt, 0.16);

    npow  = trapz(norm_pow(low:high))  / length(norm_pow(low:high))  / norm_time;
    bpow  = trapz(norm_pow(bradi:low)) / length(norm_pow(bradi:low)) / norm_time;
    tpow  = trapz(norm_pow(high:tachi))/ length(norm_pow(high:tachi))/ norm_time;

    n_pow_rel = npow / (npow + bpow + tpow);
    b_pow_rel = bpow / (npow + bpow + tpow);
    t_pow_rel = tpow / (npow + bpow + tpow);

    %% === Dominant, super- and sub-dominant band power ===
    dom_band_l = find_closest(freq_faslt, peak_f - 0.015);
    dom_band_h = find_closest(freq_faslt, peak_f + 0.015);
    
    % Adjust if indices are at the boundaries
    if dom_band_l == 1
        dom_band_l = dom_band_l + 1;
    end
    if dom_band_h == length(freq_faslt)
        dom_band_h = dom_band_h - 1;
    end

    npow_band  = sum(norm_pow(dom_band_l:dom_band_h)) / length(norm_pow(dom_band_l:dom_band_h)) / norm_time;
    bpow_band  = sum(norm_pow(1:dom_band_l))         / length(norm_pow(1:dom_band_l))         / norm_time;
    tpow_band  = sum(norm_pow(dom_band_h:end))       / length(norm_pow(dom_band_h:end))       / norm_time;

    n_pow_rel_freq = npow_band / (bpow_band + tpow_band);
    b_pow_rel_freq = bpow_band / (npow_band + tpow_band);
    t_pow_rel_freq = tpow_band / (npow_band + bpow_band);

    %% === Entropies from power spectrum ===
    norm_pow_prob = norm_pow / sum(norm_pow); 
    Spectral_entropy = -sum(norm_pow_prob(norm_pow_prob > 0) .* log2(norm_pow_prob(norm_pow_prob > 0)));
    Spectral_entropy = Spectral_entropy / log2(length(norm_pow_prob)); 

    [~, PS_tsallis_entropy] = minim_Tallis(norm_pow_prob, 500, 'off');

    %% === Continuous Signal & Energy for Nonlinear Dynamics ===
    if isempty(anomalies)
        time_pow_cont = time_pow;
        filt_sig = cheb_EGG_filt(x, fs); 
    else
        % Recalculate TFR on filled_sig to get continuous energy
        [faslt_cont, ~, ~, ~] = general_faslt_plot_no_filt(filled_sig, t, fs, freq_faslt', 0.01, 0.2, [], 'off', 'off', 'zero');
        time_pow_cont = sum((abs(faslt_cont)/2).^2 ./ freq_faslt, 1);
        filt_sig = cheb_EGG_filt(filled_sig, fs);
    end

    %% === Phase space and nonlinear dynamics on continuous signal ===
    [~, eLag, eDIm] = phaseSpaceReconstruction(filt_sig, 'MaxDim', 9, 'MaxLag', 100);
    SEn_sig   = SampEn(filt_sig, 'm', eDIm, 'tau', eLag);  SEn_sig = SEn_sig(end);
    De_sig    = DispEn(filt_sig, 'm', eDIm, 'tau', eLag, 'Norm', true);

    corrDim = correlationDimension(filt_sig, eLag, eDIm, 'NumPoints', 100);
    LLE     = lyapunovExponent(filt_sig, fs, eLag, eDIm);
    
    results = nta_AIS(filt_sig, struct('dim',eDIm,'tau',eLag,'verbose',0)); 
    AIS = results.AIS;
    results = nta_dfa(filt_sig, struct('scales',[max(2,round(length(filt_sig)/2000)), round(length(filt_sig)/300)], 'verbose', 0, 'plt', 0));
    dfa_exp = results.expo;

    %% === Phase space and entropies on instantaneous energy ===
    time_pow_cont_norm = time_pow_cont / sum(time_pow_cont);
    [~, eLagTP, eDImTP] = phaseSpaceReconstruction(time_pow_cont_norm, 'MaxDim', 9, 'MaxLag', 100);
    SEn_TP = SampEn(time_pow_cont_norm, 'm', eDImTP, 'tau', eLagTP);  SEn_TP = SEn_TP(end);
    De_TP  = DispEn(time_pow_cont_norm, 'm', eDImTP, 'tau', eLagTP, 'Norm', true);

    %% === Renyi, Hoyer and Gradient entropies from valid TFR ===
    valid_TFR = faslt(:, column_mask);
    t_valid = t(column_mask);
    R = renyi(abs(valid_TFR).^2, t_valid, freq_faslt, 3) / log2(numel(valid_TFR));
    [~, ~, M_H, ~] = tfr_measures(valid_TFR);

    Gradient_entropy = GradEn_v2(abs(valid_TFR).^2);

    valid_TFR_dom_band = valid_TFR(dom_band_l:dom_band_h, :);
    G_dom = GradEn_v2(abs(valid_TFR_dom_band).^2);
    GradEn_relative = G_dom / (Gradient_entropy);

    %% === Store all results into a struct ===
    % The order below is strictly aligned with the config_db 
    % in create_spider_plot to ensure correct index mapping.
    
    parameters = struct(...
        ... % --- FREQUENCY & GLOBAL POWER (1-4) ---
        'Dominant_frequency', peak_f, ...
        'Dominant_power', peak_pow, ...
        'Spectral_centroid', spectral_centroid, ...
        'spectral_spread', spectral_spread, ...
        ...
        ... % --- ABSOLUTE POWER PER BAND (5-7) ---
        'Total_normogastric_power', npow, ...
        'Total_bradigastric_power', bpow, ...
        'Total_tachigastric_power', tpow, ...
        ...
        ... % --- RELATIVE POWER PER BAND (8-13) ---
        'Mean_relative_normogastric_power', n_pow_rel, ...
        'Mean_relative_bradigastric_power', b_pow_rel, ...
        'Mean_relative_tachigastric_power', t_pow_rel, ...
        'Mean_relative_dominant_power', n_pow_rel_freq, ...
        'Mean_relative_Super_dominant_power', t_pow_rel_freq, ...
        'Mean_relative_Sub_dominant_power', b_pow_rel_freq, ...
        ...
        ... % --- POWER SPECTRUM ENTROPIES (14-15) ---
        'Spectral_entropy', Spectral_entropy, ...
        'PS_Tsallis_entropy', PS_tsallis_entropy, ...
        ...
        ... % --- TIME-FREQUENCY (TFR) STATISTICS (16-20) ---
        'Fraction_of_normogastria', perc, ...
        'Mean_ridge_frequencies', mean(if_curve_ridge, 'omitnan'), ...
        'Ridge_frequencies_std', std(if_curve_ridge, 'omitnan'), ...
        'Mean_energy', mean_En, ...
        'Energy_std', std_En, ...
        ...
        ... % --- TIME-FREQUENCY ENTROPIES & SPARSITY (21-24) ---
        'Hoyer_measure', M_H, ...
        'Renyi_entropy', R, ...
        'Gradient_entropy', Gradient_entropy, ...
        'Relative_Gradient_TFR_entropy', GradEn_relative, ...
        ...
        ... % --- TIME DOMAIN ENTROPIES (25-28) ---
        'Sig_dispersion_entropy', De_sig, ...
        'Sig_sample_entropy', SEn_sig, ...
        'InstEnergy_dispersion_entropy', De_TP, ...
        'InstEnergy_sample_entropy', SEn_TP, ...
        ...
        ... % --- NONLINEAR METRICS & COMPLEXITY (29-32) ---
        'corrDim', corrDim, ...
        'LLE', LLE, ...
        'AIS', AIS, ...
        'DFA_exponent', dfa_exp);
end

%% helper functions 
function idx = find_closest(vec, val)
    [~, idx] = min(abs(vec - val));
end