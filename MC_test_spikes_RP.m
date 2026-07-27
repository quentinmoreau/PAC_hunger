function [parameters_results] = MC_test_spikes_RP(signal, t, fs, nSurrogates, channel_id)
% MC_TEST_SPIKES_SHUFFLING Statistical detection of spike artifacts using Tsallis Entropy.
%
% This function evaluates the presence of spike-like artifacts in a signal by 
% comparing the Tsallis Entropy of the original data against a distribution 
% generated from shuffled surrogates (Random Permutation, RP). A significantly lower entropy in the 
% original signal indicates the presence of structured, high-energy events 
% (spikes) compared to randomized noise.
%
% USAGE:
%   [parameters_results] = MC_test_spikes_shuffling(signal, t, fs, nSurrogates)

    fprintf('--- Starting artifact detection: Tsallis Entropy & Shuffling Test ---\n');
    tic;

    %% --- Step 0: Set Random Seed for Reproducibility ---
    rng(channel_id, 'twister');
    
    %% --- Step 1: Signal Preprocessing ---
    % Apply band-pass filtering and remove edge effects
    signal = cheb_EGG_filt(signal, fs);       
    signal = signal(5:end-5);                 
    t = t(5:end-5);
    
    % Downsampling for computational efficiency
    signal = decimate(signal, 10);            
    t = downsample(t, 10);
    fs = fs / 10;
    
    %% --- Step 2: Time-Frequency Representation (Original Signal) ---
    freq_faslt = 0.016:0.002:0.16;
    % Extract instantaneous power using Fractional Adaptive Superlet Transform (FASLT)
    [~, time_pow, ~] = general_faslt_plot_mute(signal, t, fs, freq_faslt, 0.016, 0.16, [], 'off', 'off', 'zero');
    
    %% --- Step 3: Tsallis Entropy Estimation ---
    [~, entr_orig] = minim_Tallis(time_pow', 500, 'off');
    
    %% --- Step 4: Surrogate Data Generation ---
    % Generate a distribution of null-hypothesis surrogates via random shuffling
    surrogates = zeros(nSurrogates, length(signal));
    for i = 1:nSurrogates
        surrogates(i, :) = signal(randperm(length(signal)));
    end
    
    %% --- Step 5: Iterative Entropy Computation for Surrogates ---
    entr_surr = zeros(1, nSurrogates);
    
    for i = 1:nSurrogates
        if mod(i, round(nSurrogates/10)) == 0 || i == 1
            fprintf('Processing surrogate %d of %d...\n', i, nSurrogates);
        end
        
        % Compute instantaneous power for each surrogate
        [~, time_pow_surr, ~] = general_faslt_plot_mute(surrogates(i, :), t, fs, freq_faslt, 0.016, 0.16, [], 'off', 'off', 'zero');
        
        % Calculate Tsallis Entropy for the surrogate distribution
        [~, entr_surr(i)] = minim_Tallis(time_pow_surr, 500, 'off');
    end
    
    %% --- Step 6: Statistical Analysis (Left-Tailed Test) ---
    % p-value defined as the probability that a surrogate has lower entropy than the original
    p_valueL = (sum(entr_surr <= entr_orig) + 1) / (nSurrogates + 1);
    
    parameters_results = struct();
    parameters_results.TsEn.original_value   = entr_orig;
    parameters_results.TsEn.p_valueL         = p_valueL;
    parameters_results.TsEn.mean_surrogates  = mean(entr_surr);
    parameters_results.TsEn.std_surrogates   = std(entr_surr);
    
    %% --- Step 7: Results Visualization ---
    fprintf('Original Tsallis Entropy: %.4f\n', entr_orig);
    fprintf('Statistical p-value: %.4f\n', p_valueL);
    
    artifact_detected = p_valueL < 0.05;
    if artifact_detected
        fprintf('>> STATUS: Artifact detected (Low Entropy).\n');
    else
        fprintf('>> STATUS: Signal clean (No significant spikes).\n');
    end
    
    % Visualization of the Null Hypothesis Distribution
    figure; hold on;
    histogram(entr_surr, 'FaceColor', [0.4 0.7 1], 'EdgeColor', 'none', 'FaceAlpha', 0.9);
    xline(entr_orig, 'r-', 'LineWidth', 2);
    
   title(sprintf('Artifact Detection via RP Surrogates for channel %d', channel_id), 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Tsallis Entropy Value');
    ylabel('Frequency');
    
    if artifact_detected
        subtitle(sprintf('p = %.3f (Significant Artifacts)', p_valueL), 'Color', [0.8 0 0], 'FontWeight', 'bold');
    else
        subtitle(sprintf('p = %.3f (Clean Signal)', p_valueL), 'Color', [0 0.6 0], 'FontWeight', 'bold');
    end
    
    grid on; 
    ax = gca; 
    ax.GridAlpha = 0.3;
    legend({'Shuffled Surrogates', 'Original Signal'}, 'Location', 'northeast');
    hold off;
    
    toc;
end