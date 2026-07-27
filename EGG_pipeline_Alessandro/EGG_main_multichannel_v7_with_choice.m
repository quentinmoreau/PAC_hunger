%% ----------------------------------------------------------------------%%
%                   EGG Signal Processing Pipeline                        %
%%-----------------------------------------------------------------------%%

%% 0. PIPELINE SELECTION
% Create a custom modal dialog for better layout and font size control
fig = uifigure('Name', 'Pipeline Selection', ...
               'Position', [100, 100, 600, 320], ...
               'WindowStyle', 'modal');
movegui(fig, 'center');
% Initialize UserData (the script will pause until this changes)
fig.UserData = 'WAITING';
% Define the message with cleaner formatting and manual indentation
dialog_msg = {
    'Select the artifact processing mode for the GSPT:',
    '',
    ' [1] Conditional (RP Control)',
    '      Runs a Monte Carlo test on instantaneous energy. Artifact localization',
    '      only triggers if Tsallis entropy is lower than surrogates.',
    '',
    ' [2] Forced (Always Localize)',
    '      Skips the Monte Carlo check and ALWAYS runs the artifact',
    '      localization and removal algorithm.'
};
% Add the text label
uilabel(fig, 'Text', dialog_msg, ...
             'Position', [30, 80, 540, 210], ...
             'FontSize', 14, ...
             'FontName', 'Helvetica', ...
             'VerticalAlignment', 'top');
% Add the buttons
uibutton(fig, 'Position', [40, 25, 240, 40], ...
              'Text', 'Conditional (RP Control)', ...
              'FontSize', 13, ...
              'FontWeight', 'bold', ...
              'ButtonPushedFcn', @(src, event) set(fig, 'UserData', 'Conditional (RP Control)'));
uibutton(fig, 'Position', [320, 25, 240, 40], ...
              'Text', 'Forced (Always Localize)', ...
              'FontSize', 13, ...
              'FontWeight', 'bold', ...
              'ButtonPushedFcn', @(src, event) set(fig, 'UserData', 'Forced (Always Localize)'));
% Pause the script until a button is clicked OR the window is closed
waitfor(fig, 'UserData');
% Retrieve the user's choice
if isvalid(fig)
    % User clicked one of the buttons
    choice = fig.UserData;
    delete(fig); % Close the window
else
    % User closed the window using the 'X' button
    choice = ''; 
end
% Handle the case where the user cancels the selection
if isempty(choice) || strcmp(choice, 'WAITING')
    error('Selection canceled. The script has been aborted.');
end
% Assign the appropriate function handle based on the user's choice
if strcmp(choice, 'Conditional (RP Control)')
    pipeline_func = @pipeline_single_channel_v12;
    disp('Selected: pipeline_single_channel_v11 (Conditional RP Control)');
else
    pipeline_func = @pipeline_single_channel_v12_no_RP_control;
    disp('Selected: pipeline_single_channel_v11_no_RP_control (Forced Localization)');
end

%% 0.5 SD THRESHOLDS CONFIGURATION
fig_sd = uifigure('Name', 'SD Threshold Configuration', ...
               'Position', [100, 100, 660, 420], ...
               'WindowStyle', 'modal');
movegui(fig_sd, 'center');
fig_sd.UserData = 'WAITING';

dialog_msg_sd = {
    'Do you want to run the SD selection algorithm using the',
    'default thresholds?',
    '',
    'These values represent the maximum allowable difference (range) in Hz',
    'between the Dominant Frequencies (DF) across channels.',
    '',
    ' Default Thresholds:',
    '  • Strict:  0.005 Hz',
    '  • Medium:  0.010 Hz',
    '  • Loose:   0.020 Hz',
    '',
    ' ⚠️ NOTE: Modifying these values does NOT alter the core EGG feature',
    ' extraction. It only updates the channel reliability labels (Good/Bad) in',
    ' the final Excel outputs and modifies the channel subsets used to',
    ' compute the Total Correlations.',
    '',
    'Select an option below:'
};

uilabel(fig_sd, 'Text', dialog_msg_sd, ...
             'Position', [30, 80, 600, 320], ...
             'FontSize', 14, ...
             'FontName', 'Helvetica', ...
             'VerticalAlignment', 'top');

uibutton(fig_sd, 'Position', [60, 25, 240, 40], ...
              'Text', 'Use Default Thresholds', ...
              'FontSize', 13, ...
              'FontWeight', 'bold', ...
              'ButtonPushedFcn', @(src, event) set(fig_sd, 'UserData', 'DEFAULTS'));

uibutton(fig_sd, 'Position', [360, 25, 240, 40], ...
              'Text', 'Set Custom Thresholds', ...
              'FontSize', 13, ...
              'FontWeight', 'bold', ...
              'ButtonPushedFcn', @(src, event) set(fig_sd, 'UserData', 'CUSTOM'));

waitfor(fig_sd, 'UserData');

if isvalid(fig_sd)
    sd_choice = fig_sd.UserData;
    delete(fig_sd); 
else
    sd_choice = ''; 
end

if isempty(sd_choice) || strcmp(sd_choice, 'WAITING')
    error('Selection canceled. The script has been aborted.');
end

% Initialize default thresholds in Hz
sd_strict = 0.005;
sd_medium = 0.010;
sd_loose  = 0.020;

if strcmp(sd_choice, 'CUSTOM')
    fig_custom = uifigure('Name', 'Custom SD Thresholds', ...
                   'Position', [100, 100, 420, 300], ...
                   'WindowStyle', 'modal');
    movegui(fig_custom, 'center');
    fig_custom.UserData = 'WAITING';
    
    uilabel(fig_custom, 'Text', 'Enter custom SD thresholds (in Hz):', ...
                 'Position', [40, 230, 340, 30], ...
                 'FontSize', 14, 'FontName', 'Helvetica', 'FontWeight', 'bold');
    
    uilabel(fig_custom, 'Text', 'Strict (Hz):', 'Position', [60, 170, 100, 30], 'FontSize', 13, 'FontName', 'Helvetica');
    ef_strict = uieditfield(fig_custom, 'numeric', 'Position', [160, 170, 160, 30], ...
                            'Value', sd_strict, 'FontSize', 13);
    
    uilabel(fig_custom, 'Text', 'Medium (Hz):', 'Position', [60, 120, 100, 30], 'FontSize', 13, 'FontName', 'Helvetica');
    ef_medium = uieditfield(fig_custom, 'numeric', 'Position', [160, 120, 160, 30], ...
                            'Value', sd_medium, 'FontSize', 13);
    
    uilabel(fig_custom, 'Text', 'Loose (Hz):', 'Position', [60, 70, 100, 30], 'FontSize', 13, 'FontName', 'Helvetica');
    ef_loose = uieditfield(fig_custom, 'numeric', 'Position', [160, 70, 160, 30], ...
                           'Value', sd_loose, 'FontSize', 13);
    
    uibutton(fig_custom, 'Position', [140, 15, 140, 35], ...
                  'Text', 'Confirm', ...
                  'FontSize', 13, ...
                  'FontWeight', 'bold', ...
                  'ButtonPushedFcn', @(src, event) set(fig_custom, 'UserData', 'CONFIRM'));
                  
    waitfor(fig_custom, 'UserData');
    
    if isvalid(fig_custom)
        sd_strict = ef_strict.Value;
        sd_medium = ef_medium.Value;
        sd_loose  = ef_loose.Value;
        delete(fig_custom);
    else
        error('Custom threshold configuration canceled. The script has been aborted.');
    end
end
fprintf('Thresholds in use -> Strict: %.3f Hz | Medium: %.3f Hz | Loose: %.3f Hz\n', sd_strict, sd_medium, sd_loose);

%% 1. SETUP
% Output folder
% ⚠️ The pipeline will generate three Excel files with the results. 
% Please change the folder path as desired.
output_folder = 'C:\Users\aless\Documents\MATLAB\Prova_stampa_pipeline';
% ⚠️ Update the path below to match your local FieldTrip installation.
addpath("C:\Users\aless\Documents\fieldtrip-20251023\fieldtrip-20251023")
ft_defaults;
% ⚠️ Update the path below to match your local data path
cfg = [];
cfg.dataset='C:\Users\aless\Documents\Paper_EGG_drafts\dati_per_github\Example_Clean.eeg';

% ⚠️ Channel selection can be done in several ways.
% By default, all available channels are selected.
% Choose the option that best suits your data.
%----------------------------------------------
% manual channel selection
cfg.channel = {'egg1', 'egg2', 'egg3', 'egg4'}; 
%-----------------------------------------------
%select all channels
%hdr = ft_read_header(cfg.dataset);
%all_channels = hdr.label;
%cfg.channel = all_channels; 
%----------------------------------------------
% Select all channels whose names contain a common identifier (e.g., 'egg')
%egg_channels = all_channels(contains(lower(all_channels), 'egg'));
%cfg.channel = egg_channels; 
%----------------------------------------------
% Extract channel names for subsequent use during data export
channel_names = cfg.channel;
% load the datas
EGG_raw = ft_preprocessing(cfg);
nCh = numel(cfg.channel); % number of channels

%% 2. REMOVE HIGH-FREQUENCY SPIKES
disp('Removing high-frequency spikes...')
for ch = 1:nCh
    [index, cleaned_signal] = find_nyquist_spike_v5(EGG_raw.trial{1}(ch,:), EGG_raw.time{1}, 1000, 100);
    EGG_raw.trial{1}(ch,:) = cleaned_signal;
    if ~isempty(index)
        sgtitle(['Channel ' num2str(ch)]);
    end
end

%% 3. RESAMPLE TO 10 Hz
disp('Resampling...')
cfg = [];
cfg.detrend = 'no';
cfg.demean  = 'yes';
cfg.resamplefs = 10;
EGG_downsampled = ft_resampledata(cfg, EGG_raw);
fs = cfg.resamplefs;

%% 4. REMOVE EDGE EFFECTS
disp('Removing filter edge effects...')
t = EGG_downsampled.time{1}(4:end-4);
x = EGG_downsampled.trial{1}(:, 4:end-4); % matrix [nCh x nTime]

%% 5. RUN PIPELINE ON EACH CHANNEL
disp('Running pipeline on individual channels...')
freq_range = 0.016:0.0001:0.16;
% Preallocation
norm_pow = cell(1, nCh);
parameters = cell(1, nCh);
efreq = cell(1, nCh);
eamp = cell(1, nCh);
ephi = cell(1, nCh);
r = cell(1, nCh);
tp = cell(1, nCh);
ief = cell(1, nCh);
amp_curve = cell(1, nCh);
ph_curve = cell(1, nCh);
if_curve = cell(1, nCh);
fill_s = cell(1, nCh);
ridge_curve = cell(1, nCh);
parameters_results = cell(1, nCh);
for ch = 1:nCh
    fprintf('Processing channel %d...\n', ch);
    % Use the selected function handle assigned in step 0
    [~, norm_pow{ch}, parameters{ch}, efreq{ch}, eamp{ch}, ephi{ch}, ...
     r{ch}, tp{ch}, ief{ch}, amp_curve{ch}, ph_curve{ch}, if_curve{ch}, ...
     fill_s{ch}, ridge_curve{ch}, parameters_results{ch}] = ...
     pipeline_func(x(ch,:), t, fs, freq_range, ch); 
end

%% 6. PLOT POWER SPECTRA
figure; hold on;
colors = lines(nCh); % Automatically generates a distinct color map
spectra = norm_pow;
labels = arrayfun(@(ch) sprintf('Channel %d', ch), 1:nCh, 'UniformOutput', false);
% --- Initial (all channels) sd and max range ---
dom_freqs = zeros(1, nCh);
for ch = 1:nCh
    [~, peak_idx] = max(spectra{ch});
    dom_freqs(ch) = freq_range(peak_idx);
end
initial_sd = std(dom_freqs);
max_range=range(dom_freqs);
% -------------------------------------------------
for ch = 1:nCh
    plot(freq_range, spectra{ch}, 'Color', colors(ch,:), 'LineWidth', 2);
end
drawnow;
yl = ylim; y_bottom = yl(1); y_top = yl(2);
for ch = 1:nCh
    [peak_val, peak_idx] = max(spectra{ch});
    dom_freq = freq_range(peak_idx);
    plot(dom_freq, peak_val, 'o', 'Color', colors(ch,:), ...
        'MarkerFaceColor', colors(ch,:), 'MarkerSize', 6, 'HandleVisibility', 'off');
    y_offset = 0.02 * y_top;
    x_shift = 0.001 * (-1)^ch;
    text(dom_freq + x_shift, peak_val + y_offset, sprintf('%.4f Hz', dom_freq), ...
        'Color', colors(ch,:), 'FontSize', 12, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'HandleVisibility', 'off');
end
% Highlights the normogastric band
[~, inf_norm] = min(abs(freq_range - 0.033));
[~, sup_norm] = min(abs(freq_range - 0.067));
fill([freq_range(inf_norm) freq_range(sup_norm) freq_range(sup_norm) freq_range(inf_norm)], ...
     [y_bottom y_bottom y_top y_top], [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
xline(freq_range(inf_norm), '--k', 'LineWidth', 1.5);
xline(freq_range(sup_norm), '--k', 'LineWidth', 1.5);
xlabel('Frequency (Hz)', 'FontSize', 14);
ylabel('Power', 'FontSize', 14);
title('Power Spectrum Density Across Channels', 'FontSize', 16);
% --- DF variability info ---
subtitle(sprintf('Dominant Frequency Variability (All Channels) | SD: %g Hz | \\Delta_{max} (Range): %g Hz', initial_sd, max_range), 'FontSize', 13);
legend(labels, 'Location', 'northeast', 'FontSize', 12);
grid on; set(gca, 'FontSize', 12);

%% 7. RIDGE CURVE CORRELATIONS
disp('Computing ridge curve correlations...');
ridge_mat = cell2mat(cellfun(@(x) x(:), ridge_curve, 'UniformOutput', false));
ridge_cor = corr(ridge_mat, 'Rows','complete');
figure;
h = heatmap(ridge_cor);
h.XDisplayLabels = labels;
h.YDisplayLabels = labels;
title('Pairwise correlations between TFR ridge curves');

%% 9. TOTAL CORRELATION CALCULATIONS AND SUBSET SELECTIONS
% ====================================================================
% A. GENERATE MASKS FOR SUBSET SELECTIONS (MC & Greedy SD Dropping)
% ====================================================================
% 1. Monte Carlo Selection Mask (Simes Test)
all_reliable = false(1, nCh);
for ch = 1:nCh
    try
        [is_globally_significant, ~] = calculate_simes_test([parameters_results{ch}.spec_skewness.p_value_right, ...
            parameters_results{ch}.spec_sparsity.p_value_right], 0.05);
        all_reliable(ch) = is_globally_significant;
    catch
        all_reliable(ch) = false;
    end
end
% 2. SD Greedy Dropping Masks
if ~exist('dom_freqs', 'var')
    dom_freqs = zeros(1, nCh);
    for ch = 1:nCh
        [~, peak_idx] = max(norm_pow{ch});
        dom_freqs(ch) = freq_range(peak_idx);
    end
end

% user defined custom values for thresholds
keep_strict = evaluate_channels_sd(dom_freqs, sd_strict);
keep_medium = evaluate_channels_sd(dom_freqs, sd_medium);
keep_loose  = evaluate_channels_sd(dom_freqs, sd_loose);

% Define Selection Scenarios
tc_scenarios = {
    'All_Channels',          true(1, nCh);
    'Monte_Carlo_Selection', all_reliable;
    sprintf('SD_Selection_Strict_%.3fHz', sd_strict), keep_strict;
    sprintf('SD_Selection_Medium_%.3fHz', sd_medium), keep_medium;
    sprintf('SD_Selection_Loose_%.3fHz', sd_loose),   keep_loose
};

n_scenarios = size(tc_scenarios, 1);
tc_results = cell(n_scenarios, 1);
% ====================================================================
% B. COMPUTE MULTIVARIATE METRICS ACROSS SUBSETS
% ====================================================================
if nCh > 1
    disp('Computing total correlations and entropies across subsets...');
    len_sig = length(fill_s{1});                 % reference signal length

    for sc = 1:n_scenarios
        scenario_name = tc_scenarios{sc, 1};
        mask          = tc_scenarios{sc, 2};
        n_sel         = sum(mask);
        sel_names     = strjoin(channel_names(mask), ', ');

        % If fewer than 2 channels are selected, network metrics are undefined
        if n_sel < 2
            tc_results{sc} = struct('Method', scenario_name, 'N_Channels', n_sel, 'Names', sel_names, ...
                'r_tot_frac', NaN, 'Metrics', NaN(1, 9));
            continue;
        end

        % Subset data arrays
        sub_fill_s = fill_s(mask);   sub_r     = r(mask);       sub_tp    = tp(mask);
        sub_ief    = ief(mask);      sub_amp   = amp_curve(mask);
        sub_if     = if_curve(mask); sub_ph    = ph_curve(mask); sub_ridge = ridge_curve(mask);

        % Build reliability matrix for the subset
        sub_r_mat = false(n_sel, len_sig);
        for ch = 1:n_sel
            ri = sub_r{ch}(:)';
            if length(ri) >= len_sig, sub_r_mat(ch, :) = ri(1:len_sig);
            else,                     sub_r_mat(ch, 1:length(ri)) = ri; end
        end
        sub_r_tot = any(sub_r_mat, 1);
        if     length(sub_r_tot) > len_sig, sub_r_tot = sub_r_tot(1:len_sig);
        elseif length(sub_r_tot) < len_sig, sub_r_tot = [sub_r_tot, false(1, len_sig - length(sub_r_tot))]; end
        sub_r_tot_frac = mean(sub_r_tot);

        if any(~sub_r_tot)
            decim        = @(x) conditional_decimate(x, sub_r_tot);
            build_matrix = @(C) cell2mat(cellfun(@(v) decim(v(:)), C, 'UniformOutput', false));

            % Effective sample count after reliability-based decimation.
            % Varies per subject/scenario -> used to normalize TC below.
            N_eff = numel(decim(sub_fill_s{1}(:)));

            % Theoretical upper bound of matrix-based Renyi-2 total correlation:
            %   0 <= TC <= (n_sel - 1) * log2(N_eff)
            % Dividing by it removes the channel-count and sample-size scale,
            % giving a normalized integration index in [0,1] comparable across subjects.
            TC_norm_factor = (n_sel - 1) * log2(N_eff);

            % --- Total correlations (raw) ---
            try TC_ch    = total_corr_general_v3(build_matrix(sub_fill_s), 2); catch, TC_ch    = NaN; end
            try TC_en    = total_corr_general_v3(build_matrix(sub_tp), 2);     catch, TC_en    = NaN; end
            try TC_ief   = total_corr_general_v3(build_matrix(sub_ief), 2);    catch, TC_ief   = NaN; end
            try TC_eamp  = total_corr_general_v3(build_matrix(sub_amp), 2);    catch, TC_eamp  = NaN; end
            try TC_efreq = total_corr_general_v3(build_matrix(sub_if), 2);     catch, TC_efreq = NaN; end
            try TC_ephi  = total_corr_general_v3(build_matrix(sub_ph), 2);     catch, TC_ephi  = NaN; end
            try TC_ridge = total_corr_general_v3(build_matrix(sub_ridge), 2);  catch, TC_ridge = NaN; end

            % --- Normalize TC to [0,1] ---
            TC_vec  = [TC_ch, TC_en, TC_ief, TC_eamp, TC_efreq, TC_ephi, TC_ridge];
            TC_norm = TC_vec / TC_norm_factor;

            % Phase space reconstruction parameters per selected channel
            eLag = zeros(1, n_sel); eDim = zeros(1, n_sel);
            for ch = 1:n_sel
                try
                    sig = sub_fill_s{ch}(:); sig = sig(1:len_sig); sig_valid = sig(~sub_r_tot);
                    if numel(sig_valid) >= 20
                        [~, eLag(ch), eDim(ch)] = phaseSpaceReconstruction(sig_valid, 'MaxDim', 9, 'MaxLag', 100);
                    else
                        eLag(ch) = 1; eDim(ch) = 2;
                    end
                catch
                    eLag(ch) = 1; eDim(ch) = 2;
                end
            end

            % Multivariate data matrix
            data_matrix = cell2mat(cellfun(@(v) v(:), sub_fill_s, 'UniformOutput', false));
            min_len     = min(cellfun(@length, sub_fill_s));
            data_matrix = data_matrix(1:min_len, :);

            % --- Multivariate entropies (NORMALIZED versions) ---
            % Embedding (eDim/eLag) is data-driven and changes per subject, so the raw
            % values are NOT comparable. MvPermEn returns the normalized entropy as its
            % 2nd output; MvDispEn is normalized via 'Norm', true.
            try [~, MPerm] = MvPermEn(data_matrix, 'm', eDim, 'tau', eLag, 'Norm', true); catch, MPerm = NaN; end
            try     MDisp  = MvDispEn(data_matrix, 'm', eDim, 'tau', eLag, 'Norm', true); catch, MDisp = NaN; end

            metrics = [TC_norm, MPerm, MDisp];   % all in [0,1]
        else
            metrics = NaN(1, 9);
        end

        tc_results{sc} = struct('Method', scenario_name, 'N_Channels', n_sel, 'Names', sel_names, ...
            'r_tot_frac', sub_r_tot_frac, 'Metrics', metrics);
    end

    % Extract Global Baseline for the Spider Plot (Scenario 1)
    tot_cor_parf = tc_results{1}.Metrics;
    TC_ch = tot_cor_parf(1); TC_en = tot_cor_parf(2); TC_ief = tot_cor_parf(3);
    TC_eamp = tot_cor_parf(4); TC_efreq = tot_cor_parf(5); TC_ephi = tot_cor_parf(6);
    TC_ridge = tot_cor_parf(7); MPerm = tot_cor_parf(8); MDisp = tot_cor_parf(9);

    % Compute global r_tot purely for the spider plot subtitle logic
    r_mat_all = false(nCh, len_sig);
    for ch = 1:nCh
        ri = r{ch}(:)';
        if length(ri) >= len_sig, r_mat_all(ch, :) = ri(1:len_sig);
        else,                     r_mat_all(ch, 1:length(ri)) = ri; end
    end
    r_tot = any(r_mat_all, 1);
    if     length(r_tot) > len_sig, r_tot = r_tot(1:len_sig);
    elseif length(r_tot) < len_sig, r_tot = [r_tot, false(1, len_sig - length(r_tot))]; end

    %% 10. PLOT FINAL RESULTS (SPIDER PLOT)
    % Good-quality criterion (Uses the MC Simes test precalculated above)
    good_quality = all(all_reliable) && mean(r_tot) < 0.6;

    % --- Color Settings  ---
    if good_quality
        theme_color = [0 0.5 0]; data_color = [0 0.6 0];
        subtitle_str = 'Signal Quality: Reliable (Green)';
    else
        theme_color = [0.7 0 0]; data_color = [0.8 0 0];
        subtitle_str = 'Signal Quality: Unreliable/Noisy* (Red)';
    end

    % All 9 metrics are now normalized to [0,1]
    precision_linear  = repmat(3, 1, 7);
    precision_entropy = repmat(3, 1, 2);

    labels_spider = {'Raw filtered data TC', 'Time energy TC', 'Time spectral entropy TC', ...
        'IA TC', 'IF TC', 'IP TC', 'TFR ridge TC', ...
        'Mv Permutation entropy', 'Mv Dispersion entropy'};
    num_params_tc = length(labels_spider);
    dummy_labels  = repmat({''}, 1, num_params_tc);
    edge_col      = repmat('w', 1, num_params_tc);
    lbl_offset    = 0.25;

    f = figure('Color', 'w');
    f.Position(3:4) = [950, 750];

    % --- Fixed limits: everything is normalized in [0,1] ---
    my_scaling = repmat({'linear'}, 1, num_params_tc);
    min_limits = zeros(1, num_params_tc);
    max_limits = ones(1, num_params_tc);

    spider_plot(tot_cor_parf, ...
        'AxesLabels', dummy_labels, ...
        'AxesLimits', [min_limits; max_limits], ...
        'AxesScaling', my_scaling, ...
        'AxesPrecision', [precision_linear, precision_entropy], ...
        'AxesLabelsEdge', edge_col, ...
        'AxesLabelsOffset', lbl_offset, ...
        'AxesFontSize', 8, ...
        'AxesFontColor', [0.2 0.2 0.2], ...
        'AxesDisplay', 'none', ...
        'AxesInterval', 5, ...
        'FillOption', {'on'}, ...
        'FillTransparency', 0.15, ...
        'Color', data_color, ...
        'AxesColor', [0.7 0.7 0.7], ...
        'LineStyle', {'-'}, ...
        'LineWidth', 2, ...
        'Marker', {'o'}, ...
        'MarkerSize', 5);

    title(sprintf('Total correlation across %d channels', nCh), 'FontSize', 14);
    full_subtitle = sprintf('%s | Fraction of excluded signal: %.2f', subtitle_str, mean(r_tot));
    subtitle(full_subtitle, 'Color', theme_color, 'FontWeight', 'bold');

    hold on;
    theta = linspace(pi/2, pi/2 - 2*pi, num_params_tc + 1);
    theta(end) = [];
    r_start = 1.02;
    r_end   = 1 + lbl_offset - 0.02;
    for i = 1:num_params_tc
        [x_guide, y_guide] = pol2cart(theta(i), [r_start, r_end]);
        plot(x_guide, y_guide, ':', 'Color', theme_color, 'LineWidth', 0.8);
        [x_txt, y_txt] = pol2cart(theta(i), r_end + 0.05);
        if x_txt > 0.1, h_align = 'left'; elseif x_txt < -0.1, h_align = 'right'; else, h_align = 'center'; end
        text(x_txt, y_txt, labels_spider{i}, 'Color', theme_color, 'FontName', 'Helvetica', 'FontSize', 10, ...
            'FontWeight', 'bold', 'HorizontalAlignment', h_align, 'VerticalAlignment', 'middle', 'Interpreter', 'none');
    end
    hold off;

    disp('Global normalized metrics (All Channels):')
    total_correlations = struct( ...
        'Raw_filtered_data_TC', TC_ch, 'Time_energy_TC', TC_en, 'Time_spectral_entropy_TC', TC_ief, ...
        'Instantaneous_amplitude_TC', TC_eamp, 'Instantaneous_frequency_TC', TC_efreq, ...
        'Instantaneous_phase_TC', TC_ephi, 'TFR_ridge_frequencies_TC', TC_ridge, ...
        'Multivariate_permutation_entropy', MPerm, 'Multivariate_dispersion_entropy', MDisp);
    disp(total_correlations);

    if ~good_quality
        a = annotation('textbox', [0.02, 0.02, 0.96, 0.05], ...
            'String', '*TC on a subset of channels following channel selection are available in the Excel output, provided a sufficient number of channels is retained.', ...
            'EdgeColor', 'none', 'Color', [0.3 0.3 0.3], 'FontSize', 10, 'FontAngle', 'italic');
    end
else
    disp('Total correlation not computed: only one channel in the data');
end
toc
% ********************************************************************
%% 11. EXPORT PARAMETERS TO EXCEL
% ********************************************************************
disp('=== EXPORTING PARAMETERS_RESULTS TO EXCEL ===');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
excel_filename = fullfile(output_folder, 'Single_channel_parameters.xlsx');
first_struct = parameters{1};
if isstruct(first_struct)
    field_names = fieldnames(first_struct);
else
    error('parameters{1} is not a struct. Check your pipeline.');
end
n_fields = numel(field_names);
n_rows   = nCh;
table_cell = cell(n_rows+1, n_fields+5);
table_cell{1,1} = 'Channel';
for f = 1:n_fields
    table_cell{1,1+f} = field_names{f};
end
table_cell{1,n_fields+2} = 'Fraction_of_contaminated_signal'; 
table_cell{1,n_fields+3} = 'is_noise'; 
table_cell{1,n_fields+4} = sprintf('SD_Selection_Strict_%.3fHz', sd_strict); 
table_cell{1,n_fields+5} = sprintf('SD_Selection_Medium_%.3fHz', sd_medium); 
table_cell{1,n_fields+6} = sprintf('SD_Selection_Loose_%.3fHz', sd_loose); 
for ch = 1:nCh
    table_cell{ch+1,1} = channel_names{ch};
    for f = 1:n_fields
        val = parameters{ch}.(field_names{f});
        if isnumeric(val)
            if isscalar(val), table_cell{ch+1,1+f} = val;
            else, table_cell{ch+1,1+f} = mat2str(val); end
        elseif ischar(val) || isstring(val)
            table_cell{ch+1,1+f} = char(val);
        elseif isstruct(val)
            try table_cell{ch+1,1+f} = jsonencode(val); catch, table_cell{ch+1,1+f} = '[struct]'; end
        elseif iscell(val)
            try table_cell{ch+1,1+f} = jsonencode(val); catch, table_cell{ch+1,1+f} = '[cell]'; end
        else
            table_cell{ch+1,1+f} = '[unknown type]';
        end
    end
    
    if ~isempty(r{ch}), frac_contaminate = sum(r{ch}) / length(r{ch});
    else, frac_contaminate = NaN; end
    table_cell{ch+1, n_fields+2} = frac_contaminate;
    
    if all_reliable(ch)
        table_cell{ch+1,n_fields+3} = 'no';
    else
        table_cell{ch+1,n_fields+3} = 'yes';
    end
    
    if keep_strict(ch), table_cell{ch+1,n_fields+4} = 'good'; else, table_cell{ch+1,n_fields+4} = 'bad'; end
    if keep_medium(ch), table_cell{ch+1,n_fields+5} = 'good'; else, table_cell{ch+1,n_fields+5} = 'bad'; end
    if keep_loose(ch), table_cell{ch+1,n_fields+6} = 'good'; else, table_cell{ch+1,n_fields+6} = 'bad'; end
end
try writecell(table_cell, excel_filename); disp(['✔ parameters exported to: ', excel_filename])
catch ME, warning(['❌ Error writing Excel file: ', ME.message]); end
% ********************************************************************
%% 12. EXPORT TOTAL CORRELATIONS AND MULTIVARIATE ENTROPIES 
% ********************************************************************
disp('=== EXPORTING TOTAL CORRELATIONS AND MULTIVARIATE ENTROPIES ===');
excel_filename2 = fullfile(output_folder, 'Total_correlations_and_multivariate_entropies.xlsx');
labels_TC = { ...
    'Raw_filtered_data_total_corr', 'Time_energy_total_corr', 'Time_spectral_entropy_total_corr', ...
    'IA_total_corr', 'IF_total_corr', 'IP_total_corr', 'TFR_ridge_total_corr', ...
    'Multivariate_permutation_entropy', 'Multivariate_dispersion_entropy' ...
};
n_cols = 4 + length(labels_TC);
table_TC = cell(n_scenarios + 1, n_cols);
% Write Headers
table_TC(1, 1:4) = {'Selection_Method', 'N_Good_Channels', 'Selected_Channels', 'Fraction_of_excluded_signal'};
table_TC(1, 5:end) = labels_TC;
% Fill Rows dynamically from tc_results struct generated in Section 9
if nCh > 1
    for sc = 1:n_scenarios
        table_TC{sc+1, 1} = tc_results{sc}.Method;
        table_TC{sc+1, 2} = tc_results{sc}.N_Channels;
        table_TC{sc+1, 3} = tc_results{sc}.Names;
        table_TC{sc+1, 4} = tc_results{sc}.r_tot_frac;
        
        for m = 1:length(labels_TC)
            table_TC{sc+1, 4+m} = tc_results{sc}.Metrics(m);
        end
    end
else
    % Fallback if processing only 1 channel
    table_TC{2,1} = 'All_Channels';
    table_TC(2, 2:end) = {NaN};
end
try writecell(table_TC, excel_filename2); disp(['✔ Total correlations exported to: ', excel_filename2])
catch ME, warning(['❌ Error writing Excel file: ', ME.message]); end
% ********************************************************************
%% 13. EXPORT INSTANTANEOUS AMPLITUDE, PHASE, AND FREQUENCY CURVES 
% ********************************************************************
disp('=== EXPORTING INSTANTANEOUS CURVES ===');
excel_filename3 = fullfile(output_folder, 'Instantaneous_curves.xlsx');
time_vec = t(:);
nT = length(time_vec);
total_columns = 1 + 3*nCh;
table_curves = cell(nT+1, total_columns);
table_curves{1,1} = 'Time';
col_index = 2;
for ch = 1:nCh, table_curves{1, col_index} = sprintf('Amplitude_ch%d', ch); col_index = col_index + 1; end
for ch = 1:nCh, table_curves{1, col_index} = sprintf('Phase_ch%d', ch); col_index = col_index + 1; end
for ch = 1:nCh, table_curves{1, col_index} = sprintf('Frequency_ch%d', ch); col_index = col_index + 1; end
table_curves(2:end,1) = num2cell(time_vec);
col_index = 2;
for ch = 1:nCh
    amp = amp_curve{ch}(:); if length(amp) < nT, amp = [amp; NaN(nT - length(amp),1)]; end
    table_curves(2:end, col_index) = num2cell(amp); col_index = col_index + 1;
end
for ch = 1:nCh
    ph = ph_curve{ch}(:); if length(ph) < nT, ph = [ph; NaN(nT - length(ph),1)]; end
    table_curves(2:end, col_index) = num2cell(ph); col_index = col_index + 1;
end
for ch = 1:nCh
    freq = if_curve{ch}(:); if length(freq) < nT, freq = [freq; NaN(nT - length(freq),1)]; end
    table_curves(2:end, col_index) = num2cell(freq); col_index = col_index + 1;
end
% === MODIFICATO: Summary block espanso a 5 colonne ===
summary_block = cell(nCh+2, 5);
summary_block{1,1} = 'Channel Reliability Summary (Multiple Methods)';
summary_block{2,1} = 'Channel';
summary_block{2,2} = 'MC_Simes_Status';
summary_block{2,3} = sprintf('SD_Selection_Strict_%.3fHz', sd_strict);
summary_block{2,4} = sprintf('SD_Selection_Medium_%.3fHz', sd_medium);
summary_block{2,5} = sprintf('SD_Selection_Loose_%.3fHz', sd_loose);
for ch = 1:nCh
    % Inserimento del nome del canale (corregge il bug della colonna A vuota)
    if exist('channel_names', 'var')
        summary_block{2+ch, 1} = channel_names{ch};
    else
        summary_block{2+ch, 1} = sprintf('Ch %d', ch);
    end
    
    % 1. Monte Carlo Evaluation (Simes)
    if all_reliable(ch) 
        summary_block{2+ch, 2} = 'Reliable';
    else
        summary_block{2+ch, 2} = 'Unreliable: weak normogastric evidence';
    end
    
    % 2.  SD Evaluation with Greedy Dropping
    if keep_strict(ch), summary_block{2+ch, 3} = 'Kept'; else, summary_block{2+ch, 3} = 'Dropped'; end
    if keep_medium(ch), summary_block{2+ch, 4} = 'Kept'; else, summary_block{2+ch, 4} = 'Dropped'; end
    if keep_loose(ch), summary_block{2+ch, 5} = 'Kept'; else, summary_block{2+ch, 5} = 'Dropped'; end
end
[n_summary_rows, n_summary_cols] = size(summary_block);
if n_summary_cols < total_columns
    summary_block(:, n_summary_cols+1 : total_columns) = {''}; 
end
combined_table = [summary_block; cell(1,total_columns); table_curves];
try 
    writecell(combined_table, excel_filename3); 
    disp(['✔ Instantaneous curves + reliability exported to: ', excel_filename3]);
catch ME
    warning(['❌ Error writing Excel file: ', ME.message]); 
end
%% Local helper functions
function y = conditional_decimate(x, mask)
    x = x(:);
    n = min(length(x), length(mask));
    valid = ~mask(1:n);
    x = x(1:n);
    x_valid = x(valid);
    x_valid(~isfinite(x_valid)) = NaN;
    if all(isnan(x_valid)), y = []; return; end
    x_valid = fillmissing(x_valid, 'linear', 'EndValues','nearest');
    if numel(x_valid) > 50
        try y = decimate(x_valid, 2); catch, y = downsample(x_valid, 2); end
    else
        y = downsample(x_valid, 2);
    end
end
function [is_globally_significant, p_simes] = calculate_simes_test(p_values, alpha)
    if nargin < 2 || isempty(alpha), alpha = 0.05; end
    p_values = p_values(~isnan(p_values));
    m = length(p_values);
    if m == 0, is_globally_significant = false; p_simes = NaN; return; end
    p_sorted = sort(p_values(:));
    j_indices = (1:m)';
    thresholds = (j_indices .* alpha) / m;
    rejection_vector = p_sorted <= thresholds;
    is_globally_significant = any(rejection_vector);
    p_simes = min(1, max(0, min((m .* p_sorted) ./ j_indices)));
end
function keep = evaluate_channels_sd(freqs, tol_hz)
    min_fraction = 0.50; 
    n_start = length(freqs);
    keep = true(1, n_start);
    keep(isnan(freqs)) = false; 
    
    c4 = @(N) sqrt(2 / (N - 1)) * exp(gammaln(N / 2) - gammaln((N - 1) / 2));
    Phi = @(x) 0.5 * (1 + erf(x / sqrt(2)));
    d2 = @(N) integral(@(x) 1 - Phi(x).^N - (1 - Phi(x)).^N, -10, 10);
    calc_sd_threshold = @(N) tol_hz * (c4(N) / d2(N));
    
    initial_valid_channels = sum(keep);
    min_required_channels = max(3, ceil(initial_valid_channels * min_fraction));
    
    while sum(keep) > min_required_channels
        current_freqs = freqs(keep);
        n_current = length(current_freqs);
        current_sd = std(current_freqs); 
        current_threshold = calc_sd_threshold(n_current);
        if current_sd <= current_threshold, break; end
        
        active_indices = find(keep);
        best_sd_without_i = Inf;
        worst_idx = NaN;
        
        for i = 1:length(active_indices)
            idx = active_indices(i);
            test_mask = keep;
            test_mask(idx) = false;
            test_freqs = freqs(test_mask);
            if length(test_freqs) < 2, continue; end
            test_sd = std(test_freqs); 
            if test_sd < best_sd_without_i, best_sd_without_i = test_sd; worst_idx = idx; end
        end
        if (sum(keep) - 1) < min_required_channels, break; end
        keep(worst_idx) = false;
    end
    
    final_n = sum(keep);
    if final_n >= min_required_channels
        if std(freqs(keep)) > calc_sd_threshold(final_n), keep(:) = false; end
    else
        keep(:) = false; 
    end
end