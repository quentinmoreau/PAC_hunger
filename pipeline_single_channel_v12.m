function [time_pow, norm_pow, parameters, efreq, eamp, ephi, r1, tp, ief, amp_curve, ph_curve, if_curve, filled_signal, ridge_curve, parameters_results] = pipeline_single_channel_v12(x, t, fs, freq_faslt, channel_id)
% PIPELINE_SINGLE_CHANNEL_V10 - Main processing pipeline for EGG signal analysis
% Performs time-frequency analysis, artifact detection, quality check and parameter extraction
% from a single channel EGG recording
%
% Inputs:
%   x - Input signal
%   t - Time vector
%   fs - Sampling frequency
%   freq_faslt - Frequency range for time-frequency analysis
%   channel_id - Channel identifier for plotting and logging
%
% Outputs:
%   [Multiple outputs related to signal processing results]
    %% Initialization and signal preparation
    tic;
    
    % Perform time-frequency analysis
    [faslt, time_pow, norm_pow] = general_faslt_plot_v3(x, t, fs, freq_faslt, 0.01, 0.2, [], 'on', 'on', 'zero');
    sgtitle(sprintf('FASLT for channel %d', channel_id));
    
    %% Spike-like Artifact detection using surrogate testing
    parameters_results_spikes = MC_test_spikes_RP(x, t, fs, 300, channel_id);
    TsEn_p = parameters_results_spikes.TsEn.p_valueL;
    sig_kurt = kurtosis(x);
    IE_kurt = kurtosis(time_pow);
    disp(['Energy kurtosis: ', num2str(IE_kurt)])
    disp(['Signal kurtosis: ', num2str(sig_kurt)])
    
    %% Main processing branch – spike-like artifact-free vs. contaminated signals
    % After removing spike-like artifacts, the signal may still be affected by other types of noise.
    % These non-impulsive noise sources do not produce sharp, localized spikes,
    % but can still significantly degrade the quality of the signal.
    if TsEn_p >= 0.05 && sig_kurt <= 10 && IE_kurt <= 10
        % Clean channel processing
        [parameters, efreq, eamp, ephi, tp, ief, amp_curve, ph_curve, if_curve, ridge_curve, filled_signal, r1, parameters_results] = ...
            process_clean_channel(x, t, fs, faslt, channel_id);
    else
        % Contaminated channel processing
        [parameters, efreq, eamp, ephi, r1, tp, ief, amp_curve, ph_curve, if_curve, filled_signal, ridge_curve, parameters_results, norm_pow] = ...
            process_contaminated_channel(x, t, fs, freq_faslt, channel_id, TsEn_p, sig_kurt, IE_kurt);
    end
    
    toc;
end

%% Helper functions =======================================================

function [parameters, efreq, eamp, ephi, tp, ief, amp_curve, ph_curve, if_curve, ridge_curve, filled_signal, r1, parameters_results] = ...
    process_clean_channel(x, t, fs, faslt, channel_id)
    %% Initialize all output variables for clean channel
    filled_signal = x;  % Use original signal for clean channel
    r1 = zeros(size(t)); % No anomalies in clean channel
    
    %% Display channel status
    fprintf('---------- Channel %d has no significant artifacts\n', channel_id);
    
    %% Signal quality assessment
    disp('Assessing signal quality against noise:');
    %% Surrogate testing on reconstructed signal
    [parameters_results] = MC_test_IAAWFT_v9(x, t, fs, 300, [], channel_id);
   
    fprintf('p-value (spectral skewness > surrogates): %.4f\n', parameters_results.spec_skewness.p_value_right);
    fprintf('p-value (spectral sparsity > surrogates): %.4f\n', parameters_results.spec_sparsity.p_value_right);
    fprintf('p-value (dominant band power > surrogates): %.4f\n', parameters_results.dom_pow_rel.p_value_right);
    
    [is_globally_significant, p_simes] = calculate_simes_test([parameters_results.spec_skewness.p_value_right, ...
        parameters_results.spec_sparsity.p_value_right], 0.05);
    
    if is_globally_significant
        fprintf('---------- Channel %d has reliable oscillating dynamics\n', channel_id);
    else
        fprintf('---------- Channel %d dynamics likely produced by noise\n', channel_id);
    end
    
    %% Parameter extraction (Call to the new v4!)
    disp('Extracting parameters...');
     [parameters, tp, ief, ridge_curve] = ...
            par_est_TFR_v5(x, t, 'off', [], fs, faslt, x, channel_id);
    
    %% Display extracted parameters
    fprintf('Parameters for channel %d:\n', channel_id);
    disp(parameters);
    
    %% Instanteous normogastric quantities
    f_dom = get_scalar(parameters.Dominant_frequency);   
    opts  = struct('method','adaptive', 'f_dom', f_dom); 
    [efreq, eamp, ephi, ~, amp_curve, ph_curve, if_curve] = ...
    wavelet_inst_extr_v5(x, t, [], 'on', fs, opts);
    if is_globally_significant
          sgtitle(sprintf('Fundamental gastric component - Channel %d', channel_id));
    else
        % Create the sgtitle and store the handle
        sgtitle(sprintf('Fundamental gastric component - Channel %d (Low reliability)', channel_id));
        % Add warning text in the figure
        annotation('textbox', [0.15 0.89 0.7 0.05], ...
        'String', '⚠ Estimated curves may be unreliable due to high noise or weak normogastric spectral evidence', ...
        'EdgeColor', 'none', ...
        'Color', [0.8 0.1 0.1], ...
        'FontSize', 11, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');
    end
    
    %% Create spider plot visualization
    create_spider_plot(parameters, channel_id, is_globally_significant);
end

function [parameters, efreq, eamp, ephi, r1, tp, ief, amp_curve, ph_curve, if_curve, filled_signal, ridge_curve, parameters_results, norm_pow] = ...
    process_contaminated_channel(x, t, fs, freq_faslt, channel_id, TsEn_p, sig_kurt, IE_kurt)
    
    %% Display artifact information
    fprintf('-------------- Artifacts detected in channel %d (TsEn_p=%.3f, sig_kurt=%.1f, IE_kurt=%.1f)\n', ...
            channel_id, TsEn_p, sig_kurt, IE_kurt);
            
    %% Artifact removal and reconstruction
    [filled_signal, anomaly_intervals, faslt, time_pow, norm_pow] = ...
        find_anomalies_and_plot_TFR_v5(x, t, 1000, fs, 100, channel_id);
    sgtitle(sprintf('Cleaned FASLT for channel %d', channel_id));
    
    %% Create anomaly indicator vector
    r1 = interval_to_vector(anomaly_intervals, length(t));
    
    %% Check if reconstruction is possible
    if any(r1 == 0) 
        %% Parameter extraction from reconstructed signal (Call to the new v4!)
        disp('Extracting parameters from reconstructed signal...');
        [parameters, tp, ief, ridge_curve] = ...
            par_est_TFR_v5(x, t, 'off', anomaly_intervals, fs, faslt, filled_signal, channel_id);
        
        %% Surrogate testing on reconstructed signal
        [parameters_results] = MC_test_IAAWFT_v9(x, t, fs, 300, anomaly_intervals, channel_id);
        
        %% Display quality assessment
        plotTitle_s = sprintf('Artifact detection via Surrogate test for channel %d', channel_id);
        sgtitle(plotTitle_s, 'Interpreter', 'none');
        set(findall(gcf, 'Type', 'axes'), 'TitleFontSizeMultiplier', 0.7);
        
        disp('Assessing signal quality against noise:');
        fprintf('p-value (spectral skewness > surrogates): %.4f\n', parameters_results.spec_skewness.p_value_right);
        fprintf('p-value (spectral sparsity > surrogates): %.4f\n', parameters_results.spec_sparsity.p_value_right);
        fprintf('p-value (dominant band power > surrogates): %.4f\n', parameters_results.dom_pow_rel.p_value_right);
    
        [is_globally_significant, p_simes] = calculate_simes_test([parameters_results.spec_skewness.p_value_right, ...
        parameters_results.spec_sparsity.p_value_right], 0.05);
    
        if is_globally_significant
            fprintf('---------- Channel %d has reliable oscillating dynamics after reconstruction\n', channel_id);
        else
            fprintf('---------- Channel %d dynamics likely noise after reconstruction\n', channel_id);
        end
        
        %% Instanteous normogastric quantities
        f_dom = get_scalar(parameters.Dominant_frequency);
        opts  = struct('method','adaptive', 'f_dom', f_dom); 
        [efreq, eamp, ephi, ~, amp_curve, ph_curve, if_curve] = ...
        wavelet_inst_extr_v5(filled_signal, t, anomaly_intervals, 'on', fs, opts);
        if is_globally_significant
            sgtitle(sprintf('Fundamental gastric component - Channel %d', channel_id));
        else
            % Create the sgtitle and store the handle
            sgtitle(sprintf('Fundamental gastric component - Channel %d (Low reliability)', channel_id));
            % Add warning text in the figure
            annotation('textbox', [0.15 0.89 0.7 0.05], ...
            'String', '⚠ Estimated curves may be unreliable due to high noise or weak normogastric spectral evidence', ...
            'EdgeColor', 'none', ...
            'Color', [0.8 0.1 0.1], ...
            'FontSize', 11, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center');
        end
        
        %% Create spider plot visualization
        create_spider_plot(parameters, channel_id, is_globally_significant);
    else
        %% Handle unrecoverable signal
        disp('Channel too compromised by artifacts - no reliable TFR possible');
    end
end

function create_spider_plot(parameters, channel_id, is_reliable)
% CREATE_SPIDER_PLOT - Versione Aggiornata per par_est_TFR_v5 (32 parametri)
% Disegna lo spider plot allineato alla perfezione con la nuova struct in uscita.
    %% 1. Check preliminary
    if isempty(fieldnames(parameters))
        disp('No parameters available for spider plot');
        return;
    end
    
    par_names = fieldnames(parameters);
    num_params = numel(par_names);
    
    %% 2. Graphic Configuration (Grouped by Logical Domain)
    % Order: {Display Label, [Min, Max], Scaling, Precision, TickFormat}
    config_db = {
        % --- FREQUENCY & GLOBAL POWER ---
        'Dom Freq',        [0.016, 0.16],  'linear', 3, '%.3f';  
        'Dom Power',       [100, 1e10],    'log',    0, '%.0e';  
        'Centroid',        [0, 0.16],      'linear', 3, '%.3f';  
        'Spread',          [0, 0.072],     'linear', 3, '%.3f';  
        
        % --- ABSOLUTE POWER PER BAND ---
        'Tot Normo',       [1e-6, 1e5],     'log',    0, '%.0f';  
        'Tot Bradi',       [1e-6, 1e5],     'log',    0, '%.0f';  
        'Tot Tachi',       [1e-6, 1e5],     'log',    0, '%.0f';  
        
        % --- RELATIVE POWER PER BAND ---
        'Rel Normo',       [0, 1],         'linear', 2, '%.2f';  
        'Rel Bradi',       [0, 1],         'linear', 2, '%.2f';  
        'Rel Tachi',       [0, 1],         'linear', 2, '%.2f';  
        'Rel Dominant',    [1e-4, 100],    'log',    1, '%.1f';  
        'Rel Sup Dom',     [1e-4, 10],     'log',    1, '%.1f';  
        'Rel Sub Dom',     [1e-4, 10],     'log',    1, '%.1f';  
        
        % --- POWER SPECTRUM ENTROPIES ---
        'Spectral En',     [0, 1],         'linear', 2, '%.2f';  
        'PS Tsallis',      [0, 1],         'linear', 2, '%.2f';  
        
        % --- TIME-FREQUENCY (TFR) STATISTICS ---
        'Normo %',         [0, 1],         'linear', 2, '%.2f';
        'Mn Ridge',        [0, 0.16],      'linear', 3, '%.3f';  
        'Ridge SD',        [0, 0.16],      'linear', 3, '%.3f';  
        'Mn Energy',       [100, 1e10],    'log',    0, '%.0e';  
        'En SD',           [100, 1e10],    'log',    0, '%.0e';  
        
        % --- TIME-FREQUENCY ENTROPIES & SPARSITY ---
        'TFR Hoyer',           [0, 1],         'linear', 2, '%.2f';  
        'TFR Renyi',           [0, 1],         'linear', 2, '%.2f';  
        'TFR GradEn',          [0, 1],         'linear', 2, '%.2f';  
        'Rel Dom GradEn',      [0, 5],         'linear', 2, '%.2f';  
        
        % --- TIME DOMAIN ENTROPIES ---
        'Sig DispEn',      [0, 1],         'linear', 1, '%.1f';  
        'Sig SampEn',      [0, 5],         'linear', 1, '%.1f';  
        'IE DispEn',       [0, 1],         'linear', 1, '%.1f';  
        'IE SampEn',       [0, 5],         'linear', 1, '%.1f';  
        
        % --- NONLINEAR METRICS & COMPLEXITY ---
        'Corr Dim',        [-5, 5],        'linear', 1, '%.1f';  
        'Lyap Exp',        [0, 5],         'linear', 2, '%.2f';  
        'Act Inf',         [0, 5],         'linear', 2, '%.2f';  
        'DFA Exp',         [0, 2],         'linear', 2, '%.2f'   
    };
    %% 3. Preparazione Dati
    p_values = zeros(1, num_params);
    final_labels = {};      
    dummy_labels = {};      
    final_limits = [];
    final_scaling = {};
    final_precision = [];
    final_tick_format = {};
    
    for i = 1:min(num_params, size(config_db, 1))
        val = parameters.(par_names{i});
        if isnan(val) || isinf(val), val = 0; end
        p_values(i) = val;
        
        final_labels{i} = config_db{i, 1}; %#ok<AGROW>
        dummy_labels{i} = '';              %#ok<AGROW>
        final_limits = [final_limits, config_db{i, 2}']; %#ok<AGROW>
        final_scaling{i} = config_db{i, 3}; %#ok<AGROW>
        final_precision(i) = config_db{i, 4}; %#ok<AGROW>
        final_tick_format{i} = config_db{i, 5}; %#ok<AGROW>
    end
    %% 4. Impostazioni Colori
    if is_reliable
        theme_color = [0 0.5 0];      % Verde Scuro
        data_color  = [0 0.6 0];      % Verde Dati
        subtitle_str = 'Signal Quality: Reliable (Green)';
    else
        theme_color = [0.7 0 0];      % Rosso Scuro
        data_color  = [0.8 0 0];      % Rosso Dati
        subtitle_str = 'Signal Quality: Unreliable/Noisy (Red)';
    end
    
    edge_col = repmat('w', 1, numel(p_values)); 
    lbl_offset = 0.25; 
    f = figure('Color', 'w');
    f.Position(3:4) = [950, 750]; 
    
    %% 5. Plot Principale
    spider_plot(p_values, ...
        'AxesLabels', dummy_labels, ... 
        'AxesLimits', final_limits, ...
        'AxesScaling', final_scaling, ...
        'AxesPrecision', final_precision, ...
        'AxesTickFormat', final_tick_format, ...
        'AxesLabelsEdge', edge_col, ...
        'AxesLabelsOffset', lbl_offset, ... 
        'AxesFontSize', 8, ...
        'AxesFontColor', [0.2 0.2 0.2], ... 
        'AxesDisplay', 'none', ...          
        'AxesInterval', 2, ...      
        'FillOption', {'on'}, ...
        'FillTransparency', 0.15, ...       
        'Color', data_color, ...            
        'AxesColor', [0.7 0.7 0.7], ...
        'LineStyle', {'-'}, ...
        'LineWidth', 2, ...                 
        'Marker', {'o'}, ...
        'MarkerSize', 5);
    
    title(sprintf('Parameters for channel %d', channel_id), 'FontSize', 14);
    subtitle(subtitle_str, 'Color', theme_color, 'FontWeight', 'bold');
    %% 6. Disegno Manuale Etichette e Guide
    hold on;
    theta = linspace(pi/2, pi/2 - 2*pi, num_params + 1);
    theta(end) = []; 
    
    r_start = 1.02; 
    r_end = 1 + lbl_offset - 0.02; 
    
    for i = 1:num_params
        % 1. Disegna Linea Guida
        [x_guide, y_guide] = pol2cart(theta(i), [r_start, r_end]);
        plot(x_guide, y_guide, ':', 'Color', theme_color, 'LineWidth', 0.8);
        
        % 2. Disegna Testo (Manualmente)
        [x_txt, y_txt] = pol2cart(theta(i), r_end + 0.05);
        
        if x_txt > 0.1
            h_align = 'left';  
        elseif x_txt < -0.1
            h_align = 'right'; 
        else
            h_align = 'center'; 
        end
        
        text(x_txt, y_txt, final_labels{i}, ...
            'Color', theme_color, ...       
            'FontName', 'Helvetica', ...
            'FontSize', 9, ...              
            'FontWeight', 'bold', ...
            'HorizontalAlignment', h_align, ...
            'VerticalAlignment', 'middle', ...
            'Interpreter', 'none');
    end
    hold off;
end

function [is_globally_significant, p_simes] = calculate_simes_test(p_values, alpha)
% CALCULATE_SIMES_TEST Global hypothesis testing using Simes' procedure.
% (Invariato rispetto all'originale)
    if nargin < 2 || isempty(alpha)
        alpha = 0.05;
    end
    % Remove any potential NaN values
    p_values = p_values(~isnan(p_values));
    m = length(p_values);
    
    if m == 0
        is_globally_significant = false;
        p_simes = NaN;
        return;
    end
    % Sort p-values in ascending order
    p_sorted = sort(p_values(:));
    
    % Compute the Simes threshold for each ordered p-value: (j * alpha) / m
    j_indices = (1:m)';
    thresholds = (j_indices .* alpha) / m;
    
    % Check if any p-value meets the Simes criterion
    rejection_vector = p_sorted <= thresholds;
    is_globally_significant = any(rejection_vector);
    
    % Optional: Calculate the adjusted Simes p-value for reporting
    p_simes = min((m .* p_sorted) ./ j_indices);
    
    % Ensure the p-value is bounded between 0 and 1
    p_simes = min(1, max(0, p_simes));
end


function v = get_scalar(S) %#ok<DEFNU>
% Convenience: pull a scalar dominant frequency out of a struct field such
% as parameters.Dominant_frequency (which may itself be a struct). Tries a
% few common subfield names, else the first numeric scalar found.
if isnumeric(S) && isscalar(S), v = S; return; end
if isstruct(S)
    for name = {'value','val','freq','f','dominant','Hz'}
        if isfield(S, name{1}) && isscalar(S.(name{1})) && isnumeric(S.(name{1}))
            v = S.(name{1}); return;
        end
    end
    fn = fieldnames(S);
    for k = 1:numel(fn)
        if isnumeric(S.(fn{k})) && isscalar(S.(fn{k})), v = S.(fn{k}); return; end
    end
end
error('get_scalar:type', 'Could not extract a scalar dominant frequency; pass opts.f_dom explicitly.');
end