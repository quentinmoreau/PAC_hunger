function [faslt, time_pow, norm_pow] = clean_faslt_plot_no_filt_v3(x, t, fs, freq_faslt, anomalies, diagno, enlarge, last_plot)

% Compute raw FASLT without filtering
[faslt, time_pow, norm_pow, ~] = general_faslt_plot_no_filt(x, t, fs, freq_faslt, 0.01, 0.2, [], 'off', diagno, 'zero');

%if ~isempty(anomalies)
    % Zero out anomalies in FASLT
    for i = 1:2:length(anomalies)-1
        idx_start = max(anomalies(i) - enlarge, 1);
        idx_end   = min(anomalies(i+1) + enlarge, length(x));
        faslt(:, idx_start:idx_end) = 0;
    end

    % Recompute energy metrics
    faslt_abs = abs(faslt)/2;
    norm_pow  = sum(faslt_abs.^2, 2);
    freq_faslt = freq_faslt(:); % ensure column vector
    time_pow  = sum(faslt_abs.^2 ./ freq_faslt, 1);

    for i = 1:2:length(anomalies)-1
        idx_start = max(anomalies(i) - enlarge, 1);
        idx_end   = min(anomalies(i+1) + enlarge, length(x));
        time_pow(idx_start:idx_end) = NaN;
    end

    % Compute summary statistics

    peak_pow = max(norm_pow);
peak_freq = freq_faslt(norm_pow == peak_pow);
mean_time = mean(time_pow, 'omitnan');

% Entropy measures
norm_pow_prob = norm_pow / sum(norm_pow); 
[~, PS_tsallis_entropy] = minim_Tallis(norm_pow_prob, 500, 'off');

%[~, eLag, eDIm] = phaseSpaceReconstruction(norm_pow, 'MaxDim', 9, 'MaxLag', 100);


%De    = DispEn(norm_pow, 'm', eDIm, 'tau', eLag);



column_mask = any(faslt ~= 0, 1);
valid_TFR = faslt(:, column_mask);
t_valid = t(column_mask);
R = renyi(abs(valid_TFR).^2, t_valid, freq_faslt(:), 3) / log2(numel(valid_TFR));


% Format values with controlled precision


text_str = sprintf(['Dominant freq: %.4f Hz\n' ...
                    'Dominant pow: %.2e\n' ...
                    'PSD Tsallis: %.4f\n' ...
                    'Mean Energy: %.2e\n' ...
                    'Rényi entropy: %.4f'], ...
                    peak_freq, peak_pow, PS_tsallis_entropy, mean_time, R);


%end

%% DIAGNOSTIC PLOTS
if strcmp(diagno, 'on')
    figure('Name','FASLT Diagnostic','NumberTitle','off');
    subplot(2, 1, 1);
    plot(t, x);
    hold on;
    if ~isempty(anomalies)
        xline(t(anomalies), '--k');
    end
    ylabel('\muV');
    title('Raw signal');
    grid on;
    set(gca, 'FontSize', 16);

    subplot(2, 1, 2);
    tplot(faslt, t, freq_faslt);
    hold on;
    if ~isempty(anomalies)
        xline(t(anomalies), '--k');
    end
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('Cleaned FASLT');
    grid on;
    set(gca, 'FontSize', 16);
end

%% FINAL COMPACT PLOT
if strcmp(last_plot, 'on')
    figure('Name','FASLT Analysis','NumberTitle','off','Position',[100 100 1300 700]);
    tlayout = tiledlayout(5,6,'TileSpacing','compact','Padding','compact');

    %% TFR
    ax_tfr = nexttile([4 5]);
    faslt_plot = abs(faslt)/2;

    % Highlight anomalies with a uniform color
    artifact_mask = false(1, length(t));
    for k = 1:2:length(anomalies)
        artifact_mask(anomalies(k):anomalies(k+1)) = true;
    end
    max_val = max(faslt_plot(:));
    faslt_plot(:, artifact_mask) = max_val * 1.2;

    imagesc(t, freq_faslt, faslt_plot);
    set(gca, 'YDir','normal');
    cmap = 1 - slanCM('gist_rainbow'); 
    cmap_smooth = 0.5 * cmap + 0.5 * ones(size(cmap)); 
    colormap(ax_tfr, cmap_smooth);
    colorbar('eastoutside');
    ylabel('Frequency (Hz)');
    title('FASLT Time-Frequency Representation');
    set(gca, 'FontSize', 13);
    grid on;
    hold on;
    if ~isempty(anomalies)
        xline(t(anomalies), '--k', 'LineWidth', 1.5);
    end
    hold off;

    %% Power Spectrum
    ax_ps = nexttile([4 1]);
    plot(norm_pow, freq_faslt, 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
    grid on;
    ylim([min(freq_faslt) max(freq_faslt)]);
    %xlabel('Normalized Power');
    set(gca, 'YDir','normal');
    title('Power spectrum (PS)');
    set(gca, 'FontSize', 13, 'LineWidth', 1.2);

    %% Instantaneous Energy
    ax_ie = nexttile([1 5]);
    plot(t, time_pow, 'LineWidth', 2, 'Color', [0 0 0]);
    hold on;
    
    yl = ylim;
    
    for k = 1:2:length(anomalies)
        x_patch = [t(anomalies(k)) t(anomalies(k+1)) t(anomalies(k+1)) t(anomalies(k))];
        y_patch = [yl(1)-1000 yl(1)-1000 yl(2)+1000 yl(2)+1000];
        patch(x_patch, y_patch, [0.85 0.33 0.10], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    end
    
    % Riplotto la linea per sicurezza (puoi anche rimuovere questa seconda chiamata se vuoi)
    plot(t, time_pow, 'LineWidth', 2, 'Color', [0 0 0]);
    hold off;
    grid on;
    xlabel('Time (s)');
    %ylabel('Instantaneous Energy (IE)');
    title('Instantaneous Energy (IE)');
    set(gca, 'FontSize', 13, 'LineWidth', 1.2);
    xlim([t(1) t(end)]);
    hold on;
    if ~isempty(anomalies)
        xline(t(anomalies), '--r', 'LineWidth', 1.5);
    end
    hold off;

    %% Statistics box 
    ax_stats = nexttile(30); %lower right corner (tile 30 = row 5, col 6)
    axis off;
    text(-0.35, 1.2, text_str, ...
    'FontSize', 11, 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'Interpreter', 'none', 'BackgroundColor', 'w', 'EdgeColor', 'k', ...
    'Margin', 6);
end


end
