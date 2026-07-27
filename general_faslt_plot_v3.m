function [faslt, time_pow, norm_pow, buffer] = general_faslt_plot_v3(x, t, fs, freq_faslt, fmin, fmax, anomalies, ssa, plots, padding)
%GENERAL_FASLT_PLOT Calculates FASLT time-frequency representation and related powers.
%
% INPUTS:
%   x          - input signal vector
%   t          - time vector corresponding to x
%   fs         - sampling frequency
%   freq_faslt - frequencies used in FASLT calculation
%   fmin       - minimum frequency of interest (not used here, placeholder)
%   fmax       - maximum frequency of interest (not used here, placeholder)
%   anomalies  - vector of anomaly indices (start and end pairs)
%   ssa        - string 'on' or 'off' to apply SSA detrending
%   plots      - string 'on' or 'off' to enable plotting
%   padding    - padding parameter for FASLT calculation
%
% OUTPUTS:
%   faslt      - FASLT complex time-frequency representation
%   time_pow   - energy per unit time (time marginal of power)
%   norm_pow   - power spectrum (frequency marginal of power)
%   buffer     - auxiliary output from FASLT calculation (not detailed here)

%% Step 1: Optional SSA detrending
if strcmp(ssa, 'on')
    
    % In extremely clean signals SSA may misidentify the fundamental oscillation as a
    % low-frequency trend. To avoid this, we superimpose a small linear trend,
    % ensuring that SSA isolates this artificial trend instead of the true
    % oscillatory component.
    signal_amp = max(x) - min(x);
    trend_linear = (signal_amp / (t(end) - t(1))) * (t - t(1));
    x_trended = x + trend_linear;

    [LT, ~, ~] = trenddecomp(x_trended, 'ssa', 60 * fs); % Decompose signal trend with SSA
    x = x_trended - LT'; % Remove long-term trend
end

%% Step 2: Filtering the signal with Chebyshev filter tailored for EGG signals
x = cheb_EGG_filt(x, fs);

%% Step 3: FASLT calculation depending on presence of anomalies
if isempty(anomalies)
    % No anomalies - process the full clean signal
    [faslt, buffer] = faslt_complex_calc_no_filt(x, fs, freq_faslt, padding);

else
    % Anomalies present - process clean segments between artifacts
    
    faslt_all = [];
    
    % Process clean segment before first anomaly if any
    if anomalies(1) > 1
        x_segment = x(1 : anomalies(1)-1);
        [faslt_segment, buffer] = faslt_complex_calc_no_filt(x_segment, fs, freq_faslt, padding);
        faslt_all = [faslt_all, faslt_segment];
        disp(['The number of artifact processed is 1 over ', num2str(length(anomalies)/2)]);
    else
        disp(['The number of artifact processed is 1 over ', num2str(length(anomalies)/2)]);
    end
    
    % Process intermediate clean segments between anomalies
    if length(anomalies) > 2
        for i = 2 : 2 : length(anomalies) - 1
            % Extract segment between anomaly intervals
            x_segment = x(anomalies(i) : anomalies(i+1));
            [faslt_segment, buffer] = faslt_complex_calc_no_filt(x_segment, fs, freq_faslt, padding);
            
            % Fill gap between current and previous anomaly with zeros
            gap_length = length(anomalies(i-1) : anomalies(i)) - 2;
            zero_fill = zeros(length(freq_faslt), gap_length);
            
            faslt_all = [faslt_all, zero_fill, faslt_segment];
            disp(['The number of artifact processed is ', num2str(i/2 + 1), ' over ', num2str(length(anomalies)/2)]);
        end
        
        % Fill last anomaly interval with zeros
        zero_fill_last = zeros(length(freq_faslt), length(anomalies(end-1) : anomalies(end)) - 1);
        faslt_all = [faslt_all, zero_fill_last];
        
        % Process segment after last anomaly if exists
        if anomalies(end) ~= length(x)
            x_segment = x(anomalies(end) : end);
            [faslt_segment, buffer] = faslt_complex_calc_no_filt(x_segment, fs, freq_faslt, padding);
            faslt_all = [faslt_all, faslt_segment];
        else
            % Append a column of zeros to maintain size consistency
            faslt_all = [faslt_all, zeros(size(faslt_all, 1), 1)];
        end
        
    else
        % Only one anomaly interval present
        if anomalies(2) ~= length(x)
            x_segment = x(anomalies(2) : end);
            [faslt_segment, buffer] = faslt_complex_calc_no_filt(x_segment, fs, freq_faslt, padding);
            
            zero_fill = zeros(length(freq_faslt), length(anomalies(1) : anomalies(2)) - 1);
            faslt_all = [faslt_all, zero_fill, faslt_segment];
        else
            zero_fill = zeros(length(freq_faslt), length(anomalies(1) : anomalies(2)));
            faslt_all = [faslt_all, zero_fill];
        end
    end
    
    faslt = faslt_all;
end

%% Step 4: Ensure freq_faslt is a column vector for correct math
if size(freq_faslt, 1) == 1
    freq_faslt = freq_faslt';
end

%% Step 5: Calculate power spectra and energy over time

if isempty(anomalies)
    % No anomalies - use full faslt
    norm_pow = sum((abs(faslt)/2).^2, 2);              % Power spectrum (frequency marginal)
    time_pow = sum(((abs(faslt)/2).^2) ./ freq_faslt, 1); % Energy per unit time (time marginal)
else
    % Zero out faslt values within anomalous intervals
    for i = 1 : 2 : length(anomalies) - 1
        faslt(:, anomalies(i) : anomalies(i+1)) = 0;
    end
    
    norm_pow = sum((abs(faslt)/2).^2, 2);
    time_pow = sum(((abs(faslt)/2).^2) ./ freq_faslt, 1);
    
    % Mark anomalous intervals in time_pow as NaN for plotting clarity
    for i = 1 : 2 : length(anomalies) - 1
        time_pow(anomalies(i) : anomalies(i+1)) = NaN;
    end
end

%% Step 6: Calculate summary statistics for display
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


%% Step 7: Optional plotting of results

if strcmp(plots, 'on')
    figure('Name','FASLT Analysis','NumberTitle','off','Position',[100 100 1300 700]);
    %figure('Name','FASLT Analysis','NumberTitle','off','Position',[200 200 800 500]);


    tlayout = tiledlayout(5,6,'TileSpacing','compact','Padding','compact');

    % FASLT TFR - top left large area (4 rows x 5 columns)
    ax_tfr = nexttile([4 5]);

    % Prepare FASLT matrix with uniform artifact values and no discontinuities
    faslt_plot = abs(faslt)/2; % normalize

    % Create artifact mask to highlight with uniform color
    artifact_mask = false(1, length(t));
    for k = 1:2:length(anomalies)
        artifact_mask(anomalies(k):anomalies(k+1)) = true;
    end

    % Replace artifact regions with a value above the max to ensure visual contrast
    max_val = max(faslt_plot(:));
    faslt_plot(:, artifact_mask) = max_val * 1.2;

    imagesc(t, freq_faslt, faslt_plot);
    set(gca,'YDir','normal');
    
    % Use a softened reversed colormap
    cmap = 1 - slanCM('gist_rainbow'); 
    alpha = 0.5; 
    cmap_smooth = alpha * cmap + (1 - alpha) * ones(size(cmap)); 
    colormap(ax_tfr, cmap_smooth);
    
    colorbar('eastoutside');
    ylabel('Frequency (Hz)');
    title('FASLT Time-Frequency Representation');
    set(gca, 'FontSize', 13);
    grid on;

    % Add vertical dashed lines at anomaly borders
    hold on;
    for k = 1:length(anomalies)
        xline(t(anomalies(k)), '--k', 'LineWidth', 1.5);
    end
    hold off;

    % Power Spectrum - right column, aligned vertically with TFR (col 6, rows 1–4)
    ax_ps = nexttile([4 1]);
    plot(norm_pow, freq_faslt, 'LineWidth', 2, 'Color', [0 0.4470 0.7410]); % blue
    grid on;
    ax_ps.GridLineStyle = ':';
    ax_ps.GridAlpha = 0.3;
    ylim([min(freq_faslt) max(freq_faslt)]);
    %xlabel('{\mu V}^2/Hz');
    %hXLabel = get(ax_ps, 'XLabel'); 
    %set(hXLabel, 'Units', 'normalized');
    %pos = get(hXLabel, 'Position');
    %disp(pos)  
    %set(hXLabel, 'Position', [pos(1)-0.5, pos(2), 0]);
    set(gca,'YDir','normal');
    title('Power spectrum (PS)');
    set(gca, 'FontSize', 13, 'LineWidth', 1.2);

    % No vertical lines here; the y-axis is frequency
    hold(ax_ps, 'on');
    for k = 1:length(anomalies)
        yline(freq_faslt(1), 'LineStyle', 'none'); 
    end
    hold(ax_ps, 'off');

    % Instantaneous Energy (IE) - bottom plot (last row, 5 columns)
    ax_ie = nexttile([1 5]);
    plot(t, time_pow, 'LineWidth', 2, 'Color', [0 0 0]); % black
    hold on;

    % Add shaded artifact bands behind the curve
    yl = ylim;
    for k = 1:2:length(anomalies)
        x_patch = [t(anomalies(k)) t(anomalies(k+1)) t(anomalies(k+1)) t(anomalies(k))];
        y_patch = [yl(1)-1000 yl(1)-1000 yl(2)+1000 yl(2)+1000]; % wider than actual limits
        patch(x_patch, y_patch, [0.85 0.33 0.10], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    end

    % Replot curve on top to ensure visibility
    plot(t, time_pow, 'LineWidth', 2, 'Color', [0 0 0]);
    hold off;

    grid on;
    xlabel('Time (s)');
    %ylabel('{\mu V}^2');
    title('Instantaneous Energy (IE)');
    set(gca, 'FontSize', 13, 'LineWidth', 1.2);
    xlim([t(1) t(end)]);

    % Add vertical red dashed lines at anomaly borders
    hold on;
    for k = 1:length(anomalies)
        xline(t(anomalies(k)), '--r', 'LineWidth', 1.5);
    end
    hold off;


    % Stats box 
    ax_stats = nexttile(30); %lower right corner (tile 30 = row 5, col 6)
    axis off;
    text(-0.35, 1.2, text_str, ...
    'FontSize', 11, 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'Interpreter', 'none', 'BackgroundColor', 'w', 'EdgeColor', 'k', ...
    'Margin', 6);


end
