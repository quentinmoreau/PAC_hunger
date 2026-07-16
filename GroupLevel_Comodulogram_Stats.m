%% =========================================================================
%% GROUP-LEVEL COMODULOGRAM ANALYSIS UGP: GRAND AVERAGE, PLOTTING, STATS (PRE vs POST)
%% =========================================================================
% Loads the per-subject/per-condition comodulograms (MI_comod: chan x amp x
% phase) produced by the comodulogram pipeline, packages them as FieldTrip
% pseudo-TFR structures, computes grand averages, plots the comodulogram
% over chosen electrodes and as a scalp topography, and runs a paired
% cluster-based permutation test (Post vs Pre) across channels x amplitude
% frequency x phase frequency.
%
% KEY TRICK: FieldTrip has no native "freq x freq" data type. We treat the
% comodulogram as a pseudo time-frequency representation (TFR)

clear; 
% eeglab nogui;
% ft_defaults;
% 
%% =========================================================================
%% 1. CONFIGURATION
%% =========================================================================
base_folder = '/Users/qmoreau/Documents/Work/UGO_PAC/';
output_root = fullfile(base_folder, 'PAC_PreVsPost_Comparison');

dirs = [];
dirs.stats = fullfile(output_root, 'Group_Stats');   % where per-subject Comodulogram.mat files live
dirs.group = fullfile(output_root, 'Group_Level');   % new outputs go here
if ~exist(dirs.group, 'dir'), mkdir(dirs.group); end

subjects = {'P001','P002', 'P003', 'P004', 'P005', 'P006', 'P007', 'P008', 'P009', 'P010', ...
            'P011', 'P012', 'P013', 'P014', 'P015', 'P017', 'P018', 'P019', 'P021', ...
            'P022', 'P023', 'P024', 'P025', 'P026', 'P027', 'P028', 'P029', 'P030', ...
            'P031', 'P032','P033', 'P034', 'P035', 'P036', 'P037', 'P038', 'P039','P040', 'P0166', 'P0200'};

conditions = {'pre', 'post'};

% Electrodes to average over for the line-style comodulogram plot
electrodes_of_interest = {'Pz', 'Cz', 'CPz'};

% Amplitude/phase band to inspect in the topography plot (edit to match
% your effect of interest, e.g. alpha amplitude / gastric peak phase)
band_of_interest_amp   = [8 13];       % Hz (EEG amplitude axis)
band_of_interest_phase = [0.033 0.066];% Hz (EGG phase-driver axis)

% Cluster-based permutation settings
n_permutations = 1000;   % use 10000 for the final/reported run
clusteralpha   = 0.05;
alpha_level    = 0.025;  % two-tailed cluster-level alpha

%% =========================================================================
%% 2. LOAD PER-SUBJECT COMODULOGRAMS -> FIELDTRIP PSEUDO-TFR STRUCTS
%% =========================================================================
data_pre  = {};
data_post = {};
used_subjects = {};

for i = 1:length(subjects)
    subj = subjects{i};
    subj_ok = true;
    tmp_data = struct('pre', [], 'post', []);

    for c = 1:2
        cond_name = conditions{c};
        fname = fullfile(dirs.stats, sprintf('%s_%s_Comodulogram.mat', subj, cond_name));

        if ~exist(fname, 'file')
            fprintf('[!] Comodulogram missing for %s %s. Skipping subject.\n', subj, cond_name);
            subj_ok = false;
            continue;
        end

        S = load(fname, 'MI_comod', 'chan_labels', 'EEG_chanlocs', ...
                  'phase_centers', 'amp_centers');

        % --- Package as a FieldTrip "pseudo-TFR" structure ---
        % dimord 'chan_freq_time':
        %   freq  <- EEG amplitude-frequency axis
        %   time  <- EGG phase-driver-frequency axis (reused field, NOT real time!)
        ft_pac = [];
        ft_pac.label     = S.chan_labels(:);
        ft_pac.dimord    = 'chan_freq_time';
        ft_pac.freq      = S.amp_centers(:)';     % Hz, EEG amplitude
        ft_pac.time      = S.phase_centers(:)';   % Hz, EGG phase (stored as "time")
        ft_pac.powspctrm = S.MI_comod;            % chan x amp x phase -> matches dimord

        tmp_data.(cond_name) = ft_pac;
    end

    if subj_ok
        data_pre{end+1}      = tmp_data.pre;
        data_post{end+1}     = tmp_data.post;
        used_subjects{end+1} = subj;
        fprintf(' -> %s: comodulogram pair loaded.\n', subj);
    end
end

nSubj = length(used_subjects);
fprintf('\nLoaded %d complete Pre/Post comodulogram pairs.\n', nSubj);
if nSubj == 0
    error('No complete subject pairs found - check file paths/names.');
end

% save(fullfile(dirs.group, 'Comodulogram_AllSubjects.mat'), ...
%      'data_pre', 'data_post', 'used_subjects', 'layout', '-v7.3');
% 
%% =========================================================================
%% 3. GRAND AVERAGE ACROSS SUBJECTS (FOR PLOTTING ONLY)
%% =========================================================================
% keepindividual='no' collapses across subjects -> used purely for
% visualisation. The statistics in Section 5 use the per-subject cell
% arrays (data_pre / data_post) directly, as ft_freqstatistics expects.
cfg = [];
cfg.parameter = 'powspctrm';
cfg.keepindividual = 'no';

GA_pre  = ft_freqgrandaverage(cfg, data_pre{:});
GA_post = ft_freqgrandaverage(cfg, data_post{:});

% save(fullfile(dirs.group, 'GrandAverage_PrePost.mat'), 'GA_pre', 'GA_post');
% fprintf('Grand averages computed (N=%d).\n', nSubj);

%% =========================================================================
%% 4. PLOT GRAND-AVERAGE COMODULOGRAM OVER SPECIFIC ELECTRODES
%% =========================================================================
% ft_singleplotTFR averages across the selected channels and plots the 2D
% map. Axes are relabelled manually since "time" here is really phase
% frequency, not time.
cfg = [];
cfg.channel   = electrodes_of_interest;
cfg.parameter = 'powspctrm';
cfg.colorbar  = 'yes';
cfg.zlim      = 'maxabs';
cfg.layout = 'EEG1005';
ft_singleplotTFR(cfg, GA_pre);
title(sprintf('PRE - Comodulogram (%s)', strjoin(electrodes_of_interest, ', ')));
xlabel('EGG phase-driver frequency (Hz)');
ylabel('EEG amplitude frequency (Hz)');

ft_singleplotTFR(cfg, GA_post);
title(sprintf('POST - Comodulogram (%s)', strjoin(electrodes_of_interest, ', ')));
xlabel('EGG phase-driver frequency (Hz)');
ylabel('EEG amplitude frequency (Hz)');

% saveas(gcf, fullfile(dirs.group, 'Comodulogram_PrePost_SelectedElectrodes.png'));
% close(gcf);

% --- Scalp topography for a chosen amp/phase band-of-interest ---
% cfg.xlim selects along the "time" (= phase-frequency) axis,
% cfg.ylim selects along the "freq" (= amplitude-frequency) axis.
cfg = [];
cfg.parameter = 'powspctrm';
cfg.layout = 'EEG1005';
cfg.colorbar  = 'yes';
cfg.xlim      = [.03 .06];
cfg.ylim      = band_of_interest_amp;
cfg.comment   = 'no';
ft_topoplotTFR(cfg, GA_pre);
title('PRE - Topography PAC EGG (.03 - .06) ALPHA (8-13)');

ft_topoplotTFR(cfg, GA_post);
title('POST - Topography PAC EGG (.03 - .06) ALPHA (8-13)');

cfg.ylim      = [3 7];
cfg.comment   = 'no';
ft_topoplotTFR(cfg, GA_pre);
title('PRE - Topography PAC EGG (.03 - .06) THETA (3-7)');

ft_topoplotTFR(cfg, GA_post);
title('POST - Topography PAC EGG (.03 - .06) THETA (3-7)');

cfg.ylim      = [15 30];
cfg.comment   = 'no';
ft_topoplotTFR(cfg, GA_pre);
title('PRE - Topography PAC EGG (.03 - .06) BETA (15-30)');

ft_topoplotTFR(cfg, GA_post);
title('POST - Topography PAC EGG (.03 - .06) BETA (15-30)');

% saveas(gcf, fullfile(dirs.group, 'Comodulogram_PrePost_Topography.png'));
% close(gcf);

%% =========================================================================
%% 5. CLUSTER-BASED (PAIRED, POST vs PRE)
%% =========================================================================
cfg = [];
cfg.method      = 'triangulation';                          
cfg.layout      = 'EEG1005.lay';
cfg.feedback    = 'yes';          % show a neighbour plot
neighbours      = ft_prepare_neighbours(cfg, GA_post); % define neighbouring channels

cfg = [];
cfg.channel       = {'all'};
% cfg.method        = 'analytic';
% cfg.correctm      = 'bonferroni';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.numrandomization = 1000;
cfg.neighbours       = neighbours;

subj = nSubj; % adapt depending on the number of the individual included in your ga data
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

cfg.latency          = [.03 .06];
cfg.avgovertime      = 'yes';
cfg.clusterthreshold = 'nonparametric_common';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 1;
cfg.clustertail      = 1;
cfg.alpha            = 0.05;
cfg.correcttail      = 'alpha';
cfg.frequency        = [8 13];
cfg.avgoverfreq      = 'yes';
stat = ft_freqstatistics(cfg, data_pre{:}, data_post{:});

% save(fullfile(dirs.group, 'Stat_Comodulogram_Paired.mat'), 'stat');
% fprintf('\nCluster stats computed.\n');

%% =========================================================================
%% 6. PLOT SIGNIFICANT CLUSTERS
%% =========================================================================
% (a) T-value comodulogram for the electrodes of interest, with the
%     significance mask outlined.
cfg = [];
cfg.channel       = electrodes_of_interest;
cfg.parameter     = 'stat';
% cfg.maskparameter = 'mask';
% cfg.maskstyle     = 'outline';
cfg.colorbar      = 'yes';

figure('Color','w','Position',[0 0 700 550]);
ft_singleplotTFR(cfg, stat);
xlabel('EGG phase-driver frequency (Hz)');
ylabel('EEG amplitude frequency (Hz)');
title(sprintf('T-values, significant clusters outlined (%s)', ...
      strjoin(electrodes_of_interest, ', ')));
% saveas(gcf, fullfile(dirs.group, 'Stat_Comodulogram_SelectedElectrodes.png'));
% close(gcf);
%% =========================================================================
%% FINAL FIGURE: PAC TOPOGRAPHIES (PRE / POST / CLUSTER STAT) x 3 FREQ BANDS
%% =========================================================================
% Assumes the following already exist in the workspace from the group-level
% script: GA_pre, GA_post, data_pre, data_post, neighbours, design,
% n_permutations, clusteralpha, alpha_level.
%
% Layout: 3 rows (Theta / Alpha / Beta) x 3 columns
%   col 1 = PRE topography
%   col 2 = POST topography          (shares a colorbar/scale with col 1)
%   col 3 = POST vs PRE cluster t-map, significant channels highlighted
%
% Color scaling:
%   - PRE/POST (cols 1-2): MI (Tort modulation index) is non-negative, so a
%     symmetric-around-zero zlim wastes half the colormap. Scale is
%     [0, 98th-percentile] of PRE+POST values in that band's window,
%     rendered with a sequential colormap ('hot').
%   - t-map (col 3): signed statistic -> symmetric zlim, diverging
%     blue-white-red colormap so zero sits at a neutral color.
%
% WHY THIS VERSION IS DIFFERENT: ft_topoplotTFR draws its filled topography
% as scaled/indexed CData, meaning the actual pixel colors are resolved at
% render time from whichever colormap is current -- neither cfg.colormap
% nor a manual colormap(ax,...) call reliably survives once a *later* axis
% in the same figure sets a different colormap. Instead of fighting that,
% each panel's colors are "frozen" into explicit true-color RGB right after
% plotting (freeze_axis_colors, defined at the bottom of this file), which
% makes each axis's rendering permanently independent of every other axis's
% colormap. The axis's Colormap/CLim are still set afterwards purely so its
% colorbar displays the correct gradient/legend.
%
% Layout note: colorbars get their OWN dedicated narrow slots (manually
% positioned axes, not tiledlayout -- ColumnWidth requires MATLAB R2022a+,
% so plain axes('Position', ...) is used for full version independence),
% and each colorbar's Colormap is set explicitly -- MATLAB colorbars
% otherwise default to the FIGURE's colormap rather than their target
% axis's.

%% ---- 1. Config ----
band_names = {'Theta (3-7 Hz)', 'Alpha (8-13 Hz)', 'Beta (15-30 Hz)'};
band_freqs = {[3 7], [8 13], [15 30]};
phase_win  = [.03 .06];   % EGG phase-driver window, fixed across bands

mi_colormap = hot(256);
diverging_map = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1); ...
                  ones(128,1), linspace(1,0,128)', linspace(1,0,128)'];  % blue-white-red

n_bands = numel(band_names);

%% ---- 2. Run cluster stats per band (reusing neighbours/design from before) ----
stat_band = cell(1, n_bands);

for b = 1:n_bands
    cfg = [];
    cfg.channel          = {'all'};
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'ft_statfun_depsamplesT';
    cfg.correctm         = 'cluster';
    cfg.numrandomization = n_permutations;
    cfg.neighbours       = neighbours;

    cfg.design = design;   % paired design built in the group-level script
    cfg.uvar   = 1;
    cfg.ivar   = 2;

    cfg.latency          = phase_win;
    cfg.avgovertime      = 'yes';
    cfg.clusterthreshold = 'nonparametric_common';
    cfg.clusteralpha     = clusteralpha;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan        = 2;
    cfg.tail             = 1;
    cfg.clustertail      = 1;
    cfg.correcttail      = 'alpha';
    cfg.alpha            = alpha_level;
    cfg.frequency        = band_freqs{b};
    cfg.avgoverfreq      = 'yes';

    stat_band{b} = ft_freqstatistics(cfg, data_pre{:}, data_post{:});
end

%% ---- 3. Pre-compute color scales (per band, per column type) ----
mi_zlim  = cell(1, n_bands);   % non-negative, shared PRE/POST scale, per band
t_zlim   = cell(1, n_bands);   % symmetric t-map scale, per band

for b = 1:n_bands
    freq_idx = GA_pre.freq >= band_freqs{b}(1) & GA_pre.freq <= band_freqs{b}(2);
    time_idx = GA_pre.time >= phase_win(1)      & GA_pre.time <= phase_win(2);

    vals_pre  = GA_pre.powspctrm(:, freq_idx, time_idx);
    vals_post = GA_post.powspctrm(:, freq_idx, time_idx);
    pooled    = [vals_pre(:); vals_post(:)];
    pooled    = pooled(~isnan(pooled));
    cap       = prctile(pooled, 98);          % robust upper bound (avoids outlier washout)
    mi_zlim{b} = [0.5e-3, cap];                % lower bound raised to 0.6e-3 (was 0)

    tvals = stat_band{b}.stat(:);
    mt = max(abs(tvals), [], 'omitnan');
    t_zlim{b} = [-mt, mt];
end

%% ---- 4. Build the figure ----
% 5 "columns" per row: PRE | POST | MI colorbar | t-map | t colorbar.
% Positions are computed manually (normalized figure units) instead of via
% tiledlayout, so this works on any MATLAB version and topography panels
% are guaranteed identical size across all rows -- colorbars get their own
% narrow slots and can never distort a topography tile's width.
fig = figure('Color', 'w', 'Position', [50 50 1350 1250]);

margin_l = 0.05; margin_r = 0.03; margin_t = 0.14; margin_b = 0.05;
gap_col  = 0.015; gap_row = 0.09;
col_rel  = [1, 1, 0.12, 1, 0.12];              % PRE | POST | cb | t-map | cb
n_col    = numel(col_rel);

avail_w   = 1 - margin_l - margin_r - (n_col-1)*gap_col;
unit_w    = avail_w / sum(col_rel);
col_w     = unit_w * col_rel;
col_x     = margin_l + [0, cumsum(col_w(1:end-1) + gap_col)];

avail_h  = 1 - margin_t - margin_b - (n_bands-1)*gap_row;
row_h    = avail_h / n_bands;
row_y    = @(r) 1 - margin_t - r*row_h - (r-1)*gap_row;   % r = 1 is the top row

for b = 1:n_bands

    band_lbl = band_names{b};
    freq_rng = band_freqs{b};
    y        = row_y(b);

    % ---- PRE topography ----
    ax1 = axes('Parent', fig, 'Position', [col_x(1), y, col_w(1), row_h]);
    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.layout    = 'EEG1005';
    cfg.xlim      = phase_win;
    cfg.ylim      = freq_rng;
    cfg.zlim      = mi_zlim{b};
    cfg.comment   = 'no';
    cfg.colorbar  = 'no';
    cfg.figure    = ax1;
    ft_topoplotTFR(cfg, GA_pre);
    freeze_axis_colors(ax1, mi_colormap, mi_zlim{b});
    title(sprintf('PRE\n%s', band_lbl), 'FontWeight', 'normal', 'Color', 'k');

    % ---- POST topography (same scale/colormap as PRE -> directly comparable) ----
    ax2 = axes('Parent', fig, 'Position', [col_x(2), y, col_w(2), row_h]);
    cfg.figure = ax2;
    ft_topoplotTFR(cfg, GA_post);
    freeze_axis_colors(ax2, mi_colormap, mi_zlim{b});
    title(sprintf('POST\n%s', band_lbl), 'FontWeight', 'normal', 'Color', 'k');

    % MI colorbar: own dedicated slot, own explicit Colormap (colorbar
    % objects otherwise default to the FIGURE's colormap, not the axis's)
    cb1 = colorbar(ax2);
    cb1.Position = [col_x(3), y + 0.1*row_h, col_w(3), 0.8*row_h];
    cb1.Colormap = mi_colormap;
    cb1.Limits   = mi_zlim{b};
    ylabel(cb1, 'MI (a.u.)');

    % ---- Cluster stat: POST vs PRE (t-values, significant channels marked) ----
    ax3 = axes('Parent', fig, 'Position', [col_x(4), y, col_w(4), row_h]);
    stat = stat_band{b};

    sig_mask = false(numel(stat.label), 1);
    if isfield(stat, 'posclusters') && ~isempty(stat.posclusters)
        sig_clusters = find([stat.posclusters.prob] < alpha_level);
        sig_mask = ismember(stat.posclusterslabelmat, sig_clusters);
    end

    cfg = [];
    cfg.parameter        = 'stat';
    cfg.layout            = 'EEG1005';
    cfg.zlim              = t_zlim{b};
    cfg.comment           = 'no';
    cfg.colorbar          = 'no';
    cfg.highlight          = 'on';
    cfg.highlightchannel   = stat.label(sig_mask);
    cfg.highlightsymbol    = 'o';
    cfg.highlightsize      = 8;
    cfg.figure             = ax3;
    ft_topoplotTFR(cfg, stat);
    freeze_axis_colors(ax3, diverging_map, t_zlim{b});
    title(sprintf('POST vs PRE (t)\n%s', band_lbl), 'FontWeight', 'normal', 'Color', 'k');

    cb2 = colorbar(ax3);
    cb2.Position = [col_x(5), y + 0.1*row_h, col_w(5), 0.8*row_h];
    cb2.Colormap = diverging_map;
    cb2.Limits   = t_zlim{b};
    ylabel(cb2, 't-value');
end

sgtitle(sprintf('PAC (EEG amplitude x EGG phase %.3g-%.3g Hz): PRE vs POST', ...
        phase_win(1), phase_win(2)), 'FontWeight', 'bold', 'Color', 'k');

% To save:
% exportgraphics(fig, fullfile(dirs.group, 'Figure_PrePost_Cluster_3Bands.png'), 'Resolution', 300);
% exportgraphics(fig, fullfile(dirs.group, 'Figure_PrePost_Cluster_3Bands.pdf'), 'ContentType', 'vector');


%% =========================================================================
%% LOCAL FUNCTION: freeze_axis_colors
%% =========================================================================
function freeze_axis_colors(ax, cmap, clim)
% Bake every colormapped graphics object in AX into explicit true-color RGB,
% using CMAP and CLIM, so its appearance no longer depends on which
% colormap happens to be current elsewhere in the figure. Also sets the
% axis's own Colormap/CLim so a colorbar(ax) attached to it still shows the
% correct gradient.

    objs = findobj(ax, '-property', 'CData');
    for k = 1:numel(objs)
        obj = objs(k);
        cdata = get(obj, 'CData');

        % skip objects that are already true-color (RGB triplets, 3rd dim = 3)
        if ndims(cdata) == 3 && size(cdata, 3) == 3
            continue
        end

        idx = (cdata - clim(1)) / max(clim(2) - clim(1), eps);   % 0-1
        idx = round(idx * (size(cmap, 1) - 1)) + 1;
        idx = min(max(idx, 1), size(cmap, 1));

        nanmask = isnan(cdata);
        rgb = reshape(cmap(idx(:), :), [size(cdata), 3]);
        if any(nanmask(:))
            for ch = 1:3
                chan = rgb(:, :, ch);
                chan(nanmask) = NaN;
                rgb(:, :, ch) = chan;
            end
        end

        set(obj, 'CData', rgb, 'CDataMapping', 'direct');
    end

    ax.Colormap = cmap;
    ax.CLim     = clim;
end

%% =========================================================================
%% EXTRACT PER-SUBJECT PAC (MI) VALUES, AVERAGED OVER ELECTRODE ROIs PER BAND
%% =========================================================================
% Assumes the following already exist in the workspace (from the
% group-level loading script): data_pre, data_post, used_subjects.
%
% For each subject and condition (PRE/POST), averages MI over:
%   - the fixed EGG phase-driver window (.03-.06 Hz)
%   - the band's amplitude-frequency range
%   - the band's specific electrode ROI (as given, not a scalp-wide average)
%
% Output: one row per subject, one column per band x condition
% (Theta_PRE, Theta_POST, Alpha_PRE, Alpha_POST, Beta_PRE, Beta_POST),
% written to PAC_ROI_values.csv.

%% ---- 1. Config: bands, phase window, and electrode ROIs ----
band_names = {'Theta', 'Alpha', 'Beta'};
band_freqs = {[3 7], [8 13], [15 30]};
phase_win  = [.03 .06];

roi_electrodes = { ...
    {'Fp1','F7','F8','Fp2','AF7','F5','F6','AF8'}, ...                                  % Theta ROI
    {'C3','CP1','Pz','P3','P4','CP6','CP2','Cz','C4','C1','CP3','P1','P5','P6','P2','CPz','CP4','C6','C2'}, ... % Alpha ROI
    {'FC5','C3','C4','FC6','FC3','C5','C6','FC4'} ...                                   % Beta ROI
};

n_bands = numel(band_names);
n_subj  = numel(used_subjects);

%% ---- 2. Extract, per subject / band / condition ----
PAC_pre  = nan(n_subj, n_bands);
PAC_post = nan(n_subj, n_bands);

for b = 1:n_bands

    freq_rng = band_freqs{b};
    roi      = roi_electrodes{b};

    for i = 1:n_subj

        % --- PRE ---
        d = data_pre{i};
        freq_idx = d.freq >= freq_rng(1) & d.freq <= freq_rng(2);
        time_idx = d.time >= phase_win(1) & d.time <= phase_win(2);
        chan_idx = ismember(lower(d.label), lower(roi));

        missing = setdiff(lower(roi), lower(d.label));
        if ~isempty(missing)
            fprintf('[!] %s PRE band %s: channel(s) not found: %s\n', ...
                    used_subjects{i}, band_names{b}, strjoin(missing, ', '));
        end

        vals = d.powspctrm(chan_idx, freq_idx, time_idx);
        PAC_pre(i, b) = mean(vals(:), 'omitnan');

        % --- POST ---
        d = data_post{i};
        freq_idx = d.freq >= freq_rng(1) & d.freq <= freq_rng(2);
        time_idx = d.time >= phase_win(1) & d.time <= phase_win(2);
        chan_idx = ismember(lower(d.label), lower(roi));

        missing = setdiff(lower(roi), lower(d.label));
        if ~isempty(missing)
            fprintf('[!] %s POST band %s: channel(s) not found: %s\n', ...
                    used_subjects{i}, band_names{b}, strjoin(missing, ', '));
        end

        vals = d.powspctrm(chan_idx, freq_idx, time_idx);
        PAC_post(i, b) = mean(vals(:), 'omitnan');
    end
end

%% ---- 3. Assemble into a table and save ----
col_names = {};
col_data  = {};
for b = 1:n_bands
    col_names{end+1} = sprintf('%s_PRE',  band_names{b}); %#ok<SAGROW>
    col_data{end+1}  = PAC_pre(:, b);
    col_names{end+1} = sprintf('%s_POST', band_names{b}); %#ok<SAGROW>
    col_data{end+1}  = PAC_post(:, b);
end

T = table(used_subjects(:), 'VariableNames', {'Subject'});
for k = 1:numel(col_names)
    T.(col_names{k}) = col_data{k};
end

disp(T);

out_path = fullfile(pwd, 'PAC_ROI_values.csv');
writetable(T, out_path);
fprintf('\nSaved: %s\n', out_path);