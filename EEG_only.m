ft_defaults
addpath('/Users/qmoreau/Documents/Work/UGO_PAC/eeglab2025.1.0');   % le dossier racine d'EEGLAB, pas un sous-dossier
eeglab nogui

path_cleaned = '/Users/qmoreau/Documents/Work/UGO_PAC/EEG_preproc_ASR_100';
path_freq    = '/Users/qmoreau/Documents/Work/UGO_PAC/PAC_PreVsPost_Comparison/Group_Stats';

subjects = {'P001','P002', 'P003', 'P004', 'P005', 'P006', 'P007', 'P008', 'P009', 'P010', ...
    'P011', 'P012', 'P013', 'P014', 'P015', 'P017', 'P018', 'P019', 'P021', ...
    'P022', 'P023', 'P024', 'P025', 'P026', 'P027', 'P028', 'P029', 'P030', ...
    'P031', 'P032','P033', 'P034', 'P035', 'P036', 'P037', 'P038', 'P039','P040', 'P0166', 'P0200'};


conditions = {'pre', 'post'};

for i = 1:length(subjects)
    subj = subjects{i};
    fprintf('\n=== SOGGETTO %d/%d: %s ===\n', i, length(subjects), subj);

    for c = 1:length(conditions)
        cond_name = conditions{c};

        % --- A. FILE FINDING ---
        file_eeg = sprintf('%s_%s_cleaned_UGO.set', subj, cond_name);
        full_eeg = fullfile(path_cleaned, file_eeg);

        try
            if ~exist(full_eeg, 'file')
                error('File not found: %s', full_eeg);
            end

            EEG = pop_loadset('filename', file_eeg, 'filepath', path_cleaned);

            %  CONVERT TO FIELDTRIP FORMAT 
            data_ft = eeglab2fieldtrip(EEG, 'raw', 'none');

            cfg = [];
            cfg.overlap = 0;
            cfg.length  = 4;
            data_ft_seg    = ft_redefinetrial(cfg, data_ft);
            
            % Reset every trial's time axis to start at 0
            ntrials = numel(data_ft_seg.time);
            for k = 1:ntrials
                nsamples = numel(data_ft_seg.time{k});
                data_ft_seg.time{k} = (0:nsamples-1) / data_ft_seg.fsample;
            end

            % TIME-FREQUENCY ANALYSIS 
            cfg            = [];
            cfg.method     = 'mtmconvol';
            cfg.taper      = 'hanning';
            cfg.foi        = 1:.5:40;                 % frequency band of interest
            cfg.channel    = 'all';
            cfg.t_ftimwin  = ones(length(cfg.foi),1)*0.5; % 50 ms time window
            cfg.output     = 'pow';
            cfg.keeptrials  = 'no';
            cfg.toi        = 0:0.5:4;  
            freq           = ft_freqanalysis(cfg, data_ft_seg);

            % disp(freq)
            % 
            % cfg = [];
            % cfg.baselinetype = 'db';
            % cfg.baseline = [-inf inf];
            % cfg.layout = 'EEG1005';
            % cfg.showscale = 'no';
            % cfg.showcmnt = 'no';
            % ft_multiplotTFR(cfg, freq);

            % SAVE 
            out_file = fullfile(path_freq, sprintf('%s_%s_freq.mat', subj, cond_name));
            save(out_file, 'freq', '-v7.3');
            fprintf('   Saved %s\n', out_file);

        catch ME
            fprintf('   [!] Errore %s: %s\n', cond_name, ME.message);
        end
    end
end

%% cluster based permutation test

subjects = {'P001','P002', 'P003', 'P004', 'P005', 'P006', 'P007', 'P008', 'P009', 'P010', ...
    'P011', 'P012', 'P013', 'P014', 'P015', 'P017', 'P018', 'P019', 'P021', ...
    'P022', 'P023', 'P024', 'P025', 'P026', 'P027', 'P028', 'P029', 'P030', ...
    'P031', 'P032','P033', 'P034', 'P035', 'P036', 'P037', 'P038', 'P039','P040', 'P0166', 'P0200'};


%% --- A. LOAD ALL SUBJECTS, KEEP ONLY COMPLETE PRE/POST PAIRS ---
freq_pre  = {};
freq_post = {};
valid_subjects = {};

for i = 1:length(subjects)
    subj = subjects{i};
    f_pre  = fullfile(path_freq, sprintf('%s_pre_freq.mat',  subj));
    f_post = fullfile(path_freq, sprintf('%s_post_freq.mat', subj));

    if exist(f_pre, 'file') && exist(f_post, 'file')
        tmp_pre  = load(f_pre,  'freq');
        tmp_post = load(f_post, 'freq');
        valid_subjects{end+1} = subj;     
        freq_pre{end+1}  = tmp_pre.freq; 
        freq_post{end+1} = tmp_post.freq; 
    else
        fprintf('   [!] Skipping %s: missing pre and/or post freq file\n', subj);
    end
end

nsubj = numel(valid_subjects);
fprintf('\nRunning cluster stats on %d subjects with complete pre/post data\n', nsubj);
if nsubj < 2
    error('Not enough subjects with both pre and post data to run stats.');
end

cfg = [];
cfg.keepindividual = 'yes';
% cfg.foilim         = 'all';
% cfg.toilim         = 'all';
GA_pre  = ft_freqgrandaverage(cfg, freq_pre{:});
GA_post = ft_freqgrandaverage(cfg, freq_post{:});

cfg = [];
cfg.method      = 'triangulation';                          
cfg.layout      = 'EEG1005.lay';
cfg.feedback    = 'no';          % show a neighbour plot
neighbours      = ft_prepare_neighbours(cfg, GA_pre); % define neighbouring channels

design = zeros(2, 2*nsubj);
design(1, 1:nsubj)          = 1:nsubj;   % uvar: subject identity
design(1, nsubj+1:2*nsubj)  = 1:nsubj;
design(2, 1:nsubj)          = 1;         % ivar: condition (1 = pre)
design(2, nsubj+1:2*nsubj)  = 2;         %              (2 = post)

bands = struct('name', {'theta','alpha','beta'}, 'range', {[3 7], [8 13], [15 30]});
stat  = struct();
time_win = [0 4];

for b = 1:numel(bands)
    bandname = bands(b).name;
    fprintf('\n=== Cluster-based permutation test: %s (%g-%g Hz) ===\n', ...
        bandname, bands(b).range(1), bands(b).range(2));

    % --- Manually pre-average each subject's data to a single value per
    % channel (validated to give correct, non-NaN results -- see the
    % Direct_Statfun_Check.m diagnostic). ft_freqstatistics's own internal
    % avgovertime/avgoverfreq averaging was silently producing NaN for
    % every channel despite clean input data; pre-averaging ourselves and
    % feeding FieldTrip already-collapsed (singleton freq/time) data
    % sidesteps that step entirely while still using its cluster-based
    % permutation machinery (neighbours, montecarlo shuffling) unchanged. ---
    band_pre  = cell(1, nsubj);
    band_post = cell(1, nsubj);

    for i = 1:nsubj
        f = freq_pre{i};
        freq_idx = f.freq >= bands(b).range(1) & f.freq <= bands(b).range(2);
        time_idx = f.time >= time_win(1)        & f.time <= time_win(2);

        band_pre{i}.label     = f.label;
        band_pre{i}.dimord    = 'chan_freq_time';
        band_pre{i}.freq      = mean(bands(b).range);
        band_pre{i}.time      = mean(time_win);
        band_pre{i}.powspctrm = mean(f.powspctrm(:, freq_idx, time_idx), [2 3], 'omitnan');

        f = freq_post{i};
        band_post{i}.label     = f.label;
        band_post{i}.dimord    = 'chan_freq_time';
        band_post{i}.freq      = mean(bands(b).range);
        band_post{i}.time      = mean(time_win);
        band_post{i}.powspctrm = mean(f.powspctrm(:, freq_idx, time_idx), [2 3], 'omitnan');
    end

    cfg                  = [];
    cfg.channel          = 'all';
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'ft_statfun_depsamplesT';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.025;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan        = 2;
    cfg.tail             = 0;
    cfg.clustertail      = 0;
    cfg.alpha            = 0.025;            % two-tailed -> 0.025 each side
    cfg.numrandomization = 5000;
    cfg.neighbours       = neighbours;
    cfg.design           = design;
    cfg.uvar             = 1;
    cfg.ivar             = 2;
    % NOTE: no cfg.latency/cfg.frequency/cfg.avgovertime/cfg.avgoverfreq --
    % data is already collapsed to one value per channel, above.

    stat.(bandname) = ft_freqstatistics(cfg, band_pre{:}, band_post{:});

    sig = stat.(bandname).mask;
    fprintf('   Significant (channel,time) points: %d / %d\n', sum(sig(:)), numel(sig));
    fprintf('   NaN t-values: %d / %d\n', sum(isnan(stat.(bandname).stat(:))), numel(stat.(bandname).stat));
end

%% =========================================================================
%% FINAL FIGURE: SPECTRAL POWER (PRE / POST / CLUSTER STAT) x 3 FREQ BANDS
%% =========================================================================
% Assumes the following already exist in the workspace, from the
% preprocessing + cluster-stat code above:
%   GA_pre, GA_post   - ft_freqgrandaverage output, keepindividual='yes'
%                        (dimord 'subj_chan_freq_time') -- used here only
%                        to build the PRE/POST topographies, not the stats
%   stat              - struct with fields .theta / .alpha / .beta, each
%                        the ft_freqstatistics output for that band
%                        (two-tailed, cfg.latency=[0 4], cfg.avgovertime/
%                        avgoverfreq = 'yes')
%   valid_subjects, neighbours, design  - as built above
%
% Layout: 3 rows (Theta / Alpha / Beta) x 2 "blocks":
%   col 1 = PRE topography
%   col 2 = POST topography          (shares a colorbar/scale with col 1)
%   col 3 = POST vs PRE cluster t-map, significant channels highlighted
%   (+ 2 narrow dedicated colorbar columns, one per block)
%
% This mirrors Figure_PrePost_Cluster_3Bands.m (the PAC figure) and reuses
% the same fixes worked out there:
%   - power is non-negative -> sequential 'hot' colormap, zlim anchored at 0
%     (not symmetric), t-map keeps a diverging blue-white-red colormap
%   - each panel's colors are "frozen" into explicit RGB right after
%     plotting (freeze_axis_colors, bottom of file) so FieldTrip's
%     figure-wide colormap behavior can't make one band's colors bleed
%     into another
%   - colorbars get their own dedicated narrow axes (not tiledlayout, since
%     ColumnWidth needs MATLAB R2022a+) so all three rows stay the same size
%   - all titles forced to black, extra top/row margin to avoid overlap
%   - any channel with a NaN t-value is dropped from the t-map topoplot
%     (FieldTrip's internal NaN handling can crash griddata otherwise); if
%     fewer than 5 channels remain valid, the topoplot is skipped and a
%     placeholder message shown instead of letting the whole figure crash
%
% Significance: stat.(band).mask is used directly (FieldTrip already
% combines both tails of the two-tailed test into this field), so no
% manual posclusters/negclusters bookkeeping is needed here.
 
%% ---- 1. Config ----
band_names = {'Theta (3-7 Hz)', 'Alpha (8-13 Hz)', 'Beta (15-30 Hz)'};
band_keys  = {'theta', 'alpha', 'beta'};
band_freqs = {[3 7], [8 13], [15 30]};
time_win   = [0 4];        % full epoch, matches cfg.latency used for the stats
 
pow_colormap  = hot(256);
diverging_map = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1); ...
                  ones(128,1), linspace(1,0,128)', linspace(1,0,128)'];  % blue-white-red
 
n_bands = numel(band_names);
 
%% ---- 2. Collapse subject dimension for plotting (stats use freq_pre/freq_post directly, see fix above) ----
avg_pre  = GA_pre;
avg_pre.powspctrm = squeeze(mean(GA_pre.powspctrm, 1, 'omitnan'));
avg_pre.dimord    = 'chan_freq_time';
 
avg_post = GA_post;
avg_post.powspctrm = squeeze(mean(GA_post.powspctrm, 1, 'omitnan'));
avg_post.dimord    = 'chan_freq_time';
 
%% ---- 3. Pre-compute color scales (per band) ----
pow_zlim = cell(1, n_bands);   % non-negative, shared PRE/POST scale, per band
t_zlim   = cell(1, n_bands);   % symmetric t-map scale, per band
 
for b = 1:n_bands
    freq_idx = avg_pre.freq >= band_freqs{b}(1) & avg_pre.freq <= band_freqs{b}(2);
    time_idx = avg_pre.time >= time_win(1)      & avg_pre.time <= time_win(2);
 
    vals_pre  = avg_pre.powspctrm(:, freq_idx, time_idx);
    vals_post = avg_post.powspctrm(:, freq_idx, time_idx);
    pooled    = [vals_pre(:); vals_post(:)];
    pooled    = pooled(~isnan(pooled));
    cap       = prctile(pooled, 98);
    pow_zlim{b} = [0, cap];    % power is non-negative -> anchor at 0
 
    tvals = stat.(band_keys{b}).stat(:);
    mt = max(abs(tvals), [], 'omitnan');
    t_zlim{b} = [-mt, mt];
end
 
%% ---- 4. Build the figure ----
fig = figure('Color', 'w', 'Position', [50 50 1350 1250]);
 
margin_l = 0.05; margin_r = 0.03; margin_t = 0.14; margin_b = 0.05;
gap_col  = 0.015; gap_row = 0.09;
col_rel  = [1, 1, 0.12, 1, 0.12];              % PRE | POST | cb | t-map | cb
n_col    = numel(col_rel);
 
avail_w = 1 - margin_l - margin_r - (n_col-1)*gap_col;
unit_w  = avail_w / sum(col_rel);
col_w   = unit_w * col_rel;
col_x   = margin_l + [0, cumsum(col_w(1:end-1) + gap_col)];
 
avail_h = 1 - margin_t - margin_b - (n_bands-1)*gap_row;
row_h   = avail_h / n_bands;
row_y   = @(r) 1 - margin_t - r*row_h - (r-1)*gap_row;
 
for b = 1:n_bands
 
    band_lbl = band_names{b};
    freq_rng = band_freqs{b};
    bkey     = band_keys{b};
    y        = row_y(b);
 
    % ---- PRE topography ----
    ax1 = axes('Parent', fig, 'Position', [col_x(1), y, col_w(1), row_h]);
    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.layout    = 'EEG1005';
    cfg.xlim      = time_win;
    cfg.ylim      = freq_rng;
    cfg.zlim      = pow_zlim{b};
    cfg.comment   = 'no';
    cfg.colorbar  = 'no';
    cfg.figure    = ax1;
    ft_topoplotTFR(cfg, avg_pre);
    freeze_axis_colors(ax1, pow_colormap, pow_zlim{b});
    title(sprintf('PRE\n%s', band_lbl), 'FontWeight', 'normal', 'Color', 'k');
 
    % ---- POST topography (same scale/colormap as PRE -> directly comparable) ----
    ax2 = axes('Parent', fig, 'Position', [col_x(2), y, col_w(2), row_h]);
    cfg.figure = ax2;
    ft_topoplotTFR(cfg, avg_post);
    freeze_axis_colors(ax2, pow_colormap, pow_zlim{b});
    title(sprintf('POST\n%s', band_lbl), 'FontWeight', 'normal', 'Color', 'k');
 
    % power colorbar: own dedicated slot, own explicit Colormap
    cb1 = colorbar(ax2);
    cb1.Color = 'k';
    cb1.Position = [col_x(3), y + 0.1*row_h, col_w(3), 0.8*row_h];
    cb1.Colormap = pow_colormap;
    cb1.Limits   = pow_zlim{b};
    ylabel(cb1, 'Power (a.u.)', 'Color', 'k');
 
    % ---- Cluster stat: POST vs PRE (t-values, significant channels marked) ----
    ax3 = axes('Parent', fig, 'Position', [col_x(4), y, col_w(4), row_h]);
    this_stat = stat.(bkey);
 
    sig_mask  = squeeze(this_stat.mask);   % combines both tails already
    valid_chan_mask = ~isnan(squeeze(this_stat.stat));
 
    if any(~valid_chan_mask)
        fprintf('[!] %s: excluding %d/%d channel(s) with NaN t-values from the topoplot: %s\n', ...
                band_lbl, nnz(~valid_chan_mask), numel(valid_chan_mask), ...
                strjoin(this_stat.label(~valid_chan_mask), ', '));
    end
 
    if nnz(valid_chan_mask) < 5
        % Too few valid channels for griddata to interpolate a scalp map
        % (a known MATLAB bug crashes it below ~3 points) -- show a
        % placeholder instead of taking down the whole figure.
        axis(ax3, 'off');
        text(ax3, 0.5, 0.5, sprintf('Only %d/%d valid channels\n-- topoplot skipped', ...
             nnz(valid_chan_mask), numel(valid_chan_mask)), ...
             'HorizontalAlignment', 'center', 'Color', 'k', 'Parent', ax3);
        title(ax3, sprintf('POST vs PRE (t)\n%s', band_lbl), 'FontWeight', 'normal', 'Color', 'k');
        continue
    end
 
    cfg = [];
    cfg.parameter        = 'stat';
    cfg.layout            = 'EEG1005';
    cfg.channel            = this_stat.label(valid_chan_mask);   % drop NaN channels (griddata crashes otherwise)
    cfg.zlim              = t_zlim{b};
    cfg.comment           = 'no';
    cfg.colorbar          = 'no';
    cfg.highlight          = 'on';
    cfg.highlightchannel   = this_stat.label(sig_mask & valid_chan_mask);
    cfg.highlightsymbol    = '.';
    cfg.highlightcolor     = 'r';
    cfg.highlightsize      = 20;
    cfg.figure             = ax3;
    ft_topoplotTFR(cfg, this_stat);
    freeze_axis_colors(ax3, diverging_map, t_zlim{b});
    title(sprintf('POST vs PRE (t)\n%s', band_lbl), 'FontWeight', 'normal', 'Color', 'k');
 
    cb2 = colorbar(ax3);
    cb2.Color = 'k';
    cb2.Position = [col_x(5), y + 0.1*row_h, col_w(5), 0.8*row_h];
    cb2.Colormap = diverging_map;
    cb2.Limits   = t_zlim{b};
    ylabel(cb2, 't-value');
end
 
sgtitle('EEG power: PRE vs POST', 'FontWeight', 'bold', 'Color', 'k');
 
% To save:
% exportgraphics(fig, fullfile(path_freq, 'Figure_Power_PrePost_Cluster_3Bands.png'), 'Resolution', 300);
% exportgraphics(fig, fullfile(path_freq, 'Figure_Power_PrePost_Cluster_3Bands.pdf'), 'ContentType', 'vector');
 
 
%% =========================================================================
%% LOCAL FUNCTION: freeze_axis_colors (same as in the PAC figure scripts)
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
 
        if ndims(cdata) == 3 && size(cdata, 3) == 3
            continue
        end
 
        idx = (cdata - clim(1)) / max(clim(2) - clim(1), eps);
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
 