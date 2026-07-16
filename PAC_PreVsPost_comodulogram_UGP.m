%%
% =========================================================================
% COMODULOGRAM PIPELINE UGP: EGG (PHASE) -> EEG (AMPLITUDE) PAC, PRE vs POST
% =========================================================================
% For every subject and every condition (pre / post), this script:
%   1. Loads and cleans the raw EGG signal, selects the best gastric
%   channel (AI pipeline)
%   2. Loads the cleaned EEG and builds FieldTrip raw structures
%   3. Bandpass-filters the EEG into a grid of amplitude bands (0.5-40 Hz)
%      and extracts the amplitude envelope of each, per channel
%   4. Bandpass-filters the EGG into a grid of phase (driver) bands
%      (0.02-0.10 Hz) and extracts the instantaneous phase of each
%   5. Computes a phase-amplitude coupling comodulogram (Tort's Modulation
%      Index by default) for every channel x amplitude-band x phase-band
%      combination
%   6. Saves the comodulogram + a diagnostic figure per subject/condition

eeglab nogui;
% ft_defaults;

% =========================================================================
% 1. CONFIGURAZIONE
% =========================================================================

path_cleaned = '/Users/qmoreau/Documents/Work/UGO_PAC/PAC_tryouts_UGP/EEG_preproc_ASR_100';

output_root  = '/Users/qmoreau/Documents/Work/UGO_PAC/PAC_PreVsPost_Comparison';

% Raw EGG/BrainVision recordings (was path_raw_root, referenced later but
% never defined in the original script)
path_raw_root = '/Users/qmoreau/Documents/Work/EGG_Giusi/Data_Pill';

dirs = []; dirs.root = output_root;
dirs.figs = fullfile(output_root, 'Diagnostic_Figs');
dirs.stats = fullfile(output_root, 'Group_Stats');
fn = fieldnames(dirs);
for k = 1:length(fn), if ~exist(dirs.(fn{k}), 'dir'), mkdir(dirs.(fn{k})); end; end

% Parametri
% Target sample rate the cleaned EEG is resampled to before filtering/PAC.
work_srate_eeg = 125;
% alpha_band     = [8 13];      % kept only for reference/legacy, unused by comodulogram
% Common sample rate (Hz) that both the EGG phase series and the EEG
% amplitude-envelope series are resampled/interpolated to, so phase and
% amplitude time series can be directly compared sample-by-sample.
pac_calc_srate = 20;          % common analysis rate for phase & amplitude series
% Expected gastric slow-wave frequency range (Hz), ~2-4 cycles/min, used
% only to pick the best of the 4 EGG channels (Section B1).
egg_freq_range = [0.033 0.067];
% EGG signal is downsampled to this rate before anomaly detection/cleaning.
fs_egg_proc    = 10;

n_bins_MI      = 18;          % Tort MI histogram bins
edge_trim_sec  = 5;           % seconds trimmed from each end after filtering (edge artifacts)
pac_method     = 'tort';      % 'tort' | 'ozkurt' | 'canolty' | 'plv'

% --- Comodulogram frequency grids ---
% X-axis: EGG driver (phase) frequency, Hz. 
% Sliding-window grid of narrow phase (driver) bands spanning the gastric
% range and a bit beyond it: each band is phase_bw Hz wide, and
% consecutive bands are spaced phase_step Hz apart (overlapping windows).
phase_lo = 0.02;  phase_hi = 0.10;  phase_bw = 0.01;  phase_step = 0.005;
% Y-axis: EEG (amplitude) frequency, Hz, 0.5-40 
% Sliding-window grid of amplitude bands spanning essentially the whole
% EEG spectrum (delta through gamma): each band is amp_bw Hz wide, spaced
% amp_step Hz apart.
amp_lo   = 0.5;   amp_hi   = 40;    amp_bw   = 2;     amp_step   = 1;

% Build the actual [start, end] band edges and center frequencies for the
% phase-driver grid from the parameters above.
phase_starts = phase_lo : phase_step : (phase_hi - phase_bw);
phase_bands  = [phase_starts', phase_starts' + phase_bw];
phase_centers = mean(phase_bands, 2);

% Same, for the amplitude grid.
amp_starts = amp_lo : amp_step : (amp_hi - amp_bw);
amp_bands  = [amp_starts', amp_starts' + amp_bw];
amp_centers = mean(amp_bands, 2);

% Number of bands in each grid - defines the size of the comodulogram.
n_phase = size(phase_bands,1);
n_amp   = size(amp_bands,1);

% LISTA SOGGETTI COMPLETA
subjects = {'P001','P002', 'P003', 'P004', 'P005', 'P006', 'P007', 'P008', 'P009', 'P010', ...
            'P011', 'P012', 'P013', 'P014', 'P015', 'P017', 'P018', 'P019', 'P021', ...
            'P022', 'P023', 'P024', 'P025', 'P026', 'P027', 'P028', 'P029', 'P030', ...
            'P031', 'P032','P033', 'P034', 'P035', 'P036', 'P037', 'P038', 'P039','P040', 'P0166', 'P0200'};

% subjects = {'P002'};
addpath('/Users/qmoreau/Documents/Work/Matlab_EGG/codici_PAC')
addpath('/Users/qmoreau/Documents/Work/Matlab_EGG')
% toolbox with functions to compute PAC -> https://neurofractal.github.io/PACmeg/
addpath(genpath('/Users/qmoreau/Documents/Work/UGO_PAC/PACmeg'));


% =========================================================================
% 2. MAIN LOOP (Pre vs Post, per subject)
% =========================================================================
fprintf('Avvio Analisi Comodulogramma Pre vs Post...\n');

for i = 1:length(subjects)
    subj = subjects{i};
    fprintf('\n=== SOGGETTO %d/%d: %s ===\n', i, length(subjects), subj);
    % Inner loop: process the pre and post condition separately for this
    % subject. Each iteration is fully independent (its own file lookup,
    % its own try/catch), so a failure in one condition doesn't affect
    % the other.
    conditions = {'pre', 'post'};

    for c = 1:2
        cond_name = conditions{c};

        % --- A. FILE FINDING ---
        file_eeg = sprintf('%s_%s_cleaned_UGO.set', subj, cond_name);
        full_eeg = fullfile(path_cleaned, file_eeg);

        if ~exist(full_eeg, 'file')
            fprintf('   [!] EEG %s mancante. Skip.\n', cond_name); continue;
        end

        raw_dir = fullfile(path_raw_root, subj, 'BASELINE');
        if ~exist(raw_dir, 'dir'), raw_dir = fullfile(path_raw_root, subj); end

        if strcmp(cond_name, 'pre')
            d = dir(fullfile(raw_dir, '*PRECOL*.vhdr'));
        else
            d = dir(fullfile(raw_dir, '*ST*.vhdr'));
            if isempty(d), d = dir(fullfile(raw_dir, '*POST*.vhdr')); end
        end

        if isempty(d)
            fprintf('   [!] EGG Raw (%s) non trovato. Skip.\n', cond_name); continue;
        end

        try
            % ----------------------------------------------------------------
            % B1. EGG Processing (unchanged from original pipeline)
            % ----------------------------------------------------------------
            % Load the raw BrainVision recording and keep only the 4 EGG
            % channels (evalc suppresses EEGLAB's verbose loading output).
            evalc('EEG_raw = pop_loadbv(raw_dir, d(1).name);');
            EEG_raw = pop_select(EEG_raw, 'channel', {'egg1', 'egg2', 'egg3', 'egg4'});

            % Downsample the EGG to fs_egg_proc if it was recorded faster
            % (gastric rhythm is very slow, so a low sample rate suffices
            % and keeps subsequent filtering/anomaly detection cheap).
            if EEG_raw.srate > fs_egg_proc, EEG_egg = pop_resample(EEG_raw, fs_egg_proc);
            else, EEG_egg = EEG_raw; end

            fs_egg = EEG_egg.srate;
            raw_data_egg = double(EEG_egg.data);
            t_egg = reshape(EEG_egg.times/1000, 1, []);

            % For each of the 4 EGG channels: interpolate/repair artifacts
            % (find_anomalies_pure), then compute a periodogram to
            % quantify how much power falls inside the expected gastric
            % frequency range (egg_freq_range).
            peaks = zeros(1,4); cleans = zeros(size(raw_data_egg)); cell_anom=cell(1,4);
            all_pxx = {}; all_f = {};
            for ch=1:4
                [filled, anom] = find_anomalies_pure(raw_data_egg(ch,:), t_egg, 1000, fs_egg, 100);
                cleans(ch,:) = reshape(filled,1,[]); cell_anom{ch}=anom;
                [pxx, f] = periodogram(filled, hamming(length(filled)), [], fs_egg);
                idx = f>=egg_freq_range(1) & f<=egg_freq_range(2);
                if any(idx), peaks(ch)=max(pxx(idx)); end
                all_pxx{ch} = pxx; all_f{ch} = f;
            end
            % Select the channel with the strongest gastric-band peak as
            % the "best" EGG channel to use as the phase driver.
            [~, best_ch] = max(peaks);

            egg_clean_data  = cleans(best_ch,:);
            egg_time        = t_egg;
            egg_srate       = fs_egg;
            egg_best_ch     = best_ch;

            % Save the cleaned EGG signal for this subject/condition, so
            % it can be inspected or reused without re-running B1.
            out_file = fullfile(dirs.stats, sprintf('%s_%s_EGGclean.mat', subj, cond_name));
            save(out_file, 'egg_clean_data', 'egg_time', 'egg_srate', 'egg_best_ch');

            h_egg = figure('Visible','off', 'Position', [0 0 800 600]);
            subplot(2,1,1); plot(t_egg, raw_data_egg(best_ch,:), 'Color',[0.7 0.7 0.7]); hold on;
            plot(t_egg, cleans(best_ch,:), 'k');
            title([subj ' ' cond_name ' Best EGG (Ch ' num2str(best_ch) ')']); axis tight;
            subplot(2,1,2); hold on; colors={'r','g','b','m'};
            for k=1:4, plot(all_f{k}, all_pxx{k}, 'Color', colors{k}, 'LineWidth', (k==best_ch)*1.5+0.5); end
            xlim([0.01 0.1]); xline(0.033,'k--'); xline(0.067,'k--'); title('PSD Selection');
            saveas(h_egg, fullfile(dirs.figs, sprintf('%s_%s_01_EGG.png', subj, cond_name))); close(h_egg);
            close all;

            % ----------------------------------------------------------------
            % B2. Build FieldTrip raw structures for EGG (1 chan) and EEG
            % ----------------------------------------------------------------
            % Wrap the cleaned single-channel EGG signal into a minimal
            % FieldTrip "raw" data structure (one trial spanning the
            % whole recording) so it can be passed to ft_preprocessing.
            ftEGG = [];
            ftEGG.label   = {'EGG'};
            ftEGG.trial{1} = egg_clean_data;
            ftEGG.time{1}  = egg_time;
            ftEGG.fsample  = egg_srate;
            % Load the cleaned, multi-channel EEG and resample it to the
            % common working sample rate (125 Hz) if it isn't already.
            EEG = pop_loadset('filename', file_eeg, 'filepath', path_cleaned);
            if round(EEG.srate) ~= work_srate_eeg
                EEG = pop_resample(EEG, work_srate_eeg);
            end
            % Convert the EEGLAB EEG structure into a FieldTrip raw
            % structure so ft_preprocessing can be used for filtering.
            ftEEG = eeglab2fieldtrip(EEG, 'raw', 'none');
            eeg_time    = ftEEG.time{1};
            chan_labels = ftEEG.label;
            n_chans     = numel(chan_labels);

            % Common analysis time grid (covers the overlap of both recordings)
            % EGG and EEG recordings may not start/end at exactly the same
            % time or sample rate, so build one shared time vector at
            % pac_calc_srate spanning only the overlapping duration.
            t0 = max(egg_time(1), eeg_time(1));
            t1 = min(egg_time(end), eeg_time(end));
            t_common = t0 : 1/pac_calc_srate : t1;
            % Trim edge_trim_sec seconds from each end of this common grid
            % to discard samples likely contaminated by filter edge
            % artifacts (ringing at the start/end of the filtered signal).
            trim_samps = round(edge_trim_sec * pac_calc_srate);
            valid_idx  = (trim_samps+1) : (length(t_common)-trim_samps);

            % ----------------------------------------------------------------
            % B3. Precompute EEG amplitude envelopes for every amp band
            %     (independent of phase band -> computed once, reused)
            % ----------------------------------------------------------------
            % For each amplitude band in the grid: bandpass-filter the EEG
            % (4th-order Butterworth) and take the Hilbert amplitude
            % envelope (cfg.hilbert = 'abs'), for every channel at once.
            amp_series = cell(n_amp,1);
            for a = 1:n_amp
                cfg = [];
                cfg.bpfilter   = 'yes';
                cfg.bpfreq     = amp_bands(a,:);
                cfg.bpfilttype = 'but';
                cfg.bpfiltord  = 4;
                cfg.hilbert    = 'abs';
                out = ft_preprocessing(cfg, ftEEG);
                tmp = out.trial{1};                      % channels x time (work_srate_eeg)
                
                % Interpolate each channel's envelope from the EEG's native
                % time base onto the shared t_common grid.
                tmp_common = zeros(n_chans, length(t_common));
                for ch = 1:n_chans
                    tmp_common(ch,:) = interp1(eeg_time, tmp(ch,:), t_common, 'pchip', 'extrap');
                end
                % Keep only the trimmed (artifact-free) portion, store for
                % reuse across all phase bands in Section B5.
                amp_series{a} = tmp_common(:, valid_idx);
                fprintf('   Amp band %d/%d (%.2f-%.2f Hz) filtered.\n', a, n_amp, amp_bands(a,1), amp_bands(a,2));
            end

            % ----------------------------------------------------------------
            % B4. Precompute EGG phase for every phase (driver) band
            % ----------------------------------------------------------------
            % For each phase-driver band in the grid: bandpass-filter the
            % EGG (3rd-order Butterworth, narrower band since these are
            % very low frequencies) and take the instantaneous phase via
            % Hilbert transform (cfg.hilbert = 'angle').
            phase_series = cell(n_phase,1);
            for p = 1:n_phase
                cfg = [];
                cfg.bpfilter   = 'yes';
                cfg.bpfreq     = phase_bands(p,:);
                cfg.bpfilttype = 'but';   % narrow low-freq band: check for ringing/instability;
                cfg.bpfiltord  = 3;       % if you see NaNs/instability, lower this further or switch to 'firws'
                cfg.hilbert    = 'angle';
                out = ft_preprocessing(cfg, ftEGG);
                tmp = out.trial{1};                       % 1 x time (egg_srate)
                % Interpolate phase safely by interpolating sin/cos
                % separately then recombining with atan2 - this avoids
                % the discontinuity/wrap-around problems of interpolating
                % a raw angle signal directly across the +-pi boundary.
                phase_c = atan2( ...
                    interp1(egg_time, sin(tmp), t_common, 'pchip', 'extrap'), ...
                    interp1(egg_time, cos(tmp), t_common, 'pchip', 'extrap'));
                % Keep only the trimmed portion, store for reuse across all
                % amplitude bands and channels in Section B5.
                phase_series{p} = phase_c(valid_idx);
                fprintf('   Phase band %d/%d (%.3f-%.3f Hz) filtered.\n', p, n_phase, phase_bands(p,1), phase_bands(p,2));
            end

            % ----------------------------------------------------------------
            % B5. Comodulogram: MI(channel, amp_freq, phase_freq)
            % ----------------------------------------------------------------
            % Preallocate the full comodulogram: one PAC value per
            % channel x amplitude band x phase band.
            MI_comod = zeros(n_chans, n_amp, n_phase);

            % Triple loop over phase bands, amplitude bands, and channels:
            % for every combination, compute the chosen PAC metric between
            % that phase-band's EGG phase series and that amplitude-band's
            % EEG envelope at that channel.
            for p = 1:n_phase
                ph = phase_series{p};
                for a = 1:n_amp
                    ampmat = amp_series{a};
                    for ch = 1:n_chans
                        switch pac_method
                            case 'tort'
                                MI_comod(ch,a,p) = calc_MI_tort(ph(:), ampmat(ch,:)', n_bins_MI);
                            case 'ozkurt'
                                MI_comod(ch,a,p) = calc_MI_ozkurt(ph(:), ampmat(ch,:)');
                            case 'canolty'
                                MI_comod(ch,a,p) = calc_MI_canolty(ph(:), ampmat(ch,:)');
                            case 'plv'
                                MI_comod(ch,a,p) = cohen_PLV(ph(:), ampmat(ch,:)');
                        end
                    end
                end
            end

            % ----------------------------------------------------------------
            % B6. Save comodulogram results for this subject/condition
            % ----------------------------------------------------------------
            EEG_chanlocs = EEG.chanlocs;
            out_pac = fullfile(dirs.stats, sprintf('%s_%s_Comodulogram.mat', subj, cond_name));
            save(out_pac, 'MI_comod', 'chan_labels', 'EEG_chanlocs', ...
                 'phase_bands', 'phase_centers', 'amp_bands', 'amp_centers', ...
                 'pac_method', 'n_bins_MI', 'pac_calc_srate', '-v7.3');

            % Quick diagnostic figure: comodulogram averaged across channels
            MI_avg = squeeze(mean(MI_comod, 1));   % amp x phase
            h = figure('Visible','on','Position',[0 0 700 550]);
            imagesc(phase_centers, amp_centers, MI_avg); axis xy;
            xlabel('EGG driver frequency (Hz)'); ylabel('EEG frequency (Hz)');
            title([subj ' ' cond_name ' - Comodulogram (channel average)']);
            colorbar; colormap(jet);
            saveas(h, fullfile(dirs.figs, sprintf('%s_%s_02_Comodulogram_avg.png', subj, cond_name)));
            close(h);

            fprintf('   [OK] %s %s: comodulogram computed (%d chans x %d amp x %d phase bins).\n', ...
                subj, cond_name, n_chans, n_amp, n_phase);

        catch ME
            fprintf('   [!] Errore %s: %s\n', cond_name, ME.message);
        end
    end
end

fprintf('\nAnalisi completata.\n');