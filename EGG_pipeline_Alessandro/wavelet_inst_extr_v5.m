function [efreq, eamp, ephi, tfrtype, amp_curve, ph_curve, if_curve, anom_new] = ...
    wavelet_inst_extr_v5(x1, t, anom_intervals, plots, fs, opts)
%WAVELET_INST_EXTR_V5  Instantaneous frequency/amplitude/phase of the gastric
% component via ridge extraction on the WT/WFT (Iatsenko framework).
%
% Drop-in successor of WAVELET_INST_EXTR_V4. It keeps the SAME segmentation
% logic (curves are estimated only on the non-anomalous segments, so the
% filled/anomalous intervals stay NaN and are effectively removed from the
% output), and adds a choice of the FREQUENCY BAND on which the component is
% tracked, WITHOUT ever changing the time-frequency RESOLUTION (f0). Keeping
% f0 fixed across subjects/segments is what makes the extracted features
% comparable; the band only decides where we look.
%
%   opts.method  : 'normo'    - fixed normogastric band [0.033 0.067] Hz
%                               (reproduces v4; default if opts omitted).
%                  'adaptive' - band = [f_dom*bf(1), f_dom*bf(2)] centred on
%                               the subject/segment dominant frequency,
%                               clamped to the gastric range (METHOD a).
%                  'wideband' - full gastric band + heavier ridge penalty,
%                               letting the path optimisation reject noise
%                               instead of narrowing the band (METHOD b).
%   opts.f_dom       : dominant frequency in Hz (REQUIRED for 'adaptive').
%                      Extract it from your pipeline's parameters struct, e.g.
%                          fd = get_scalar(parameters.Dominant_frequency);
%   opts.band_factor : [lo hi] multiplicative window for 'adaptive'
%                      (default [0.5 2], i.e. half to double f_dom).
%   opts.gastric_band: [fmin fmax] full gastric range (default [0.016 0.16]).
%   opts.penal       : [alpha beta] ecurve Method-2 penalty weights on the
%                      frequency deviation and its time-derivative
%                      (ecurve default [1 1]; 'wideband' default [2 2]).
%   opts.f0_wft      : WFT resolution parameter (default 50, as v4).
%   opts.f0_wt       : WT  resolution parameter (default f0_wft/20 = 2.5).
%
% The WT-vs-WFT selection (checktype) is preserved. See also WT, WFT,
% ECURVE, RECTFR, BESTEST (Iatsenko et al.).

% -------------------- options with backward-compatible defaults ----------
if nargin < 6 || isempty(opts), opts = struct(); end
opts = set_default(opts, 'method',       'normo');
opts = set_default(opts, 'gastric_band', [0.016 0.16]);
opts = set_default(opts, 'band_factor',  [1/1.5 1.5]);
opts = set_default(opts, 'f0_wft',       50);
opts = set_default(opts, 'f0_wt',        opts.f0_wft/20);
switch lower(opts.method)
    case 'wideband', opts = set_default(opts, 'penal', [2 2]);
    otherwise,       opts = set_default(opts, 'penal', [1 1]);
end

% -------------------- decide the extraction band -------------------------
switch lower(opts.method)
    case 'normo'
        band = [0.033 0.067];
    case 'adaptive'
        if ~isfield(opts,'f_dom') || isempty(opts.f_dom) || ~isfinite(opts.f_dom)
            error('wavelet_inst_extr_v5:fdom', ...
                  'opts.f_dom (dominant frequency, Hz) is required for method ''adaptive''.');
        end
        band = [opts.f_dom*opts.band_factor(1), opts.f_dom*opts.band_factor(2)];
        band = [max(band(1), opts.gastric_band(1)), min(band(2), opts.gastric_band(2))];
    case 'wideband'
        band = opts.gastric_band;
    otherwise
        error('wavelet_inst_extr_v5:method', 'Unknown opts.method ''%s''.', opts.method);
end
fmin = band(1); fmax = band(2);

% -------------------- initialise outputs ---------------------------------
if_curve  = nan(1, numel(t));
amp_curve = nan(1, numel(t));
ph_curve  = nan(1, numel(t));
anom_new  = anom_intervals;

% -------------------- compute both TFRs ONCE (f0 fixed) ------------------
% Resolution (f0) is held fixed regardless of the band -> comparable features.
[WT,  freq,     wopt]     = wt (x1, fs, 'fmin', fmin, 'fmax', fmax, 'f0', opts.f0_wt,  'Display', 'off');
[WFT, freq_wft, wopt_wft] = wft(x1, fs, 'fmin', fmin, 'fmax', fmax, 'f0', opts.f0_wft, 'Display', 'off');

pen = opts.penal;   % [alpha beta] penalty weights passed to ecurve

% -------------------- segmented extraction (identical to v4) -------------
if ~isempty(anom_intervals)
    % before the first anomaly
    if anom_intervals(1) > 1
        idx = 1:anom_intervals(1);
        [eamp, ephi, efreq, tfrtype] = process_segment(WT, freq, wopt, WFT, freq_wft, wopt_wft, fs, idx, pen);
        [if_curve, amp_curve, ph_curve] = assign(if_curve, amp_curve, ph_curve, idx, efreq, eamp, ephi);
    end
    % between consecutive anomalies
    for i = 2:2:length(anom_intervals) - 1
        idx = anom_intervals(i):anom_intervals(i+1);
        [eamp, ephi, efreq, tfrtype] = process_segment(WT, freq, wopt, WFT, freq_wft, wopt_wft, fs, idx, pen);
        [if_curve, amp_curve, ph_curve] = assign(if_curve, amp_curve, ph_curve, idx, efreq, eamp, ephi);
    end
    % after the last anomaly
    if anom_intervals(end) < length(t)
        idx = anom_intervals(end):length(t);
        [eamp, ephi, efreq, tfrtype] = process_segment(WT, freq, wopt, WFT, freq_wft, wopt_wft, fs, idx, pen);
        [if_curve, amp_curve, ph_curve] = assign(if_curve, amp_curve, ph_curve, idx, efreq, eamp, ephi);
    end
else
    % whole signal, no anomalies
    tfsupp = ecurve(WT, freq, fs, 'Display', 'off', 'Method', 2, 'Param', pen, 'PathOpt', 'on');
    [eamp, ephi, efreq] = bestest(tfsupp, WT, freq, wopt, 'Display', 'off');
    tfrtype = checktype(fs, eamp, efreq, 'off');
    if strcmpi(tfrtype, 'WFT')
        tfsupp = ecurve(WFT, freq_wft, fs, 'Display', 'off', 'Method', 2, 'Param', pen, 'PathOpt', 'on');
        [eamp, ephi, efreq] = bestest(tfsupp, WFT, freq_wft, wopt_wft, 'Display', 'off');
    end
    if_curve = efreq;  amp_curve = eamp;  ph_curve = ephi;
end

% -------------------- plotting (unchanged look) --------------------------
if strcmp(plots, 'on')
    plot_curves(t, amp_curve, ph_curve, if_curve, anom_intervals, opts, band);
end
end

% =========================================================================
function [eamp, ephi, efreq, tfrtype] = process_segment(WT, freq, wopt, WFT, freq_wft, wopt_wft, fs, idx, pen)
% Ridge extraction on one segment. Uses ecurve Method 2 with path
% optimisation and the [alpha beta] penalty weights in PEN: larger weights
% penalise frequency deviations/jumps more strongly, so the ridge stays on
% the smooth component instead of hopping onto transient noise peaks.
tfsupp = ecurve(WT(:, idx), freq, [fs,1], 'Display', 'off', 'Method', 2, 'Param', pen, 'PathOpt', 'on');
[eamp, ephi, efreq] = rectfr(tfsupp, WT(:, idx), freq, wopt, 'ridge');
tfrtype = checktype(fs, eamp, efreq, 'off');
if strcmpi(tfrtype, 'WFT')
    tfsupp = ecurve(WFT(:, idx), freq_wft, [fs,1], 'Display', 'off', 'Method', 2, 'Param', pen, 'PathOpt', 'on');
    [eamp, ephi, efreq] = rectfr(tfsupp, WFT(:, idx), freq_wft, wopt_wft, 'ridge');
end
end

% =========================================================================
function [ifc, ampc, phc] = assign(ifc, ampc, phc, idx, efreq, eamp, ephi)
ifc(idx) = efreq;  ampc(idx) = eamp;  phc(idx) = ephi;
end

% =========================================================================
function s = set_default(s, field, val)
if ~isfield(s, field) || isempty(s.(field)), s.(field) = val; end
end

% =========================================================================
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

% =========================================================================
function plot_curves(t, amp_curve, ph_curve, if_curve, anom_intervals, opts, band)
figure; set(gcf, 'Position', [100, 100, 1000, 700]);
anomaly_color = [1 0 0]; anomaly_alpha = 0.15;
color_signal = [0.0 0.45 0.74]; color_amp = [0.85 0.33 0.1];
color_phase  = [0.2 0.7 0.4];   color_freq = [0.85 0.65 0.13];

subplot(3,1,1)
plot(t, amp_curve .* cos(ph_curve), 'Color', color_signal, 'DisplayName', 'Reconstructed Signal', 'LineWidth', 1.5)
xlim([t(1), t(end)]); hold on
plot(t, amp_curve, 'Color', color_amp, 'DisplayName', 'Amplitude Envelope', 'LineWidth', 1.5)
drawnow; add_anomaly_patches(t, anom_intervals, anomaly_color, anomaly_alpha);
grid on; ylabel('Amplitude (\muV)', 'FontSize', 12)
title(sprintf('Reconstructed Signal and Amplitude Envelope  [%s band: %.3f-%.3f Hz]', ...
      opts.method, band(1), band(2)), 'FontSize', 13)
legend('Location', 'best'); hold off

subplot(3,1,2)
plot(t, wrapToPi(ph_curve), 'Color', color_phase, 'DisplayName', 'Phase', 'LineWidth', 1.5)
xlim([t(1), t(end)]); hold on; drawnow
add_anomaly_patches(t, anom_intervals, anomaly_color, anomaly_alpha);
grid on; yticks([-pi, -pi/2, 0, pi/2, pi]); yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
ylabel('Phase (rad)', 'FontSize', 12); title('Instantaneous Phase', 'FontSize', 14); hold off

subplot(3,1,3)
plot(t, if_curve, 'Color', color_freq, 'DisplayName', 'Instantaneous Frequency', 'LineWidth', 1.5)
xlim([t(1), t(end)]); hold on; drawnow
add_anomaly_patches(t, anom_intervals, anomaly_color, anomaly_alpha);
grid on; ylabel('Frequency (Hz)', 'FontSize', 12); xlabel('Time (s)', 'FontSize', 12)
title('Instantaneous Frequency', 'FontSize', 14); hold off
end

% =========================================================================
function add_anomaly_patches(t, anom_intervals, c, a)
if isempty(anom_intervals), return; end
yl = ylim;
for i = 1:2:length(anom_intervals)
    xs = t(anom_intervals(i)); xe = t(anom_intervals(i+1));
    patch([xs xe xe xs], [yl(1) yl(1) yl(2) yl(2)], c, ...
          'FaceAlpha', a, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    xline(xs, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xline(xe, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
end
end

% -------------------- Helper function: determine optimal TFR type --------------------
function tfrtype = checktype(fs, iamp, ifreq, DispMode)

    Perc = 0.75;
    DLev = 1.1;
    L = length(iamp);
    i1 = round((0.5 - Perc/2) * L);
    i2 = round((0.5 + Perc/2) * L);

    % Estimate time-derivatives of amplitude and frequency
    tamp = (iamp(1:end-2) + iamp(2:end-1) + iamp(3:end)) / 3;
    dtamp = fs * (iamp(3:end) - iamp(1:end-2)) / 2;
    tfreq = (ifreq(1:end-2) + ifreq(2:end-1) + ifreq(3:end)) / 3;
    dtfreq = fs * (ifreq(3:end) - ifreq(1:end-2)) / 2;

    gamp1 = sort(abs(hilbert(dtamp./tfreq - mean(dtamp./tfreq))));
    gamp2 = sort(abs(hilbert((dtamp - mean(dtamp)) / mean(tfreq))));
    gfreq1 = sort(abs(hilbert(dtfreq./tfreq - mean(dtfreq./tfreq))));
    gfreq2 = sort(abs(hilbert((dtfreq - mean(dtfreq)) / mean(tfreq))));

    Vamp = (gamp1(i2) - gamp1(i1)) / (gamp2(i2) - gamp2(i1));
    Vfreq = (gfreq1(i2) - gfreq1(i1)) / (gfreq2(i2) - gfreq2(i1));

    if isnan(Vamp), Vamp = 1; end
    if isnan(Vfreq), Vfreq = 1; end

    U = 1 / (1 + Vamp) + 1 / (1 + Vfreq);

    % Decide TFR type
    if U < DLev
        tfrtype = 'WFT'; cstr = '<';
    else
        tfrtype = 'WT';  cstr = '>';
    end

    if ~strcmpi(DispMode, 'off')
        fprintf(['Optimal TFR type was determined to be ', tfrtype, ' (']);
        fprintf('Va=%0.3f, Vf=%0.3f, 1/(1+Va)+1/(1+Vf)=%0.3f', Vamp, Vfreq, U);
        fprintf([cstr, num2str(DLev), ')\n']);
    end
end