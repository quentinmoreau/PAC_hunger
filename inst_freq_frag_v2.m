function [if_curve, perc]=inst_freq_frag_v2(faslt_frag, freq_faslt, anom_intervals, t, smooth, plots, ist, numPart)

if_curve=nan(1, size(faslt_frag, 2));

if ~isempty(anom_intervals)

if anom_intervals(1)>1
    %[tfsupp] = ecurve(faslt_frag(:, 1:anom_intervals(1)), freq_faslt, 10);
    %[tfsupp] = ecurve(faslt_frag(:, 1:anom_intervals(1)), freq_faslt, [10, 1], 'Method', 1);
    %if_curve(1:anom_intervals(1))=tfsupp(1, :);
    %[Cs, Es] = curve_ext_multi(faslt_frag(:, 1:anom_intervals(1)), log2(freq_faslt), 1, smooth, 1);
    [Cs, Es] = curve_ext_MATLAB(faslt_frag(:, 1:anom_intervals(1)), log2(freq_faslt), smooth);
    %[Cs, Es] = exridge_SAMD(faslt_frag(:, 1:anom_intervals(1)), smooth, smooth, 5);
    %[c e] = exridge(Tx, lambda)
    if_curve(1:anom_intervals(1))=freq_faslt(Cs);
end


for i=2:2:length(anom_intervals)-1

%[tfsupp] = ecurve(faslt_frag(:, anom_intervals(i):anom_intervals(i+1)), freq_faslt, 10);
%[tfsupp] = ecurve(faslt_frag(:, anom_intervals(i):anom_intervals(i+1)), freq_faslt, [10, 1], 'Method', 1);
%if_curve(anom_intervals(i):anom_intervals(i+1))=tfsupp(1, :);
[Cs, Es] = curve_ext_MATLAB(faslt_frag(:, anom_intervals(i):anom_intervals(i+1)), log2(freq_faslt), smooth);
%[Cs, Es] = curve_ext_multi(faslt_frag(:, anom_intervals(i):anom_intervals(i+1)), log2(freq_faslt), 1, smooth, 1); 
%[Cs, Es] = exridge_SAMD(faslt_frag(:, anom_intervals(i):anom_intervals(i+1)), smooth, smooth, 5); 
    if_curve(anom_intervals(i):anom_intervals(i+1))=freq_faslt(Cs);
end

if anom_intervals(end)<size(faslt_frag, 2)
    %[tfsupp] = ecurve(faslt_frag(:, anom_intervals(end):end), freq_faslt, 10);
    %[tfsupp] = ecurve(faslt_frag(:, anom_intervals(end):end), freq_faslt, [10, 1], 'Method', 1);
    %if_curve(anom_intervals(end):end)=tfsupp(1, :);
    [Cs, Es] = curve_ext_MATLAB(faslt_frag(:, anom_intervals(end):end), log2(freq_faslt), smooth);
    %[Cs, Es] = curve_ext_multi(faslt_frag(:, anom_intervals(end):end), log2(freq_faslt), 1, smooth, 1); 
    %[Cs, Es] = exridge_SAMD(faslt_frag(:, anom_intervals(end):end), smooth, smooth, 5); 
    if_curve(anom_intervals(end):end)=freq_faslt(Cs);
end

else

%[Cs, Es] = curve_ext_multi(faslt_frag, log2(freq_faslt), 1, smooth, 1); 
%[Cs, Es] = exridge_SAMD(faslt_frag, smooth, smooth, 5); 
    %if_curve=freq_faslt(Cs);
    [Cs, Es] = curve_ext_MATLAB(faslt_frag, log2(freq_faslt), smooth);
    %size(Cs)
    %Cs(end)
    %Cs(1)
    if_curve=freq_faslt(Cs);
    %[tfsupp] = ecurve(faslt_frag, freq_faslt, [10, 1], 'Method', 1);
    %if_curve=tfsupp(1, :);

end

%{
figure;
plot(t, if_curve)
xline(t(anom_intervals), '--k')
title('Isntantaneous frequencies from clean TFR')
%}


norm=find(if_curve<0.067 & if_curve>0.033);
tac=find(if_curve>0.067);
brad=find(if_curve<0.033);
perc=length(norm)/(length(if_curve)-sum(isnan(if_curve)));


if strcmp(plots, 'on')

    faslt_plot=faslt_frag;

    for i=1:2:length(anom_intervals)-1
        falst_plot(:, anom_intervals(i):anom_intervals(i+1))=NaN;
    end

    t(1)
    t(end)
    %freq_faslt
    size(t)
    size(freq_faslt)
    size(faslt_plot)

    figure;
    imagesc(t, freq_faslt, abs(faslt_plot))

    figure;
    imagesc(abs(faslt_plot))

    figure;
tplot(faslt_plot, t, freq_faslt);

figure
plot(t, if_curve, 'o', 'MarkerSize',1, 'Color', [1, 0.8, 0]);

figure;
tplot(faslt_plot, t, freq_faslt);
hold on
plot(t, if_curve, 'o', 'MarkerSize',1, 'Color', [1, 0.8, 0]);
hold off

figure;
tplot(faslt_plot, t, freq_faslt);
hold on
%plot(t, if_curve, 'o', 'MarkerSize',1, 'Color', 'green');
plot(t, if_curve, 'o', 'MarkerSize',1, 'Color', [1, 0.8, 0]);
%[0.2, 1, 0.2]
xlabel('Time (seconds)')
ylabel('Frequency (Hz)')
hold on
plot(t(tac), if_curve(tac), 'o', 'MarkerSize',1, 'Color', 'red');
plot(t(brad), if_curve(brad), 'o', 'MarkerSize',1, 'Color', 'red');
if ~isempty(anom_intervals)
xline(t(anom_intervals), '--k')
end
%title(['ridge points (normo blue, red tachi or bradi) - % of normogastric ridges: ', num2str(perc)]);
title(['FASLT for ch', numPart, ' with time-frequency ridges - (yellow - normogastric | red - tachi or bradigastric )'])
subtitle(['Fraction of normogastric ridges: ', sprintf('%.3f', perc), ' - mean ridge freq: ', sprintf('%.4f',mean(if_curve, "omitnan")), ' - std ridge freq: ', sprintf('%.4f',std(if_curve, "omitnan")) ])
hold off

end

if strcmp(ist, 'on')

%size(if_curve)
A_e_iphi=nan(1, length(t));
for i=1:length(t)
    if ~isnan(if_curve(i))
        ind=find( abs(if_curve(i)-freq_faslt) == min(abs(if_curve(i)-freq_faslt)) );
        A_e_iphi(i)=faslt_frag(ind, i);
    end
end

iamp=abs(A_e_iphi);
iphi=angle(A_e_iphi);

%figure;
%plot(t, iphi)

figure;
plot(t, (2*pi*iphi))
if ~isempty(anom_intervals)
xline(t(anom_intervals), '--k')
end
title('Instantaneous phase of the signal')

figure;
plot(t, iamp.*cos(unwrap(2*pi*iphi)))
hold on
plot(t, iamp)
if ~isempty(anom_intervals)
xline(t(anom_intervals), '--k')
end
title('Reconstructed gastric signal with its instantaneous amplitude')

end

end
