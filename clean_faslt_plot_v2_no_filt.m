function [faslt, time_pow, norm_pow] =clean_faslt_plot_v2_no_filt(x, t,fs, freq_faslt, anomalies, diagno, enlarge, last_plot)

%[faslt, time_pow, norm_pow]=faslt_plot_complete(x,t, fs, freq_faslt, 0.01, 0.2, 'predictive', 'faslt', 'on', 'abs_val', diagno);
 [faslt, time_pow, norm_pow, buf_tot]=general_faslt_plot_no_filt(x,t, fs, freq_faslt, 0.01, 0.2, [], 'off', diagno, 'zero');

if ~isempty(anomalies)

    for i=1:2:length(anomalies)-1
        faslt(:, max(anomalies(i)-enlarge, 1):min(anomalies(i+1)+enlarge, length(x)))=0;
    end

    norm_pow= (sum((abs(faslt)/2).^2, 2));
    %time_pow = (sum((abs(faslt)/2).^2, 1));
    if size(freq_faslt, 1) == 1
        freq_faslt=freq_faslt';
    end
    time_pow = sum( ((abs(faslt)/2).^2)./freq_faslt, 1);

    for i=1:2:length(anomalies)-1
        time_pow(max(anomalies(i)-enlarge, 1):min(anomalies(i+1)+enlarge, length(x)))=NaN;
    end

%mean_norm= mean(norm_pow); 
peak_pow=max(norm_pow);
peak=freq_faslt(find(norm_pow==max(norm_pow)));
var_norm= var(norm_pow);

mean_time= mean(time_pow, 'omitnan');
var_time= var(time_pow, 'omitnan');

text_str = sprintf('Peak freq: %e \nPeak pow: %e \nVar pow sp: %e \nMean En: %e \nVar En: %e', peak, peak_pow, var_norm, mean_time, var_time);


%%%%%%%%%%%%----------------------------------------------

%{
figure;
tplot(clean_faslt_fragmented, t, freq_faslt);
hold on
xline(t(zeri), '--r');
if last_zero~=0
    xline(t(last_zero), '--r');
end
title('Clen FASLT on fragmented signal')

figure;
plot(t, x);

%}

if strcmp(diagno, 'on')

figure;
subplot(2, 1, 1);

%signal
plot(t, x);
hold on
xline(t(anomalies), '--k')
%{
for i = 1:2:length(times)-1
  if i+1 <= length(times)
      disp(i);
      disp(times(i));
      disp(times(i+1));
      disp(min(x)-500);
      disp(max(x)+500);
    %p=fill([times(i), times(i+1), times(i+1),times(i)], [min(x)-500, min(x)-500, max(x)+500, max(x)+500],'black','FaceAlpha',0.1);
    %p.EdgeColor = 'none';
    hold on
  end
end
%per riempire l'ultimo intervallo in caso di artefatto che non termina
%prima della fine del segnale
if mod(length(zeri), 2) == 1
    p=fill([times(end), t(end), t(end), times(end)], [min(x)-500, min(x)-500, max(x)+500, max(x)+500],'black','FaceAlpha',0.1);
    p.EdgeColor = 'none';
    hold on
end
%}
ax = gca;
ax.YLim            = [min(x)-500, max(x)+500];
ax = gca;
ax.XLim            = [t(1) t(end)];
ax.FontSize        = 16;
ax.YLabel.String   = 'µV';
ax.YLabel.FontSize = 17;
title('Raw signal Vs clean FASLT ')


subplot_2_Pos = get(gca, 'Position');
subplot_2_Pos(3) = subplot_2_Pos(3) + 0.14; 
subplot_2_Pos(1) = subplot_2_Pos(1) -0.07;
set(gca, 'Position', subplot_2_Pos);

%-----------------------------
%plot FASLT

subplot(2, 1, 2);

%adjusting position and height of box
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1), pos(2)+0.01, pos(3), pos(4)+0.02]); % Modifica l'altezza relativa

tplot(faslt, t, freq_faslt);
hold on
xline(t(anomalies), '--k');
%grid on;



ax = gca;
ax.XLim            = [t(1) t(end)];
ax.FontSize        = 16;
ax.XLabel.String =  'Time (s)';
ax.YLabel.String   = 'Frequency (Hz)';
ax.YLabel.FontSize = 17;
%yticks([-pi, 0, pi]);
%yticklabels({'-\pi', '0', '\pi'});
subplot_2_Pos = get(gca, 'Position');
subplot_2_Pos(3) = subplot_2_Pos(3) + 0.14; 
subplot_2_Pos(1) = subplot_2_Pos(1) -0.07;
set(gca, 'Position', subplot_2_Pos);


end
%size(time_pow)
%size(t)

%%%%%----------------------------------------------------

if strcmp(last_plot, 'on')

figure;
%------------------------------------------------
subplot(6, 6, [6, 12, 18, 24, 30])

plot(norm_pow, freq_faslt);
grid on;
ylim([min(freq_faslt) max(freq_faslt)]);
f_axes = gca; % Get the current axes object of subplot 2
y_ticks_freq = get(f_axes, 'YTick');
set(f_axes, 'FontSize', 12);
title('Power spectrum');
%------------------------------------------------
subplot(6, 6, [31:35])

plot(t, time_pow);
grid on;
xline(t(anomalies), '--k');
xlim([t(1) t(end)])
t_axes = gca; % Get the current axes object of subplot 2
x_ticks_time = get(t_axes, 'XTick');
xlabel('Time (s)');
set(t_axes, 'FontSize', 12);
title('Energy per unit time');

subplot(6,6,[1:5, 7:11, 13:17, 19:23, 25:29]);

tplot(faslt, t, freq_faslt);
xline(t(anomalies), '--k');
xticks(x_ticks_time);
xticklabels([]);
yticks(y_ticks_freq);
yticklabels([]);
grid on;
ylabel('Frequency (Hz)');
title('FASLT TFR');

subplot(6, 6, 36, 'align')
annotation('textbox','String',text_str,'FontSize', 9, 'Position', [0.78 0.11 0.1 0.1]); % Adjust font size as desired
axis off;
%sgtitle('FASLT after artifacts cancellation')

%%%--------------------------------------------------------------

end

end



end