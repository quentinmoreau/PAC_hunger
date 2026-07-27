function [faslt, time_pow, norm_pow, buffer]=general_faslt_plot_mute(x,t, fs, freq_faslt, fmin, fmax, anomalies, ssa, plots, padding)

if strcmp(ssa, 'on')
    [LT, ST, R]=trenddecomp(x, 'ssa', 60*fs);
    
    %{
if strcmp(plots, 'on')
    figure;
    plot(t, x);
    hold on
    plot(t, LT, 'LineWidth', 1.5, 'Color', 'r');
    title('Evaluating the quality of ssa detrending');
end
    %}

    x=x-LT';
end

x=cheb_EGG_filt(x, fs);



%x= keep_stomach_range(x, fs, fmin, fmax, 4);


% qui sto filtrando il segnale senza padding, spero che sia giusto

%{
L=length(x);
%Filtering
    fx=fft(x,L); % Fourier transform of a signal
    Nq=ceil((L+1)/2); ff=[(0:Nq-1),-fliplr(1:L-Nq)]*fs/L; ff=ff(:); % frequencies in Fourier transform
    fx(abs(ff)<=max([fmin,fs/L]) | abs(ff)>=fmax)=0; % filter signal in a chosen frequency domain
    x=ifft(fx);
%}

    %x = keep_stomach_range(x, fs, 0.01, 0.2, 2);


if isempty(anomalies)

    % ---------------- The signal is clean

    [faslt, buffer] = faslt_complex_calc_no_filt_mute(x, fs, freq_faslt, padding);

else

    % ---------------- There are artifact to remove

    faslt_all=[];

    if anomalies(1)>1

        %disp('siamo qui')
        %anomalies(1)
    x_new= x(1:anomalies(1)-1);
    
    [faslt_all, buffer] = faslt_complex_calc_no_filt(x_new, fs, freq_faslt,padding);
    
        disp(['The number of artifact processed is 1 over ', num2str(length(anomalies)/2)]);
    else
        disp(['The number of artifact processed is 1 over ', num2str(length(anomalies)/2)]);
    end

    if length(anomalies)>2
        for i=2:2:length(anomalies)-1
            %i
            %anomalies(i)
            %anomalies(i+1)
    
            x_new=x(anomalies(i):anomalies(i+1));

            [faslt_step, buffer] = faslt_complex_calc_no_filt(x_new, fs, freq_faslt, padding);
    
            fill=zeros(length(freq_faslt), length(anomalies(i-1):anomalies(i))-2);
            faslt_all=[faslt_all, fill, faslt_step];
    
            disp(['The number of artifact processed is ', num2str(i/2+1), ' over ', num2str(length(anomalies)/2)]);
    
        end

        %riempio l'ulitmo intervallo di zeri
        fill=zeros(length(freq_faslt), length(anomalies(end-1):anomalies(end))-1);
        faslt_all=[faslt_all, fill];



        if anomalies(end)~=length(x)
            x_new= x(anomalies(end):end);
            [faslt_step, buffer] = faslt_complex_calc_no_filt(x_new, fs, freq_faslt, padding);
    
            faslt_all= [faslt_all, faslt_step];
    
        else
            size(faslt_all)
            faslt_all=[faslt_all, zeros(size(faslt_all, 1), 1)];
            size(faslt_all)
    
        end
        %anomalies

    else
        if anomalies(2)~=length(x)
            x_new= x(anomalies(2):end);
            [faslt_step, buffer] = faslt_complex_calc_no_filt(x_new, fs, freq_faslt, padding);
            fill=zeros(length(freq_faslt), length(anomalies(1):anomalies(2))-1);
            faslt_all=[faslt_all fill faslt_step];
            
        else
            fill= zeros(length(freq_faslt), length(anomalies(1):anomalies(2)));
            faslt_all=[faslt_all fill];
       
        end
    
    end

faslt=faslt_all;



end

% --------------------- end of falst calulation

%size(freq_faslt)
if size(freq_faslt, 1) == 1
    freq_faslt=freq_faslt';
end
%size(freq_faslt)
%size(faslt)

if isempty(anomalies)

    norm_pow= (sum((abs(faslt)/2).^2, 2));
    %time_pow = (sum((abs(faslt)/2).^2, 1));
    %TFR= abs(faslt).^2;
    time_pow = sum( ((abs(faslt)/2).^2)./freq_faslt, 1);

else

     for i=1:2:length(anomalies)-1
        faslt(:, anomalies(i):anomalies(i+1))=0;
    end

    norm_pow= (sum((abs(faslt)/2).^2, 2));
    %time_pow = (sum((abs(faslt)/2).^2, 1));
    time_pow = sum( ( (abs(faslt)/2).^2)./freq_faslt, 1);

    for i=1:2:length(anomalies)-1
        %time_pow(max(anomalies(i)-enlarge, 1):min(anomalies(i+1)+enlarge, length(x)))=NaN;
        time_pow(:, anomalies(i):anomalies(i+1))=NaN;
    end
    

end



%mean_norm= mean(norm_pow); 
peak_pow=max(norm_pow);
peak=freq_faslt(find(norm_pow==max(norm_pow)));
var_norm= var(norm_pow);

mean_time= mean(time_pow, 'omitnan');
var_time= var(time_pow, 'omitnan');

text_str = sprintf('Peak freq: %e \nPeak pow: %e \nVar pow sp: %e \nMean En: %e \nVar En: %e', peak, peak_pow, var_norm, mean_time, var_time);


%size(time_pow)
%size(t)


if strcmp(plots, 'on')

figure;
%------------------------------------------------
subplot(6, 6, [6, 12, 18, 24, 30])

plot(norm_pow, freq_faslt);
grid on;
ylim([min(freq_faslt) max(freq_faslt)]);
f_axes = gca; % Get the current axes object of subplot 2
y_ticks_freq = get(f_axes, 'YTick');
set(f_axes, 'FontSize', 12);
title('Power spectrum - PS');
%------------------------------------------------
subplot(6, 6, [31:35])

plot(t, time_pow);
if~isempty(anomalies)
    xline(t(anomalies), '--k')
end
grid on;
xlim([t(1) t(end)])
t_axes = gca; % Get the current axes object of subplot 2
x_ticks_time = get(t_axes, 'XTick');
xlabel('Time (s)');
set(t_axes, 'FontSize', 12);
title('Energy per unit time - E(t)');

subplot(6,6,[1:5, 7:11, 13:17, 19:23, 25:29]);

tplot(faslt, t, freq_faslt);

if~isempty(anomalies)
    xline(t(anomalies), '--k')
end

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
end

end