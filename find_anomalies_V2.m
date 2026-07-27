function [filt_buff, filtered_signal, an_edges_filt_opt]= find_anomalies_V2(x, t, wind, fs)

%wind
%filtered_signal = keep_stomach_range(x, fs, 0.016, 0.16, 2);
filtered_signal=cheb_EGG_filt(x, fs);

if isrow(filtered_signal)
    % Se è un vettore riga, trasponi per ottenere un vettore colonna
    filtered_signal = filtered_signal';
end
%L=length(wind);
%L/2

filt_buff=buffer(filtered_signal, wind,wind-1 , 'nodelay');
%filt_buff_cento=buffer(filtered_signal, 100,99 , 'nodelay');
%size(filt_buff)

%r_tot=isanomaly_FE(filtered_signal, 'liberal tukey', Inf);
r_window=isanomaly_FE(filtered_signal, 'liberal tukey', wind);

if size(r_window, 1)>size(r_window, 2)
    r_window=r_window';
end

r_tot=isanomaly_FE(filtered_signal, 'liberal tukey', Inf);

if size(r_tot, 1)>size(r_tot, 2)
    r_tot=r_tot';
end
%size(r_tot)


% Crea un vettore logico per i limiti dell'intervallo
boundary_mask = (1:length(filtered_signal) <= floor(wind/2)) | (1:length(filtered_signal) >= length(filtered_signal) - floor(wind/2));
%size(boundary_mask)

% Applica la condizione in modo vettoriale
r1 = (r_tot == 1 & boundary_mask) | (r_window == 1);

if size(r1, 1)>size(r1, 2)
    r1=r1';
end
%size(r1)

%figure;
%plot(r1, 'Color', 'b')
%hold on

%r1(1:50)

r1 = merge_short_gaps(r1, 20);
%plot(r1, 'Color', 'r')
%hold off

%r1(1:50)

if size(r1, 1)>size(r1, 2)
    r1=r1';
end

an_edges_filt_tot = find(diff([0, r1, 0]) ~= 0);

    % Rimuovi eventuali indici che vanno oltre la lunghezza della serie
    % siccede se l'aomalia finisce alla fine della serie
    an_edges_filt_tot(an_edges_filt_tot > length(x)) = length(x);


    [Mdl,tf,scores] = rrcforest(filtered_signal,ContaminationFraction=sum(r1)/length(r1), ...
        CollusiveDisplacement="average", NumLearners=100, NumObservationsPerLearner=256);
 
   
    r_window_sc=isanomaly_FE(scores, 'conservative tukey', Inf);


if size(r_window_sc, 1)>size(r_window_sc, 2)
    r_window_sc=r_window_sc';
end

r_window_sc = merge_short_gaps(r_window_sc, 20);

if size(r_window_sc, 1)>size(r_window_sc, 2)
    r_window_sc=r_window_sc';
end

an_edges_filt_sc = find(diff([0, r_window_sc, 0]) ~= 0);

    % Rimuovi eventuali indici che vanno oltre la lunghezza della serie
    % siccede se l'aomalia finisce alla fine della serie
    an_edges_filt_sc(an_edges_filt_sc > length(x)) = length(x);


    %r1_updated = anomaly_intersection(r1, r_window_sc);
    r1_updated = anomaly_intersection(r_window_sc, r1);

    an_edges_filt_opt = find(diff([0, r1_updated, 0]) ~= 0);

    % Rimuovi eventuali indici che vanno oltre la lunghezza della serie
    % siccede se l'aomalia finisce alla fine della serie
    an_edges_filt_opt(an_edges_filt_opt > length(x)) = length(x);





%{


   figure;
   subplot(3, 1, 1)
plot(t, filtered_signal)
if ~isempty(an_edges_filt_sc)
    xline(t(an_edges_filt_sc), '--k');
end
xlim([0, t(end)])
title('filtered signal with random cut forest anomalies')
subplot(3, 1, 2)
plot(t, filtered_signal)
if ~isempty(an_edges_filt_tot)
    xline(t(an_edges_filt_tot), '--k');
end
xlim([0, t(end)])
title('filtered signal with tukey anomalies')


subplot(3, 1, 3)
plot(t, filtered_signal)
if ~isempty(an_edges_filt_opt)
    xline(t(an_edges_filt_opt), '--k');
end
xlim([0, t(end)])
title('filtered signal with selected anomalies')

%}


%{

figure;

% Margine per estendere i limiti verticali
margin_factor = 0.05; % 5% del range verticale

% Subplot 1: Random Cut Forest anomalies
subplot(3, 1, 1)
plot(t, filtered_signal, 'b')
hold on
y_limits = ylim; % Ottieni i limiti dell'asse Y
range_y = y_limits(2) - y_limits(1);
extended_y_limits = [y_limits(1) - margin_factor * range_y, y_limits(2) + margin_factor * range_y];
if ~isempty(an_edges_filt_sc)
    for i = 1:2:length(an_edges_filt_sc)-1
        % Aggiunge una banda rossa per l'anomalia
        fill([t(an_edges_filt_sc(i)), t(an_edges_filt_sc(i)), t(an_edges_filt_sc(i+1)), t(an_edges_filt_sc(i+1))], ...
             [extended_y_limits(1), extended_y_limits(2), extended_y_limits(2), extended_y_limits(1)], ...
             'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end
end
hold off
xlim([t(1), t(end)])
ylim(extended_y_limits) % Imposta i nuovi limiti Y
title('Filtered signal with random cut forest anomalies')

% Subplot 2: Tukey anomalies
subplot(3, 1, 2)
plot(t, filtered_signal, 'b')
hold on
y_limits = ylim; % Ottieni i limiti dell'asse Y
range_y = y_limits(2) - y_limits(1);
extended_y_limits = [y_limits(1) - margin_factor * range_y, y_limits(2) + margin_factor * range_y];
if ~isempty(an_edges_filt_tot)
    for i = 1:2:length(an_edges_filt_tot)-1
        % Aggiunge una banda rossa per l'anomalia
        fill([t(an_edges_filt_tot(i)), t(an_edges_filt_tot(i)), t(an_edges_filt_tot(i+1)), t(an_edges_filt_tot(i+1))], ...
             [extended_y_limits(1), extended_y_limits(2), extended_y_limits(2), extended_y_limits(1)], ...
             'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end
end
hold off
xlim([t(1), t(end)])
ylim(extended_y_limits) % Imposta i nuovi limiti Y
title('Filtered signal with Tukey anomalies')

% Subplot 3: Selected anomalies
subplot(3, 1, 3)
plot(t, filtered_signal, 'b')
hold on
y_limits = ylim; % Ottieni i limiti dell'asse Y
range_y = y_limits(2) - y_limits(1);
extended_y_limits = [y_limits(1) - margin_factor * range_y, y_limits(2) + margin_factor * range_y];
if ~isempty(an_edges_filt_opt)
    for i = 1:2:length(an_edges_filt_opt)-1
        % Aggiunge una banda rossa per l'anomalia
        fill([t(an_edges_filt_opt(i)), t(an_edges_filt_opt(i)), t(an_edges_filt_opt(i+1)), t(an_edges_filt_opt(i+1))], ...
             [extended_y_limits(1), extended_y_limits(2), extended_y_limits(2), extended_y_limits(1)], ...
             'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end
end
hold off
xlim([t(1), t(end)])
ylim(extended_y_limits) % Imposta i nuovi limiti Y
title('Filtered signal with selected anomalies')

%}



 %{

figure;

subplot(3, 1, 1)
plot(t, filtered_signal, 'b')
hold on
if ~isempty(an_edges_filt_sc)
    for i = 1:size(an_edges_filt_sc, 1)
        fill([t(an_edges_filt_sc(i, 1)), t(an_edges_filt_sc(i, 2)), t(an_edges_filt_sc(i, 2)), t(an_edges_filt_sc(i, 1))], ...
             [min(filtered_signal), min(filtered_signal), max(filtered_signal), max(filtered_signal)], ...
             'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end
end
hold off
xlim([0, t(end)])
title('Filtered signal with random cut forest anomalies')

subplot(3, 1, 2)
plot(t, filtered_signal, 'b')
hold on
if ~isempty(an_edges_filt_tot)
    for i = 1:size(an_edges_filt_tot, 1)
        fill([t(an_edges_filt_tot(i, 1)), t(an_edges_filt_tot(i, 2)), t(an_edges_filt_tot(i, 2)), t(an_edges_filt_tot(i, 1))], ...
             [min(filtered_signal), min(filtered_signal), max(filtered_signal), max(filtered_signal)], ...
             'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end
end
hold off
xlim([0, t(end)])
title('Filtered signal with Tukey anomalies')

subplot(3, 1, 3)
plot(t, filtered_signal, 'b')
hold on
if ~isempty(an_edges_filt_opt)
    for i = 1:size(an_edges_filt_opt, 1)
        fill([t(an_edges_filt_opt(i, 1)), t(an_edges_filt_opt(i, 2)), t(an_edges_filt_opt(i, 2)), t(an_edges_filt_opt(i, 1))], ...
             [min(filtered_signal), min(filtered_signal), max(filtered_signal), max(filtered_signal)], ...
             'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end
end
hold off
xlim([0, t(end)])
title('Filtered signal with selected anomalies')

%}   









end