function anomaly_vector = interval_to_vector(artifact_intervals, signal_length)
    % Funzione per convertire un vettore di intervalli in un vettore binario di anomalie.
    %
    % artifact_intervals: vettore con coppie di valori consecutivi che rappresentano [inizio, fine]
    % signal_length: lunghezza totale del segnale
    %
    % Restituisce:
    % anomaly_vector: vettore binario di lunghezza signal_length con 1 per anomalie e 0 per normalità
    
    % Inizializza il vettore di anomalie a zero
    anomaly_vector = zeros(signal_length, 1);

    % Assicura che il numero di elementi in artifact_intervals sia pari
    if mod(length(artifact_intervals), 2) ~= 0
        error('artifact_intervals deve contenere un numero pari di valori (inizio e fine per ogni intervallo).');
    end

    % Imposta a 1 i valori negli intervalli di artefatti
    for i = 1:2:length(artifact_intervals)
        start_idx = artifact_intervals(i);
        end_idx = artifact_intervals(i + 1);
        
        % Assicurati che gli indici siano entro i limiti
        start_idx = max(1, start_idx);
        end_idx = min(signal_length, end_idx);
        
        % Imposta l'intervallo a 1
        anomaly_vector(start_idx:end_idx) = 1;
    end
end