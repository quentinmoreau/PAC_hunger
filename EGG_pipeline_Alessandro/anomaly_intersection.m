function r1_updated = anomaly_intersection(r1, r_window_sc)
    % Verifica che i due vettori abbiano la stessa lunghezza
    if length(r1) ~= length(r_window_sc)
        error('I vettori devono avere la stessa lunghezza.');
    end

    r1=r1(:)';
    
    % Trova gli indici di inizio e fine dei gruppi di 1 consecutivi in r1
    d_r1 = diff([0, r1, 0]);
    start_indices = find(d_r1 == 1);
    end_indices = find(d_r1 == -1) - 1;
    
    % Inizializza r1_updated come copia di r1
    r1_updated = r1;
    
    % Loop su ciascun gruppo di 1 consecutivi
    for i = 1:length(start_indices)
        start_idx = start_indices(i);
        end_idx = end_indices(i);
        
        % Controlla se c'è un'intersezione con r_window_sc
        if ~any(r_window_sc(start_idx:end_idx) == 1)
            % Se non c'è intersezione, rimpiazza il gruppo con 0
            r1_updated(start_idx:end_idx) = 0;
        end
    end
end