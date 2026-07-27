function modified_vector = merge_short_gaps(binary_vector, L)
    % Funzione per unire gruppi di 1 separati da gruppi di 0 più piccoli di L
    % Input:
    %   binary_vector: vettore di 0 e 1
    %   L: lunghezza massima del gruppo di 0 da trasformare in 1
    % Output:
    %   modified_vector: vettore binario modificato

    % Trova le posizioni degli zeri e degli uni
    binary_vector = binary_vector(:); % Assicurati che sia un vettore colonna
    diff_vector = [diff(binary_vector)];
    start_zeros = find(diff_vector == - 1)+1; % Inizio dei gruppi di zeri
    end_zeros = find(diff_vector == 1); % Fine dei gruppi di zeri

    % Gestisci il caso in cui il gruppo di zeri sia all'inizio o alla fine
    if binary_vector(1) == 0
        start_zeros = [1; start_zeros];
    end
    if binary_vector(end) == 0
        end_zeros = [end_zeros; length(binary_vector)];
    end

    % Unisci i gruppi di 1 separati da gruppi di 0 corti
    for i = 1:length(start_zeros)
        if (end_zeros(i) - start_zeros(i) + 1) <= L
            binary_vector(start_zeros(i):end_zeros(i)) = 1;
        end
    end

    modified_vector = binary_vector;
end
