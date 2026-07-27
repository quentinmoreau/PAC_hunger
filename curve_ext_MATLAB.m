function [FreqOut, EnergyOut] = curve_ext_MATLAB(TFR, fs, lambda)
    % CURVE_EXT - Estrae una curva massima in energia e minima in curvatura
    % Input:
    % Txr, Txi - Parte reale e immaginaria del TFR (matrici na x N)
    % fs - Log2 delle frequenze (vettore 1 x na)
    % lambda - Parametro di regolarizzazione >= 0
    % Output:
    % FreqOut - Indici della curva (vettore 1 x N)
    % EnergyOut - Energia minima della curva trovata (scalare)

    % Costanti
    epsVal = 1e-8;

    Txr=real(TFR);
    Txi=imag(TFR);

    % Dimensioni
    [na, N] = size(TFR);

    % Calcola l'energia (somma dei quadrati delle parti reale e immaginaria)
    Energy = zeros(na, N);
    if ~isempty(Txr)
        Energy = Energy + Txr.^2;
    end
    if ~isempty(Txi)
        Energy = Energy + Txi.^2;
    end
    sumEnergy = sum(Energy, 'all');
    Energy = -log(Energy / sumEnergy + epsVal); % Converti in energia negativa logaritmica

    % Inizializza le matrici per la programmazione dinamica
    FVal = inf(na, N); % Costo cumulativo
    FVal(:, 1) = Energy(:, 1); % Costi iniziali
    PrevIdx = zeros(na, N); % Indici precedenti per tracciare il percorso

    % Programmazione dinamica: calcolo del costo minimo
    for i = 2:N
        for j = 1:na
            % Calcola il costo minimo di transizione per il tempo i
            dist = (fs(j) - fs).^2; % Distanza al quadrato
            transCost = FVal(:, i-1) + lambda * dist';
            [minCost, minIdx] = min(transCost);
            FVal(j, i) = minCost + Energy(j, i); % Aggiorna il costo cumulativo
            PrevIdx(j, i) = minIdx; % Salva l'indice del percorso minimo
        end
    end

    % Trova la posizione della minima energia totale nell'ultimo passo
    [EnergyOut, lastIdx] = min(FVal(:, N));

    % Backtracking per ottenere il percorso della curva
    FreqOut = zeros(1, N);
    FreqOut(N) = lastIdx;
    for i = N-1:-1:1
        FreqOut(i) = PrevIdx(FreqOut(i+1), i+1);
    end

    % MATLAB usa indici da 1, quindi aggiungiamo 1
    %FreqOut = FreqOut + 1;
end


