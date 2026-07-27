function filtered_signal = cheb_EGG_filt(signal, fs)
    % keep_frequency_range Mantiene solo un intervallo di frequenze specifico in un segnale.
    % INPUTS:
    %   signal        - Vettore contenente il segnale di input.
    %   Fs            - Frequenza di campionamento del segnale in Hz.
    %
    % OUTPUT:
    %   filtered_signal - Segnale filtrato, mantenendo solo l'intervallo di frequenze specificato.

    low=0.18;
    high=0.014;
    deg1=10;
    deg2=25;
    deg11=6;
    deg22=22;

    [b,a]=cheby2(deg1,deg2,low*2/fs,'low');
    x1=filtfilt(b,a,signal);


    [b,a] = cheby2(deg11,deg22,high*2/fs,'high');
    filtered_signal=filtfilt(b,a,x1);



end
