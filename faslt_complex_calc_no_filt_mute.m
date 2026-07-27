function [faslt, buffer] = faslt_complex_calc_no_filt_mute(x, fs, freq_faslt, padding)


wrong_pad_flag=1;
while wrong_pad_flag==1

  if ~strcmp(padding, 'predictive') && ~strcmp(padding, 'zero') && ~strcmp(padding, 'none')
       disp('Wrong padding scheme, please select predictive or zero');
       padding = input('Write down your choice: ');
  else
       wrong_pad_flag=0;
  end
end

if size(x, 1)>size(x, 2)
    x=x';
end

 %{
L=length(x);

    %26.01389 è il tf supp0rt per superlet a 5 cicli con f_0=0.016 e
    %epsilon=0.001 
    N_e=2^(nextpow2(length(x)+2*26.0139*fs));
    N_pad=floor((N_e-length(x))/2);

 


if strcmp(padding, 'predictive')


    % su algoritmo originale,   %w=2.^(-(L/fs-(1:L)/fs)/(wp.t2h-wp.t1h));
    % wp.t2h si intende qui il \tau(0.5) cioè per \epsilon=0.5 e ottengo
    % 5.3323
    %w=2.^(-(L/fs-(1:L)/fs)/(26.01389*2)); %weigths questo è per eps=10^-3
    w=2.^(-(L/fs-(1:L)/fs)/(5.3323*2)); %weigths  

 fprintf('padding to the right - ');

 f_sig_right = predict(x,fs,floor(N_pad),[0.016 0.16],L, w, 'on'); %forecasting signal
 L_right= length(f_sig_right);

 disp(' ');

 %padding to the left
 fprintf('padding to the left - ');
 
 f_sig_left = predict(fliplr(x),fs,floor(N_pad),[0.016 0.16],L, w, 'on'); %forecasting signal
 f_sig_left=fliplr(f_sig_left);
 L_left= length(f_sig_left);

 x_pad=[f_sig_left x f_sig_right];

 figure;
 plot([zeros(1, L_left) x zeros(1, L_right)])
 hold on
 plot(x_pad)

 figure;
 plot(w)

elseif strcmp(padding, 'zero');
    
    right_pad=zeros(1, floor(N_pad));
    left_pad=zeros(1, floor(N_pad));
    L_right= length(right_pad);
    L_left= length(left_pad);
    
    

    x_pad=[left_pad x right_pad];

elseif strcmp(padding, 'none');

    x_pad=x;
    L_right=0;
    L_left=0;
    figure;
    plot(x_pad)
    title('no padding')

end

%}

x_pad=x;
    

%disp(' ');
 %disp('Calculating Fractional adaptive superlet transform');   
%[faslt] = faslt_complex(x_pad, fs,freq_faslt, 3, [3 5], 1);
[faslt, buffer] = faslt_complex_predictive(x_pad, fs,freq_faslt, 3, [3 5], 1, padding);
%[faslt, buffer] = faslt_complex_predictive(x_pad, fs,freq_faslt, 2, [8 8], 0, padding);
%faslt_complex_predictive

%faslt = faslt(:, L_left+1:end-L_right);




% --------------------------  auxiliar functions


% computes the complex Morlet wavelet for the desired center frequency Fc
% with Nc cycles, with a sampling frequency Fs.
function w = cxmorlet(Fc, Nc, Fs)
    %we want to have the last peak at 2.5 SD
    sd  = (Nc / 2) * (1 / Fc) / 2.5;
    wl  = 2 * floor(fix(6 * sd * Fs)/2) + 1;
    w   = zeros(wl, 1);
    gi  = 0;
    off = fix(wl / 2);
    
    for i = 1 : wl
        t       = (i - 1 - off) / Fs;
        w(i)    = bw_cf(t, sd, Fc);
        gi      = gi + gauss(t, sd);
    end
    
    w = w ./ gi;
return

% compute the complex wavelet coefficients for the desired time point t,
% bandwidth bw and center frequency cf
function res = bw_cf(t, bw, cf)
    cnorm   = 1 / (bw * sqrt(2 * pi));
    exp1    = cnorm * exp(-(t^2) / (2 * bw^2));
    res     = exp(2i * pi * cf * t) * exp1;
return;

% compute the gaussian coefficient for the desired time point t and
% standard deviation sd
function res = gauss(t, sd)
    cnorm   = 1 / (sd * sqrt(2 * pi));
    res     = cnorm * exp(-(t^2) / (2 * sd^2));
return;

% tell me if a number is an integer or a fractional
function res = is_fractional(x)
    res = fix(x) ~= x;
return;