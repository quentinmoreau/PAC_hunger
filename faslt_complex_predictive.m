%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   FRACTIONAL ADAPTIVE SUPERRESOLUTION WAVELET (SUPERLET) TRANSFORM 
% 
%   AUTHOR:         Harald Bârzan
%   DATE:           April 2019
%   DESCRIPTION:
%
%   Computes the adaptive superresolution wavelet (superlet) transform on 
%   input data to produce a time-frequency representation. For each 
%   frequency of interest, the closest integer order from the order 
%   interval will be chosen to produce each superlet. A superlet is a set 
%   of wavelets with the same center frequency but different number of 
%   cycles.
%
%   REFERENCE:
%   
%   Time-frequency super-resolution with superlets
%   Moca, V.V., Nagy-Dãbâcan, A., Bârzan, H., Mure?an, R.C.
%   https://www.nature.com/articles/s41467-020-20539-9
%   
%   NOTES:
%
%   If the input data consists of multiple buffers, a wavelet spectrum will
%   be computed for each of the buffers and averaged to produce the final 
%   result.
%   If the order parameter (ord) is empty, this function will return the
%   standard CWT (one wavelet per frequency of interest).
%
%   INPUT:
%   > input         - [buffers x samples] matrix
%   > Fs            - sampling frequency in Hz
%   > F             - frequency-of-interest buffer
%   > c1            - number of initial wavelet cycles
%   > o             - [1 x 2] interval of superresolution orders (optional)
%   > mult          - specifies the use of multiplicative superresolution
%                     (0 - additive, != 0 - multiplicative)
%
%   OUTPUT:
%   > wtresult      - [frequencies x samples] superlet spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wtresult, buffer] = faslt_complex_predictive(input, Fs, F, c1, o, mult, padding)

% check frequency of interest parameter
if (isempty(F))
    error('frequencies not defined'); 
end

% check order parameter and initialize the order used at each frequency. 
% if empty, go with an order of 1 for each frequency (single wavelet per
% set)
if (~isempty(o))
    % fractional order is whatever the user put in
    order_frac  = linspace(o(1), o(2), numel(F));
    % integer order is the real number of wavelets in the set
    order_int   = ceil(order_frac);
else
    order_frac  = ones(numel(F), 1);
    order_int   = order_frac;
end

% validate input buffer
if (isempty(input))
    error('input is empty'); 
end

% if input is a column vector, turn it into a row vector instead
if (size(input, 2) == 1 && size(input, 1) > 1)
    input = input'; 
end

% get the input size
[Nbuffers, Npoints] = size(input);

% the padding will be size of the lateral zero-pads, which serve to avoid
% border effects during convolution
padding_size = 0;
weights_size=Inf;

% the wavelet sets
wavelets = cell(numel(F), max(order_int));
        
% initialize wavelet sets for either additive or multiplicative
% superresolution

% Inizializza variabili per memorizzare la wavelet più lunga e gli indici
max_wavelet_len = Inf;        % Per tenere traccia della lunghezza massima della wavelet
longest_wavelet = [];       % Per memorizzare la wavelet più lunga
longest_i_freq = 0;         % Indice della frequenza per la wavelet più lunga
longest_i_ord = 0;          % Indice dell'ordine per la wavelet più lunga

%size(order_int)
for i_freq = 1 : numel(F)
    for i_ord = 1 : order_int(i_freq)
        
        % compute the number of cycles (additive or multiplicative)
        if (mult ~= 0)
            n_cyc = i_ord * c1;
        else
            n_cyc = i_ord + c1;
        end
        
        % add the wavelet to the set
        wavelets{i_freq, i_ord} = cxmorlet(F(i_freq), n_cyc, Fs);
        
        % the margin will be the half-size of the largest wavelet
        padding_size = max(padding_size, fix(numel(wavelets{i_freq, i_ord}) / 2));
        weights_size = min(weights_size, fix(numel(wavelets{i_freq, i_ord}) / 2));
        %padding=0;

         % Verifica la lunghezza della wavelet
        wavelet_len = numel(wavelets{i_freq, i_ord});
        
        % Se la lunghezza della wavelet è la più lunga, aggiorna le variabili
        if wavelet_len < max_wavelet_len
            max_wavelet_len = wavelet_len;
            longest_wavelet = wavelets{i_freq, i_ord};
            longest_i_freq = i_freq;
            longest_i_ord = i_ord;
        end
    end
end

%weights_size
%size(input)

%figure;
%plot(real(longest_wavelet))
%hold on
%plot(abs(longest_wavelet))

%longest_i_freq
%longest_i_ord


% the zero-padded buffer
buffer = zeros(Npoints + 2 * padding_size, 1);
%padding_size

if strcmp(padding, 'predictive')

L=length(input(1, :));
%weights_size

[T1_half, T2_half] = calculate_wavelet_time_supports(longest_wavelet, F(longest_i_freq), 0.5, longest_i_ord*c1, Fs);
 w=exp(-2*(L/Fs-(1:L)/Fs)/(T2_half-T1_half)); %weigths  
  %w=2.^(-(L/Fs-(1:L)/Fs)/(weights_size)); %weigths  
   %w=2.^(-(L/Fs-(1:L)/Fs)/(10));
 %T2_half-T1_half
 %padding_size


 figure;
 plot(w)

 fprintf('padding to the right - ');

 %[T1_eps, T2_eps] = calculate_wavelet_time_supports(longest_wavelet, F(longest_i_freq), 0.001, longest_i_ord*c1, Fs)
 %T2_eps-T1_eps
 f_sig_right = predict(input(1, :),Fs,floor(padding_size),[0.016 0.16],L, w, 'on'); %forecasting signal
 %L_right= length(f_sig_right);

 %disp(' ');
 fprintf('\n');

 %padding to the left
 fprintf('padding to the left - ');
 
 f_sig_left = predict(fliplr(input(1, :)),Fs,floor(padding_size),[0.016 0.16],L, w, 'on'); %forecasting signal
 f_sig_left=fliplr(f_sig_left);
 %L_left= length(f_sig_left);

 fprintf('\n');


% the output scalogram
wtresult = zeros(numel(F), Npoints);

% convenience indexers for the zero-padded buffer
bufbegin    = padding_size + 1;
bufend      = padding_size + Npoints;

buffer(1:bufbegin-1)=f_sig_left;
buffer(bufend+1 : end)=f_sig_right;
figure;
plot(buffer)

%buffer = zeros(Npoints + 2 * padding_size, 1);

end

if strcmp(padding, 'zero')

% convenience indexers for the zero-padded buffer
bufbegin    = padding_size + 1;
bufend      = padding_size + Npoints;

% the output scalogram
wtresult = zeros(numel(F), Npoints);

buffer = zeros(Npoints + 2 * padding_size, 1);
end

% loop over the input buffers
for i_buf = 1 : Nbuffers
    for i_freq = 1 : numel(F)
        % pooling buffer, starts with 1 because we're doing geometric mean
        temp = ones(1, Npoints);
        
        % fill the central part of the buffer with input data
        buffer(bufbegin : bufend) = input(i_buf, :);
        
        % get the number of integer wavelets
        n_wavelets = floor(order_frac(i_freq));
        
        % compute the convolution of the buffer with each wavelet in the
        % current set (integer wavelets)
        for i_ord = 1 : n_wavelets
            % restricted convolution (input size == output size)
            tempcx = conv(buffer, wavelets{i_freq, i_ord}, 'same');
            
            % accumulate the magnitude (times 2 to get the full spectral
            % energy), pool with exponent = 1
            temp = temp .*tempcx(bufbegin : bufend)';
        end
        
        % handle fractional exponent
        if (is_fractional(order_frac(i_freq)) && ...
            ~isempty(wavelets{i_freq, order_int(i_freq)}))
            % set the order index
            i_ord = order_int(i_freq);
            
            % the exponent is the fractional remainder
            exponent = order_frac(i_freq) - fix(order_frac(i_freq));
            
             % restricted convolution (input size == output size)
            tempcx = conv(buffer, wavelets{i_freq, i_ord}, 'same');
            
            % accumulate the magnitude (times 2 to get the full spectral
            % energy), pool with exponent = 1
            temp = temp .*tempcx(bufbegin : bufend)'.^ exponent;
        end
            
        % compute the order of the geometric mean
        root = 1 / order_frac(i_freq);
        temp = temp .^ root;
        
        % accumulate the current FOI to the result spectrum
        wtresult(i_freq, :) = wtresult(i_freq, :) + temp;
    end
end

% scale the output by the number of input buffers
wtresult = wtresult ./ Nbuffers;

return


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