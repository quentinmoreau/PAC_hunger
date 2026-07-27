function [faslt, buffer] = faslt_complex_calc_no_filt(x, fs, freq_faslt, padding)
% Calculates the Fractional Adaptive Superlet Transform (FASLT)
% without additional filtering or zero/predictive padding.
%
% INPUTS:
%   x          - input signal (1D vector)
%   fs         - sampling frequency
%   freq_faslt - frequencies for FASLT
%   padding    - 'predictive', 'zero', or 'none'
%
% OUTPUTS:
%   faslt  - FASLT transform
%   buffer - auxiliary buffer from FASLT calculation

%% Validate padding input
wrong_pad_flag = 1;
while wrong_pad_flag == 1
    if ~strcmp(padding, 'predictive') && ~strcmp(padding, 'zero') && ~strcmp(padding, 'none')
        disp('Wrong padding scheme, please select predictive, zero, or none');
        padding = input('Write down your choice: ');
    else
        wrong_pad_flag = 0;
    end
end

%% Ensure signal is a row vector
if size(x, 1) > size(x, 2)
    x = x';
end

%% -- Padding (currently commented out) --
%{
% Original code for predictive/zero/none padding (not active)
% L = length(x);
% N_e = 2^(nextpow2(length(x) + 2*26.0139*fs));
% N_pad = floor((N_e - length(x)) / 2);
% if strcmp(padding, 'predictive')
%     % Predictive padding code
% elseif strcmp(padding, 'zero')
%     % Zero padding code
% elseif strcmp(padding, 'none')
%     % No padding
% end
%}

x_pad = x; % currently no padding applied

%% Calculate Fractional Adaptive Superlet Transform
disp('Calculating Fractional Adaptive Superlet Transform...');

[faslt, buffer] = faslt_complex_predictive(...
    x_pad, fs, freq_faslt, 3, [3 5], 1, padding);

% If padding were applied:
% faslt = faslt(:, L_left+1:end-L_right);

end

%% -------------------------- Auxiliary functions --------------------------

function w = cxmorlet(Fc, Nc, Fs)
% Compute complex Morlet wavelet for center frequency Fc
% with Nc cycles at sampling frequency Fs
    sd  = (Nc / 2) * (1 / Fc) / 2.5; 
    wl  = 2 * floor(fix(6 * sd * Fs)/2) + 1;
    w   = zeros(wl, 1);
    gi  = 0;
    off = fix(wl / 2);
    
    for i = 1 : wl
        t     = (i - 1 - off) / Fs;
        w(i)  = bw_cf(t, sd, Fc);
        gi    = gi + gauss(t, sd);
    end
    
    w = w ./ gi;
end

function res = bw_cf(t, bw, cf)
% Compute complex wavelet coefficient for time t, bandwidth bw, and center frequency cf
    cnorm = 1 / (bw * sqrt(2 * pi));
    exp1  = cnorm * exp(-(t^2) / (2 * bw^2));
    res   = exp(2i * pi * cf * t) * exp1;
end

function res = gauss(t, sd)
% Compute Gaussian coefficient for time t with standard deviation sd
    cnorm = 1 / (sd * sqrt(2 * pi));
    res   = cnorm * exp(-(t^2) / (2 * sd^2));
end

function res = is_fractional(x)
% Returns true if x is fractional, false if integer
    res = fix(x) ~= x;
end