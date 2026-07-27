function [MDisp, RDE] = MvDispEn(Data, varargin)
% MvDispEn  estimates the multivariate dispersion entropy of a multivariate dataset.
%
%   [MDisp, RDE] = MvDispEn(Data)
%
%   Returns the multivariate dispersion entropy estimate (``MDisp``) and
%   the reverse dispersion entropy (``RDE``) for the M multivariate sequences 
%   in ``Data`` using the default parameters:
%   embedding dimension = 2*ones(M,1), time delay = ones(M,1), # symbols = 3, 
%   algorithm method = 'v1' (see below), data transform = normalised cumulative density function (ncdf)
%   logarithm = natural, entropy normalization = true,
%  
%   .. important::
%      By default, ``MvDispEn`` uses the method termed ``mvDEii`` in [1],
%      which follows the original multivariate embedding approach of Ahmed & Mandic [2].
%      The ``v1`` method therefore returns a singular entropy estimate.
%    
%      If the ``v2`` method is selected (``Methodx=='v2'``), the main method
%      outlined in [1] termed ``mvDE`` is applied. In this case, entropy is estimated
%      using each combination of multivariate delay vectors with lengths 1:max(m),
%      with each entropy value returned accordingly. See [1] for more info.
%
%   [MDisp, RDE] = MvDispEn(Data, name, value, ...)
%
%   Returns the multivariate dispersion entropy estimate (``MDisp``) for the M
%   multivariate data sequences in ``Data`` using the specified name/value pair arguments:
%
%       :Data:     - Multivariate dataset, NxM matrix of N (>10) observations (rows) and M (cols) univariate data sequences
%       :m:        - Embedding Dimension, a vector of M positive integers
%       :tau:      - Time Delay, a vector of M positive integers
%       :c:        - Number of symbols in transform, an integer > 1
%       :Methodx:  - The method of multivariate dispersion entropy estimation as outlined in [1], either
%              :'v1': - employs the method consistent with the  original multivariate embedding approach of Ahmed & Mandic [2], termed ``mvDEii`` in [1]. (default)
%              :'v2': - employs the main method derived in [1],  termed ``mvDE``.
%       :Typex:    - Type of data-to-symbolic sequence transform, one  of the following {``linear``, ``kmeans``, ``ncdf``, ``equal``}
%                   See the `EntropyHub Guide` for more info on these transforms.
%       :Logx:  - Logarithm base, a positive scalar
%       :Norm:  - Normalisation of ``MDisp`` and ``RDE`` values, a boolean:
%                :false:   no normalisation (default)
%                :true:    normalises w.r.t number of possible dispersion patterns (``c^m``).
%
%   See also:
%       DispEn, DispEn2D, MvSampEn, MvFuzzEn, MvPermEn, MvMSEn, cMvMSEn,
%
%   References:
%       [1] H Azami, A Fernández, J Escudero
%           "Multivariate Multiscale Dispersion Entropy of Biomedical Times Series"
%           Entropy 2019, 21, 913.
% 
%       [2] Ahmed Mosabber Uddin, Danilo P. Mandic
%           "Multivariate multiscale entropy: A tool for complexity
%           analysis of multichannel data."
%           Physical Review E 84.6 (2011): 061918.
%
%       [3] Mostafa Rostaghi and Hamed Azami,
%            "Dispersion entropy: A measure for time-series analysis." 
%            IEEE Signal Processing Letters 
%            23.5 (2016): 610-614.
% 
%       [4] Hamed Azami and Javier Escudero,
%            "Amplitude-and fluctuation-based dispersion entropy." 
%            Entropy 
%            20.3 (2018): 210.
% 
%       [5] Li Yuxing, Xiang Gao and Long Wang,
%            "Reverse dispersion entropy: A new complexity measure for sensor signal." 
%            Sensors 
%            19.23 (2019): 5203.
% 

narginchk(1,15)
p = inputParser;
Data = squeeze(Data);
Dn = size(Data,2);
N = size(Data,1);
Chk = @(x) isnumeric(x) && isvector(x) && (length(x)==Dn) && (min(x)>0) && all(mod(x,1)==0);
Chk2 = @(x) isnumeric(x) && isscalar(x) && (x > 0);
Chk3 = ["linear","kmeans","ncdf","equal"];

addRequired(p,'Data',@(x) isnumeric(x) && ismatrix(x) && (Dn>1) && N>10);
addParameter(p,'m',2*ones(Dn,1),Chk);
addParameter(p,'tau',ones(Dn,1),Chk);
addParameter(p,'c',3, @(x) isscalar(x) && isnumeric(x) && (x>1) && (x<100) && (mod(x,1)==0));
addParameter(p,'Methodx','v1', @(x) (isstring(x) || ischar(x)) && any(strcmpi(string(x),["v1","v2"])));
addParameter(p,'Logx',exp(1),Chk2);
addParameter(p,'Norm',false,@(x) islogical(x));
addParameter(p,'Typex','ncdf',@(x) (isstring(x) || ischar(x)) && any(strcmpi(string(x),Chk3)));
parse(p,Data,varargin{:})
m = p.Results.m(:); tau = p.Results.tau(:); c = p.Results.c;
Typex = p.Results.Typex; Logx = p.Results.Logx; Norm = p.Results.Norm;
Methodx = p.Results.Methodx;

Sx = zeros(size(Data), 'int8');
for q = 1:Dn
    Sig = Data(:,q);

    switch lower(Typex)
        case 'linear'
            Zi = discretize(Sig, linspace(min(Sig),max(Sig),c+1));

        case 'kmeans'
            [Zx,Clux] = kmeans(Sig, c, 'MaxIter', 200);
            [~,xx] = sort(Clux);        Zi = zeros(1,N);
            for k = 1:c;  Zi(Zx==xx(k)) = k;    end
            clear Clux Zx xx

        case 'ncdf'
            Zx = normcdf(Sig,mean(Sig),std(Sig,1));
            Zi = discretize(Zx,linspace(0,1,c+1));

        case 'equal'
            [~,ix] = sort(Sig);
            xx = round(linspace(0,N,c+1),'Tiebreaker','even');
            Zi = zeros(1,N);
            for k = 1:c
                Zi(ix(xx(k)+1:xx(k+1))) = k;
            end
            clear ix xx
    end

    Sx(:,q) = Zi;
end

Nx = N - max((m-1).*tau);
Vex = zeros(Nx,sum(m), 'int8');
q = 1;
for k = 1:Dn
    for p=1:m(k)
        Vex(:,q) = Sx(1+(p-1)*tau(k):Nx+(p-1)*tau(k),  k);
        q = q+ 1;
    end
end

if strcmpi(Methodx, 'v1')
    Px = unique(Vex, 'rows');
    Counter = zeros(1, size(Px,1));
    for n = 1:size(Px,1)
        Counter(n) = sum(~any(Vex - Px(n,:),2));
    end
    Counter(Counter==0)=[];
    Ppi = Counter/sum(Counter);
    assert(round(sum(Ppi),5)==1, 'Potential error with probability calculation')

    MDisp = -sum(Ppi.*(log(Ppi)/log(Logx)));
    RDE = sum((Ppi - (1/(c^sum(m)))).^2);
    if Norm
        MDisp = MDisp/(log(c^sum(m))/log(Logx));
        RDE = RDE/(1 - (1/(c^sum(m))));
    end

elseif strcmpi(Methodx, 'v2')
    P = sum(m);
    for k = 1:max(m)
        fprintf('. ')
        qi = 1;
        Vez = zeros(Nx*nchoosek(P,k), k, 'int8');
        for q = nchoosek(1:P,k)'
            Vez(1+(qi-1)*Nx:qi*Nx,:) = Vex(:,q);
            qi = qi+1;
        end

        Px = unique(Vez, 'rows');
        for n = 1:size(Px,1); Counter(n) = sum(~any(Vez - Px(n,:),2)); end
   
        Counter(Counter==0)=[];
        Ppi = Counter/sum(Counter);
        assert(round(sum(Ppi),5)==1, 'Potential error with probability calculation')

        MDisp(k) = -sum(Ppi.*(log(Ppi)/log(Logx)));
        RDE(k) = sum((Ppi - (1/(c^k))).^2);    
        if Norm
            MDisp(k) = MDisp(k)/(log(c^k)/log(Logx));
            RDE(k) = RDE(k)/(1 - (1/(c^k)));
        end
    end
end

end


%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub