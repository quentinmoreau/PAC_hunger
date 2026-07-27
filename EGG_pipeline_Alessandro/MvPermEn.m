function [MPerm, MPnorm] = MvPermEn(Data, varargin)
% MvPermEn  estimates the multivariate permutation entropy of a multivariate dataset.
%
%   [MPerm MPnorm] = MvPermEn(Data)
%
%   Returns the multivariate permutation entropy estimate (``MPerm``) and
%   the normalized permutation entropy for the M multivariate sequences in
%   ``Data`` using the default parameters:
%   embedding dimension = 2*ones(M,1), time delay = ones(M,1), 
%   logarithm = 2, normalisation = w.r.t #symbols (sum(``m-1``))
%
%   .. attention::
%         The multivariate permutation entropy algorithm implemented here uses
%         multivariate embedding based on Takens' embedding theorem, and follows
%         the methods for multivariate entropy estimation through shared spatial 
%         reconstruction as originally presented by Ahmed & Mandic [1]. 
%        
%         This function does **NOT** use the multivariate permutation entropy 
%         algorithm of Morabito et al. (Entropy, 2012) where the entropy values of 
%         individual univariate sequences are averaged because such methods do not
%         follow the definition of multivariate embedding and therefore do not
%         consider cross-channel statistical complexity.
%        
%   .. note::
%         To maximize the number of points in the embedding process, this
%         algorithm uses N-max(tau*m) delay vectors and **not** N-max(m)*max(tau)
%         as employed in [1].
% 
%   [MPerm, MPnorm] = MvPermEn(Data, name, value, ...)
%
%   Returns the multivariate permutation entropy estimate (``MPerm``) for
%   the M multivariate data sequences in ``Data`` using the specified name/value pair arguments:
%
%       :Data:  - Multivariate dataset, NxM matrix of N (>10) observations (rows) and M (cols) univariate data sequences
%       :m:     - Embedding Dimension, a vector of M positive integers
%       :tau:   - Time Delay, a vector of M positive integers
%       :Typex: - Permutation entropy variation, can be one of the following strings: {``'modified'``, ``'ampaware'``, ``'weighted'``, ``'edge'``, ``'phase'``}
%                See the `EntropyHub guide <https://github.com/MattWillFlood/EntropyHub/blob/main/EntropyHub%20Guide.pdf>`_ for more info on PermEn variants.    
%       :tpx:   - Tuning parameter for associated permutation entropy variation.
%              :ampaware:  ``tpx`` is the A parameter, a value in range [0 1] (default = 0.5)
%              :edge:      ``tpx`` is the r sensitivity parameter, a scalar > 0 (default = 1)
%              :phase:     ``tpx`` is the option to unwrap the phase angle of Hilbert-transformed signal, either [] or 1 (default = 0)
%       :Norm:  - Normalisation of MPnorm value, a boolean operator:
%              :false: -  normalises w.r.t log(# of permutation symbols [sum(m)-1]) - default
%              :true:  -  normalises w.r.t log(# of all possible permutations [sum(m)!])
%       :Logx:  - Logarithm base, a positive scalar (defualt = 2; enter 0 for  natural logarithm)
%
%   See also:
%       PermEn, PermEn2D, XPermEn, MvMSEn, MvFuzzEn, MvSampEn, MvDispEn, cMvMSEn.
%
%   References:
%     [1] Ahmed Mosabber Uddin, Danilo P. Mandic
%           "Multivariate multiscale entropy: A tool for complexity
%           analysis of multichannel data."
%           Physical Review E 84.6 (2011): 061918.
%
%     [2] Christoph Bandt and Bernd Pompe, 
%           "Permutation entropy: A natural complexity measure for time series." 
%           Physical Review Letters,
%           88.17 (2002): 174102.
% 
%     [3] Chunhua Bian, et al.,
%           "Modified permutation-entropy analysis of heartbeat dynamics."
%           Physical Review E
%           85.2 (2012) : 021906
% 
%     [4] Bilal Fadlallah, et al.,
%           "Weighted-permutation entropy: A complexity measure for time 
%           series incorporating amplitude information." 
%           Physical Review E 
%           87.2 (2013): 022911.
% 
%     [5] Hamed Azami and Javier Escudero,
%           "Amplitude-aware permutation entropy: Illustration in spike 
%           detection and signal segmentation." 
%           Computer methods and programs in biomedicine,
%           128 (2016): 40-51.
% 
%     [6] Zhiqiang Huo, et al.,
%           "Edge Permutation Entropy: An Improved Entropy Measure for 
%           Time-Series Analysis," 
%           45th Annual Conference of the IEEE Industrial Electronics Soc,
%           (2019), 5998-6003
% 
%     [7] Maik Riedl, Andreas Müller, and Niels Wessel,
%           "Practical considerations of permutation entropy." 
%           The European Physical Journal Special Topics 
%           222.2 (2013): 249-262.
%
%     [8] Kang Huan, Xiaofeng Zhang, and Guangbin Zhang,
%          "Phase permutation entropy: A complexity measure for nonlinear time
%          series incorporating phase information."
%          Physica A: Statistical Mechanics and its Applications
%          568 (2021): 125686.
% 

narginchk(1,13)
p = inputParser;
Data = squeeze(Data);
Dn = size(Data,2);
N = size(Data,1);
Chk = @(x) isnumeric(x) && isvector(x) && (length(x)==Dn) && (min(x)>0) && all(mod(x,1)==0);
Chk2 = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
Chk3 = ["ampaware","weighted","edge","phase","modified"]; % 'uniquant','finegrain'};

addRequired(p,'Data',@(x) isnumeric(x) && ismatrix(x) && (Dn>1) && N>10);
addParameter(p,'m',2*ones(Dn,1),Chk);
addParameter(p,'tau',ones(Dn,1),Chk);
addParameter(p,'Typex','none',@(x) (ischar(x) || isstring(x)) && any(strcmpi(string(x),Chk3)));
addParameter(p,'tpx',[],@(x) (isnumeric(x) && isscalar(x) && (x > 0)) || isempty(x));
addParameter(p,'Logx',2,Chk2);
addParameter(p,'Norm',false,@(x) islogical(x))
parse(p,Data,varargin{:})
m = p.Results.m(:); tau = p.Results.tau(:);
Logx = p.Results.Logx; Norm = p.Results.Norm;
Typex = p.Results.Typex; tpx = p.Results.tpx;

if Logx == 0; Logx = exp(1); end
if strcmpi(Typex, "phase")
    Data = angle(hilbert(Data));
    if tpx == 1; Data = unwrap(Data);  end
end

Nx = N - max((m-1).*tau);
Sx = zeros(Nx,sum(m));
q = 1;
for k = 1:Dn
    for p=1:m(k)
        Sx(:,q) = Data(1+(p-1)*tau(k):Nx+(p-1)*tau(k),  k);
        q = q+ 1;
    end
end

% [TempUQ, Temp] = sort(Sx, 2,'ascend');
[~, Temp] = sort(Sx, 2,'ascend');
Px = unique(Temp, 'rows');
Counter = zeros(size(Px,1),1);
switch lower(Typex)
    case 'ampaware'
        if isempty(tpx)
            tpx = 0.5;
        elseif tpx<0 || tpx>1
            error('The A parameter (tpx) must be in the range [0 1]')
        end
        AA = sum(abs(Sx),2);
        AB = sum(abs(diff(Sx,[],2)),2);
        Ax = (tpx*AA/sum(m)) + ((1-tpx)*AB/(sum(m)-1));
        for n = 1:size(Px,1)
            Counter(n) = sum(Ax(~any(Temp-Px(n,:),2)));
        end
        Counter(Counter==0)=[];
        Ppi = Counter/sum(Ax);
        clear AA AB Ax

    case 'weighted'
        Wj = var(Sx,1,2);
        for n = 1:size(Px,1)
            Counter(n) = sum(Wj(~any(Temp - Px(n,:),2)));
        end
        Counter(Counter==0)=[];
        Ppi = Counter/sum(Wj);
        clear Wj

    case 'modified'
        Tx = (diff(sort(Sx,2),[],2)==0);
        for km = 1:sum(m)-1
            Temp(Tx(:,km),km+1) = Temp(Tx(:,km),km);
        end
        Px = unique(Temp,'rows');
        Counter = zeros(1,size(Px,1));
        for n = 1:size(Px,1)
            Counter(n) = sum(~any(Temp - Px(n,:),2));
        end
        Counter(Counter==0)=[];
        Ppi = Counter/sum(Counter);
        clear Tx km

    case 'edge'
        if isempty(tpx);    tpx = 1;
        elseif tpx <=0;  error('r sensitivity parameter (tpx) must be > 0 for edge permutation entropy')
        end

        for n = 1:size(Px,1)
            Tx = diff(Sx(~any(Temp - Px(n,:),2),:),[],2);
            Counter(n) = sum(mean(hypot(Tx,1),2).^tpx);
        end
        Counter(Counter==0)=[];
        Ppi = Counter/sum(Counter);
        clear Tx

    otherwise
        for n = 1:size(Px,1)
            Counter(n) = sum(~any(Temp - Px(n,:),2));
        end
        Counter(Counter==0)=[];
        Ppi = Counter/sum(Counter);
end

if round(sum(Ppi),5)~=1
    warning('Potential error with probability calculation')
end

MPerm = -sum(Ppi.*(log(Ppi)/log(Logx)));

if Norm;   MPnorm = MPerm/(log(factorial(sum(m)))/log(Logx));
else    MPnorm = MPerm/(sum(m)-1);
end

clear Temp Counter Nx Px Ppi
end

%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub