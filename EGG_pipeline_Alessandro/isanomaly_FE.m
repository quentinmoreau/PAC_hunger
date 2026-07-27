% +------------------------------------------------------+
% |    Outliers Detection with MATLAB Implementation     | 
% |                                                      |
% | Author: Ph.D. Eng. Hristo Zhivomirov        05/27/23 | 
% +------------------------------------------------------+
%
% function: outliers = isanomaly(x, method, context, gamma)
%
% Input:
% x - a vector with the input data sequence
% method - a method for outliers detection. Choose among:
%          'Wright' - Wright (Laiyite) criterion;
%          'modified Z' - Iglewicz-Hoaglin modified Z-score;
%          'MAD' - Huber-Miller MAD rule;
%          'liberal Tukey' - Tukey criterion using 1.5IQR fence;
%          'conservarive Tukey' - Tukey criterion using 3IQR fence;
%          'Romanowski' - Romanowski criterion (at confidence level gamma);
%          'Chebyshev' - Chebyshev criterion (at confidence level gamma).
% context - a context of the outliers detection. Choose among:
%           Inf - the estimation is performed globally using all data. This
%           setting is helpful if the data are uncorrelated or stationary.
%           finite positive integer number - the estimation is performed in
%           the vicinity of a predefined neighborhood (given as samples).
%           The data sample under estimation is situated in the middle of
%           this span. A context different from Inf is proper if the data
%           are corrleted and nonstationary.
% gamma - confidence level (e.g., 0.95, 0.99) for the hypothesis that the
%         detected suspicious values are really outliers. The default value
%         is 0.95. This input parameter is used only when the method for
%         outliers detection is set to 'Romanowski' or 'Chebyshev'
%         otherwise it can be omitted.
%
% Output:
% outliers - logical vector whose elements are true when an outlier is 
%            detected in the corresponding element of x.
%
% User guide:
% 1) Examine the modality of the data via histogram plot (not mandatory).
% 2) Choose an appropriate method for outliers detection from the list below...
%
% Method:                               Assumption:     Set size:   Masking:
% Wright (Laiyite) criterion            Normality,      Large set,  Sensitive       
% Iglewicz-Hoaglin modified Z-score     Normality,      Large set,  Not sensitive   
% Huber-Miller MAD rule                 Unimodality,    Arbitrary,  Not sensitive  
% Tukey criterion                       Unimodality,    Large set,  Not sensitive  
% Romanowski criterion                  Normality,      Arbitrary,  Sensitive  
% Chebyshev criterion                   Multimodality,  Large set,  Sensitive 
%
% Other methods exist too, but they are not so powerful, or their usage is
% not straightforward (e.g., Grubbs's test, Dixon's test, Chauvenet's
% criterion, etc.). Compared with the built-in Matlab function isoutlier,
% the following correspondences exist - 'Wright' == 'mean'; 'modified Z' ==
% 'median'; 'liberal Tukey' == 'quartiles'. Moreover, the 'grubs' and
% 'gesd' methods are not implemented in the present function, while 'MAD',
% 'Romanowski', and 'Chebyshev' methods are not implemented in the built-in
% function.
% 
% 3) Detected (potential) outlier(s) may be due to erroneous data or may
% indicate that the data are correct but highly unusual. A domain expert
% must perform further investigation of the cause for outlier(s) occurance.
%
% References:
% [1] S. Seo. A Review and Comparison of Methods for Detecting Outliers 
% in Univariate Data Sets (Master's Thesis). Pittsburgh, University of
% Pittsburgh, 2006. (Unpublished)
% [2] B. Iglewicz, D. Hoaglin. ASQC Basic References in Quality Control
% Vol. 16: How To Detect And Handle Outliers. Milwaukee, ASQC Quality
% Press, 1993.
% [3] C. Leys, C. Ley, O. Klein, P. Bernard, L. Licata. Detecting outliers:
% Do not use standard deviation around the mean, use absolute deviation
% around the median. Journal of Experimental Social Psychology, Vol. 49,
% No. 4, pp. 764-766, 2013.
% [4] B. Amidan, T. Ferryman, S. Cooley. "Data outlier detection using the
% Chebyshev theorem". IEEE Aerospace Conference, pp. 3814-3819, 2005.
function outliers = isanomaly_FE(x, method, context, gamma)
% input validation
if nargin == 3, gamma = 0.95; end
validateattributes(x, {'single', 'double'}, ...
                      {'vector', 'real', 'nonempty', 'finite'}, ...
                      '', 'x', 1)
validateattributes(method, {'char'}, ...
                           {'scalartext'}, ...
                           '', 'method', 2)
validateattributes(context, {'single', 'double'}, ...
                            {'scalar', 'nonnan', 'nonempty', 'positive'}, ...
                            '', 'context', 3)
validateattributes(gamma, {'single', 'double'}, ...
                          {'scalar', 'real', 'nonnan', 'nonempty', 'positive', ...
                           '>', 0, '<', 1}, ...
                           '', 'gamma', 4)
% data shaping 
if isinf(context)
    % form a contextual matrix and define the actual data
    X = x(:);
    xhat = x(:)';
    correction = 0;
elseif context < length(x)
    % form a contextual matrix and define the actual data
    [X, ~] = buffer(x, round(context), round(context-1), 'nodelay');
    xhat = X(ceil(size(X, 1)/2), :);
    correction = ceil(size(X, 1)/2) - 1;
else
    % throw an error message
    error('The context must be Inf or a number smaller than the data length!')
end
% outliers detection
switch lower(method)
    case 'wright'
        %% detect outliers via Wright (Laiyite) criterion (using the Z-score)
        % Note: be aware that the built-in zscore fucntion does not omit
        % NaNs so do not use it explicitly here!
        Zx = (xhat - mean(X, 'omitnan'))./std(X, 'omitnan');
        outlrs_ind = find(abs(Zx) > 3);   
        
    case 'modified z'
        %% detect outliers via the Iglewicz-Hoaglin modified Z-score (using the MAD)
        % Note: MAD stends for "median of the absolute deviations about the
        % median". The built-in mad function can be used explicitly here.
        % Be aware that, if more than 50% of data have identical values,
        % the MAD will equal zero and so the modified Z-score will be Inf
        % yealding a false possitive outlier alarm. For this reason all Inf
        % values of the Z-score are set to NaN.
        Mx = median(X, 'omitnan');
        MADx = median(abs(X - Mx), 'omitnan');
        Zxm = 0.6745*(xhat - Mx)./MADx;
        Zxm(isinf(Zxm)) = NaN;
        outlrs_ind = find(abs(Zxm) > 3.5);
    case 'mad'
        %% detect outliers via Huber-Miller MAD rule (using the MAD)
        % Note: MAD stends for "median of the absolute deviations about the
        % median". The built-in mad function can be used explicitly here.
        % The quantile function treats NaNs as missing values and removes
        % them.
        Mx = median(X, 'omitnan');
        Q75 = quantile(X, 0.75);
        MADex = 1./Q75.*median(abs(X - Mx), 'omitnan');
        outlrs_ind = find(abs(xhat) > 2.5*(Mx + MADex));
    case 'liberal tukey'
        %% detect outliers via Tukey criterion (using the interquartile range)
        % Note: The quantile function treats NaNs as missing values and
        % removes them. The term 'liberal' means that one uses the 1.5*IQR
        % limit as a outlier threshold.
        Qx = quantile(X, [0.25, 0.5, 0.75], 1);
        IQRx = Qx(3, :) - Qx(1, :);
        outlrs_ind = find(xhat < Qx(1, :)-1.5*IQRx | xhat > Qx(3, :)+1.5*IQRx);
        
    case 'conservative tukey'
        %% detect outliers via Tukey criterion (using the interquartile range)
        % Note: The quantile function treats NaNs as missing values and
        % removes them. The term 'conservative' means that one uses the
        % 3*IQR limit as a outlier threshold.
        Qx = quantile(X, [0.25, 0.5, 0.75], 1);
        IQRx = Qx(3, :) - Qx(1, :);
        outlrs_ind = find(xhat < Qx(1, :)-3*IQRx | xhat > Qx(3, :)+3*IQRx);
    case 'romanowski'
        %% detect outliers via Romanowski criterion (using the t-statistics)
        Zx = (xhat - mean(X, 'omitnan'))./std(X, 'omitnan');
        St = tinv(1-(1-gamma)/2, size(X, 1)-1);
        outlrs_ind = find(abs(Zx) > St*sqrt(1+1/size(X, 1)));
    case 'chebyshev'
        %% detect outliers via Chebyshev criterion (using the Chebyshev Theorem)
        mx = mean(X, 'omitnan');
        sx = std(X, 'omitnan');
        k = 1/sqrt(1-gamma);
        outlrs_ind = find(abs(xhat) > mx + k*sx);
        
    otherwise
        % throw an error message
        error('Choose a valid method for outliers detection (see the help)!')
end
% form the output of the function
outlrs_ind = outlrs_ind + correction;
outliers = false(size(x));
outliers(outlrs_ind) = true;
end