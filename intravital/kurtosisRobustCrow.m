function [robustKurtosis] = kurtosisRobustCrow( data, varargin )
% A Robust estimate of kurtosis 
%
% References:
% 
%     Tae-Hwan Kim, Halbert White, "On more robust estimation of skewness and kurtosis", 
%     Finance Research Letters, Volume 1, Issue 1, March 2004, Pages 56-73
%
% Author: Deepak Roy Chittajallu (2013)
%
%
    p = inputParser;
    p.addRequired('data', @(x) (isnumeric(x)) );
    p.addParamValue('alpha', 0.025, @(x) (isscalar(x) && x > 0 && x < 1) );
    p.addParamValue('beta', 0.25, @(x) (isscalar(x) && x > 0 && x < 1) );
    p.parse( data, varargin{:} );
    
    alpha = p.Results.alpha;
    beta = p.Results.beta;
    
    qFunc = @(q) ( quantile(data(:), q) );
    
    robustKurtosis = ((qFunc(1-alpha) + qfunc(alpha)) / (qfunc(1 - beta) - qfunc(beta))) - 2.91;

end