function [robustSkewness] = skewnessRobustHinkley( data, varargin )
% Hinkley/Bowley's Estimate of skewness robust outliers
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
    p.addOptional('alpha', 0.25, @(x) (isscalar(x) && x > 0 && x < 0.5) );
    p.parse( data, varargin{:} );
    
    alpha = p.Results.alpha;
    
    Q1 = quantile( data(:), alpha );
    Q2 = median( data(:) );
    Q3 = quantile( data(:) , 1-alpha );
    
    robustSkewness = (Q3 + Q1 - 2 * Q2) / (Q3 - Q1);    
    
end