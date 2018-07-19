function [pearsonCoff] = corrPearson( x, y )
% Computes the pearson correlation coefficient
%
% Author: Deepak Roy Chittajallu
%
    covMat = cov( x(:), y(:) );    
    pearsonCoff = covMat(2) / (std(x(:)) * std(y(:)));

end