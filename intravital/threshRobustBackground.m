function [threshVal] = threshRobustBackground(im, chopPercent)
% The threshold is calculated by trimming the top and bottom 5% of
% pixels (orered by intensity) off the image, then calculating the mean 
% and standard deviation of the remaining image. The threshold is then 
% set at 2 (empirical value) standard deviations above the mean.
%
% Required Input Arguments:
% 
%              im: ND image
% 
% Optional Input Arguments:
% 
%     chopPercent: percentage of pixels to chop at the top 
%                  and bottom
% 
% Output:
% 
%       threshVal: intensity value of the threshold
% 
% Author: Deepak Roy Chittajallu
%

   if ~exist('chopPercent', 'var')
        chopPercent = 5; % percentage of pixels to chop at the top and bottom
   end
   
   grayvals = sort( im(:) );
   i1 = round((chopPercent/100) * numel(im));
   i2 = numel(im) - i1;
   
   mu = mean( grayvals(i1:i2) );
   sigma = std( grayvals(i1:i2) );
   
   threshVal = mu + 2 * sigma;

end