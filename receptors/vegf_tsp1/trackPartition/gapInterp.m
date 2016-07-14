function [xOut,yOut] = gapInterp(x,y)
%GAPINTERP Interpolates particle position during gaps in a track (NaN's)
%
%   [xOut,yOut] = gapInterp(x,y)
%
%   Inputs: 
%       x,y:        vectors containing track coordinates; must be the same
%                   size
%
%   Outputs:
%       xOut,yOut:  vectors with linearly interpolated coordinates during
%                   track gaps
%
%Kevin Nguyen, July 2016
assert(sum(isnan(x) == isnan(y)) == numel(x),'x coordinates and y coordinates do not have the same gaps')
gaps = find(isnan(x));
xOut = x;
yOut = y;
xOut(gaps) = interp1(1:numel(x),x,gaps);
yOut(gaps) = interp1(1:numel(x),y,gaps);
end