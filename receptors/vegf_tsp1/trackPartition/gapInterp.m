function [xOut,yOut] = gapInterp(x,y)
assert(sum(isnan(x) == isnan(y)) == numel(x),'x coordinates and y coordinates do not have the same gaps')
gaps = find(isnan(x));
xOut = x;
yOut = y;
xOut(gaps) = interp1(1:numel(x),x,gaps);
yOut(gaps) = interp1(1:numel(x),y,gaps);
end