function [bandpassFreqFilter] = buildSimoncelliBandpassFreqFilter(levels,sz)

if(nargin < 2)
    sz = 1024;
end

if(isscalar(sz))
    sz = [sz sz];
end

[X, Y] = meshgrid( (0:sz(2)-1) - floor(sz(2)/2), (0:sz(1)-1) - floor(sz(1)/2) );
[theta, f] = cart2pol(X / floor(sz(2)/2) / 2,Y / floor(sz(1)/2) / 2);

lowpassFreqFilter = ones([sz levels+1]);
highpassFreqFilter = ones([sz levels+1]);

for level=1:levels
    [lowpassFreqFilter(:,:,level), highpassFreqFilter(:,:,level)] = ...
        pyramid.raisedCosLogFreqFilter(f*pi*2^(level-1));
end

bandpassFreqFilter = lowpassFreqFilter(:,:,1:end-1).*highpassFreqFilter(:,:,2:end);
bandpassFreqFilter = cat(3,highpassFreqFilter(:,:,1),bandpassFreqFilter,lowpassFreqFilter(:,:,end));


end