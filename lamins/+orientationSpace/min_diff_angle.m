function [ maxima_d_min ] = min_diff_angle( maxima, period, dim )
%min_diff_angle Finds the minimum nearest angle-angle difference
%
% INPUT
% maxima - orientation angles as from OrientationSpace.getRidgeOrientationLocalMaxima;
% period - periodic boundary of orientations in maxima (default: pi)
% dim - dimension overwhich to take difference (default: 3);
%
% OUTPUT
% minimum angle-angle difference map
if(nargin < 2)
    period = pi;
end
if(nargin < 3)
    dim = 3;
end

maxima_s = sort(maxima,dim);
[X,Y] = meshgrid(1:size(maxima,2),1:size(maxima,1));
maxima_c = sum(~isnan(maxima),dim);
maxima_last = maxima_s(sub2ind(size(maxima),Y,X,maxima_c));
maxima_d = diff(cat(dim,maxima_last,maxima_s),1,dim);
maxima_d = min(abs(maxima_d),period-abs(maxima_d));
maxima_d_min = min(maxima_d,[],dim);

end

