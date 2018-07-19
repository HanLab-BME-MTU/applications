function xCorrMap = nanXcorrMaps(xmap, ymap, varargin)
% nanXcorrMaps Compute cross correlations between two activity maps. The cross
% correlations at lag h are Corr(xmap_{t+h}, ymap_t).
% This function can handle many NaN's due to using nanXcorrelatioon.m
% function.
%
% Jungsik Noh, 2016/10/xx

%%
tmax0 = size(xmap, 2);

ip = inputParser;
ip.addRequired('xmap',@(x)(ismatrix(x) && min(size(x)) > 1));
ip.addRequired('ymap',@(x)(ismatrix(x) && min(size(x)) > 1));

ip.addParameter('lagMax', round(tmax0/4), @isnumeric);
ip.parse(xmap, ymap, varargin{:});
p = ip.Results;

%%
map1 = xmap' - repmat(mean(xmap', 1, 'omitnan'), size(xmap', 1), 1);
map2 = ymap' - repmat(mean(ymap', 1, 'omitnan'), size(ymap', 1), 1);

%lagMax = 50;
lagMax = p.lagMax;
numWindows = size(map1, 2);

xCorrMap = zeros(numWindows, 2*lagMax+1);

for w=1:numWindows
    %xCorrMap(w, :) = xcorr(map1(:, w), map2(:, w), lagMax, 'coeff');    
    xCorrMap(w, :) = nanXcorrelation(map1(:, w), map2(:, w), lagMax);
end

 

