function CMLags = xcorrMatToCMLag(xcorrMatCell)
% xcorrMatToCMLag Extract a feature of correlation maxima and their lags
% from the xcorrMat per movie.
%   OUTPUT:     r0, lag0 per layer
% 
% Jungsik Noh, 02/07/2017

% Updated 2017/03/31: 
%           Instead of mean(), it uses smoothingSpline.

numLayer = numel(xcorrMatCell);
lagMax = (size(xcorrMatCell{1}, 2) - 1)/2;

CMLags = cell(numLayer, 1);
for l = 1:numLayer
    xcmat = xcorrMatCell{l};
    %xcmeanVec = mean(xcmat, 1, 'omitnan');
    xcmeanVec = smoothingSplineCorMap(xcmat);
    [r0, i] = max(abs(xcmeanVec));
    l0 = i - lagMax - 1;
    r1 = xcmeanVec(i);
    CMLags{l} = [r1, l0];
end


end
