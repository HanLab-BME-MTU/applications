function [map] = getScatterQuantification(xData,yData,xBins,yBins)
%% 
% Input: xData,yData - data (scatter). xBins, yBins - bins to create heat map (scatter quantization)
% Output: heat map
%
% Assaf Zaritsky, June 2014
N = length(xData);
sizeX = length(xBins) - 1;
sizeY = length(yBins) - 1;
map = zeros(sizeY,sizeX);
for x = 1 : sizeX
    for y = 1 : sizeY
        count = xData >= xBins(x) & xData < xBins(x+1) & yData >= yBins(y) & yData < yBins(y+1);
        count = sum(count) ./ double(N);
        map(sizeY - y + 1,x) = count;
    end
end
end