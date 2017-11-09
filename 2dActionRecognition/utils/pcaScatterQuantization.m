%% Input: data for pca, labels, parameters
% output: scatter-quantization maps
function [out] = pcaScatterQuantization(data,labels,uniqueLabelsStr,uniqueLabelsValues,params)
if nargin < 4
    params.percentile = 2;
    params.nBins = 20;
end

[out.coeff,out.score,out.latent] = pca(data');
out.accVariance = cumsum(out.latent)./sum(out.latent);

out.uniqueLabelsStr = uniqueLabelsStr;
out.uniqueLabelsValues = uniqueLabelsValues;

nLabels = length(out.uniqueLabelsStr);

pctl1Low = prctile(out.score(:,1),params.percentile);
pctl1High = prctile(out.score(:,1),100-params.percentile);
pctl2Low = prctile(out.score(:,2),params.percentile);
pctl2High = prctile(out.score(:,2),100-params.percentile);

out.xBins = pctl1Low : (pctl1High - pctl1Low)/(params.nBins-1) : pctl1High;
out.yBins = pctl2Low : (pctl2High - pctl2Low)/(params.nBins-1) : pctl2High;


out.maps = cell(1,nLabels);
out.Ns = nan(1,nLabels);
for curInd = 1 : nLabels
        xData = out.score(labels==uniqueLabelsValues(curInd),1);
        yData = out.score(labels==uniqueLabelsValues(curInd),2);        
        [map] = getScatterQuantization(xData,yData,out.xBins,out.yBins);
        out.maps{curInd} = map./sum(map(:));
        out.Ns(curInd) = length(yData);        
end

out.diffMaps = {};
out.diffMatrix = nan(nLabels,nLabels);
for i = 1 : nLabels
    for j = i+1 : nLabels
        map1 = out.maps{i};
        map2 = out.maps{j};
        out.diffMaps{i,j} = abs(map1-map2);%./max(map1,map2);  
        tmpDiff = out.diffMaps{i,j};
        out.diffMatrix(i,j) = sum(tmpDiff(:));
    end
end

end


function [map] = getScatterQuantization(xData,yData,xBins,yBins)
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