%% Input: data, labels, parameters
% output: histograms + diff hostograms
function [out] = histogramQuantification(data,labels,uniqueLabelsStr,uniqueLabelsValues, params)
if nargin < 4
    params.percentile = 2;
    params.nBins = 20;
end

out.uniqueLabelsStr = uniqueLabelsStr;
out.uniqueLabelsValues = uniqueLabelsValues;

nLabels = length(out.uniqueLabelsStr);

lowPrctl = prctile(data,params.percentile);
highPrctl = prctile(data,100-params.percentile);

out.bins = lowPrctl : (highPrctl - lowPrctl)/(params.nBins-1) : highPrctl;

out.hists = cell(1,nLabels);
out.Ns = nan(1,nLabels);
for curInd = 1 : nLabels
        data4Label = data(labels==curInd);        
        [histogram] = getHistogram(data4Label,out.bins);
        out.hists{curInd} = histogram./sum(histogram(:));
        out.Ns(curInd) = length(data4Label);        
end

out.diffHists = {};
out.diffMatrix = nan(nLabels,nLabels);
for i = 1 : nLabels
    for j = i+1 : nLabels
        hist1 = out.hists{i};
        hist2 = out.hists{j};
        out.diffHists{i,j} = abs(hist1-hist2);%./max(map1,map2);        
        out.diffMatrix(i,j) = sum(out.diffHists{i,j});
    end
end

end


function [histogram] = getHistogram(data,bins)
N = length(data);
M = length(bins) - 1;

histogram = zeros(1,M);
for m = 1 : M
    count = data >= bins(m) & data < bins(m+1);
    count = sum(count) ./ double(N);
    histogram(m) = count;
end
end