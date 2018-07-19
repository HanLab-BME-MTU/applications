function [medianSpacing,probDensFct] = spacings(obj,display)
nSpacings = cellfun(@(a) numel(a)-1,obj.data.modelProj);
spacings = zeros(sum(nSpacings),1);
cumNSpacings = cumsum(nSpacings);
idxStart = [1;cumNSpacings(1:end-1)+1];
idxEnd = cumNSpacings;
k = 0;
for i=1:obj.data.nClusters
    k = k + 1;
    t = obj.data.modelProj{i};
    t = sort(t);
    spacings_i = arrayfun(@(tStart,tEnd) lengthBezier(obj.data.modelBezCP{i},tStart,tEnd),t(2:end),t(1:end-1));
    spacings(idxStart(k):idxEnd(k)) = spacings_i;
end

% Display histogram
nBins = round(sqrt(numel(spacings)))*10;
[counts,xout] = hist(spacings,nBins);

if nargin > 1
    if strcmp(display,'on')
        bar(xout,counts);
    end
else
    bar(xout,counts);
end
medianSpacing = median(spacings);

% Compute the probability density function
probDensFct = [xout;counts/sum(counts)]';

end