function maxGap(obj)
% Compute the max gap for each model
bigClusterIdx = 1:obj.data.nClusters;
maxGaps = zeros(numel(bigClusterIdx),1);
clusterCount = 0;
maxGapProb = zeros(numel(bigClusterIdx),1);

% Loop throught all the models
for c=bigClusterIdx
    % Determine the maximum gap of the curve
    clusterCount = clusterCount + 1;
    t = sort(obj.data.modelProj{c});
    startTs = t(2:end);
    endTs = t(1:end-1);
    fun = @(startTs,endTs) lengthBezier(obj.data.modelBezCP{c},startTs,endTs);
    gaps = arrayfun(fun,startTs,endTs);
    [maxGap,idx] = max(gaps);
    startT = startTs(idx);
    endT = endTs(idx);
    maxGaps(clusterCount) = maxGap;
    
    % Compute the probability of the gaps
    curveLength = lengthBezier(obj.data.modelBezCP{c});
    p = maxGap/curveLength;
    n = numel(obj.data.modelProj{c});
    k = 0;
    maxGapProb(clusterCount) = binopdf(k,n,p);
end

% Display the max gap histogram
subplot(2,1,1)
nBins = round(sqrt(numel(maxGaps)));
hist(maxGaps,2*nBins);
xlabel('Max gap length [nm]');
ylabel('Number of gaps');

% Display the max gap probability histogram
subplot(2,1,2)
hist(maxGapProb,2*nBins);
xlabel('Max gap probability');
ylabel('Number of gaps');
end