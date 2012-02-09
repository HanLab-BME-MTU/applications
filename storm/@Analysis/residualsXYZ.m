function [sigX,sigY,sigZ] = residualsXYZ(obj,modelIdx)
if nargin > 1 % Only for the specified model
    errors = obj.data.error(obj.data.clusters{modelIdx},:);
    normRes = reshape(obj.data.modelRes{modelIdx},size(errors,1),size(errors,2))./errors;
else % For all models
    nNormRes = cellfun(@numel,obj.data.clusters);
    normRes = zeros(sum(nNormRes),3);
    cumNNormRes = cumsum(nNormRes);
    idxStart = [1;cumNNormRes(1:end-1)+1];
    idxEnd = cumNNormRes;
    k = 0;
    for i=1:obj.data.nClusters
        k = k + 1;
        errors = obj.data.error(obj.data.clusters{i},:);
        normRes_i = reshape(obj.data.modelRes{i},size(errors,1),size(errors,2))./errors;
        normRes(idxStart(k):idxEnd(k),:) = normRes_i;
    end
end
% Display histogram
nBins = 2*round(sqrt(size(normRes,1)));

maxVal = max(max(abs(normRes)));
binSize = 2*maxVal/nBins;
binCenters = -maxVal + (0:nBins-1)*binSize + binSize/2;

hX = subplot(3,1,1);
hist(normRes(:,1),binCenters);
sigX = sqrt(1/numel(normRes(:,1))*sum(normRes(:,1).^2));
xlabel('Normalized residuals X');
ylabel('Number of data points');

hY = subplot(3,1,2);
hist(normRes(:,2),binCenters);
sigY = sqrt(1/numel(normRes(:,2))*sum(normRes(:,2).^2));
xlabel('Normalized residuals Y');
ylabel('Number of data points');

hZ = subplot(3,1,3);
hist(normRes(:,3),binCenters);
sigZ = sqrt(1/numel(normRes(:,3))*sum(normRes(:,3).^2));
xlabel('Normalized residuals Z');
ylabel('Number of data points');

% Set the new limits
newLim = repmat(max([abs(get(hX,'XLim')) abs(get(hY,'XLim')) abs(get(hZ,'XLim'))]),1,2); newLim(1) = -newLim(1);
set(hX,'XLim',newLim); set(hY,'XLim',newLim); set(hZ,'XLim',newLim);
newLim = repmat(max([abs(get(hX,'YLim')) abs(get(hY,'YLim')) abs(get(hZ,'YLim'))]),1,2); newLim(1) = 0;
set(hX,'YLim',newLim); set(hY,'YLim',newLim); set(hZ,'YLim',newLim);

disp('Analyzer: Normalized residuals histogram displayed!');
end