function rms = residuals(obj,modelIdx)
if nargin > 1 % Only for the specified model
    errors = obj.data.error(obj.data.clusters{modelIdx},:);
    normRes = obj.data.modelRes{modelIdx}./errors;
else % For all models
    normRes = cellfun(@(a,b) a./obj.data.error(b,:),obj.data.modelRes,obj.data.clusters,'UniformOutput',false);
    normRes = vertcat(normRes{:});
end
normDist = sqrt(sum(normRes.^2,2));
% Display histogram
rms = sqrt(sum(normDist.^2)/numel(normDist));
nBins = round(sqrt(numel(normDist)));
hist(normDist,nBins);
xlabel('Normalized distances');
ylabel('Number of data points');
disp('Analyzer: Normalized residuals histogram displayed!');
end