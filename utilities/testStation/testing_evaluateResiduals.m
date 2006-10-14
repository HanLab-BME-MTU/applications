function testing_evaluateResiduals(dbData)

% plot residuals as function of source, target

% dbData has a field trackResults
if ~isfield(dbData,'trackResults') || isempty(dbData,'trackResults')
    error('wrong function/wrong input')
end
trackResults = dbData.trackResults;
clear dbData

sizeResults = size(trackResults);
nTags = sizeResults(2);

% plotData is nTimepoints*nSources-by-nTags+2, and has cols
% sourceT, targetT, res2 of tags 1:nTags
% plotDataInit is the same, except that the intial residuals are collected

plotData = zeros(prod(sizeResults([1,3])),nTags+2);
plotDataInit = plotData;

% loop through timepoints, tags, sources to fill plotData
dataCt = 0;

for t = 1:sizeResults(1)
    for iSource = 1:sizeResults(3)
        
        if ~isempty(trackResults(t,1,iSource).info)
        
        sourceT = trackResults(t,1,iSource).info(3);
        targetT = trackResults(t,1,iSource).info(4);
        
        % collect residuals
        residuals = cat(1,trackResults(t,:,iSource).sigma0);
        
        dataCt = dataCt+1;
        
        % write plotData
        plotData(dataCt,1) = sourceT;
        plotData(dataCt,2) = targetT;
        plotData(dataCt,3:end) = residuals(:,2)';
        
        % write plotDataInit
        plotDataInit(dataCt,1) = sourceT;
        plotDataInit(dataCt,2) = targetT;
        plotDataInit(dataCt,3:end) = residuals(:,1)';
        
        
        
        end
        
    end
end
            
plotData(dataCt+1:end,:) = [];
plotDataInit(dataCt+1:end,:) = [];


% ind=sub2ind([200,200],plotData(:,1),plotData(:,2));
% m=zeros(200);
% m(ind) = plotData(:,3);
% figure,surf(m)

            