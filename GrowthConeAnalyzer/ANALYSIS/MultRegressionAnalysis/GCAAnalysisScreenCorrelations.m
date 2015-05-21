function [output,dataSetFinal] = GCAAnalysisScreenCorrelations(outgrowthDeltas,toPlot)
%  currently only assumes 1 group in toPlot but will loop through all the
%  params

% Get all the parameters to correlate to global function
paramsToCorrelate = fieldnames(toPlot);
% take out the info field
paramsToCorrelate(cellfun(@(x) strcmpi(x,'info'),paramsToCorrelate)) = [];
outgrowthDeltas = outgrowthDeltas./10; % convert to um/min 
varName = paramsToCorrelate; 
% varName = cellfun(@(x) strrep(x,'filo','F'),varName,'uniformoutput',0); 
% varName = cellfun(@(x) strrep(x,'Length','L'),varName,'uniformoutput',0); 
% varName = cellfun(@(x) strrep(x,'Intensity','I'),varName,'uniformoutput',0); 
% varName = cellfun(@(x) strrep(x,'retractionAnalysis_perTime','RPers'),varName,'uniformoutput',0); 
% varName = cellfun(@(x) strrep(x,'protrusionAnalysis_perTime','PPers'),varName,'uniformoutput',0); 
% varName = cellfun(@(x) strrep(x,'protrusionAnalysis_mednVel','PVel'),varName,'uniformoutput',0); 
% varName = cellfun(@(x) strrep(x,'retractionAnalysis_mednVel','RVel'),varName,'uniformoutput',0); 

% initiate a predictor/response Mat r equals number of observations (ie
% neurites in this case) and c is the number of variables - the last column
% will be the response variable which is in this case the neurite
% outgrowth.
forDataSetArray = nan(length(outgrowthDeltas),numel(paramsToCorrelate)+1);
forDataSetArray(:,end) = outgrowthDeltas;


for iParam = 1:numel(paramsToCorrelate)
    
    output(iParam).name = varName{iParam};
    % extract the current dataMat : it is format of rows = individual observations for
    % whole movie (if parameter has a calculation per frame (for example:
    % filopodia density) the rows will be equal to one.
    dataMat = toPlot.(paramsToCorrelate{iParam}){1};
    dataMat = real(dataMat);
    if strcmpi(paramsToCorrelate{iParam},'retractionAnalysis_persTime') || strcmpi(paramsToCorrelate{iParam},'protrusionAnalysis_persTime')
        nMovies = size(dataMat,2);
        % get the 75th percentile of each movie for the persistence time
        % /retraction time % set anything lower to this value to zero.
        prctile75 = arrayfun(@(i) prctile(dataMat(:,i),75),1:nMovies);
        toRemove = arrayfun(@(i) dataMat(:,i)<prctile75(i),1:nMovies,'uniformoutput',0);
        toRemoveMat = horzcat(toRemove{:});
        dataMat(toRemoveMat)=NaN; % set to NaN;
        
    end
    
    % to simplify take the median value of the parameter per movie.
    valuesC = nanmedian(dataMat,1); % take the median over the rows- will get one value per movie
    % (column) in list.
    
    forDataSetArray(:,iParam) = valuesC';
    
    if sum(isnan(valuesC))>0
        problemProj = find(isnan(valuesC));
        
        display(['NaN found for ' num2str(problemProj) 'and' paramsToCorrelate{iParam}]);
        display('Removing Value : Check Project') ;
        
        outgrowthDeltasTrunc = outgrowthDeltas(~isnan(valuesC));
        valuesC = valuesC(~isnan(valuesC));
        output(iParam).values= [outgrowthDeltasTrunc valuesC'];
    else
        output(iParam).values = [outgrowthDeltas valuesC'];
    end
    
    
    [r,p] = corrcoef(output(iParam).values);
    
    
    output(iParam).r = r;
    output(iParam).p = p;
    
    
    
end

dataSetFinal = num2cell(forDataSetArray);
ObsNames = toPlot.info.projList{1}(:,2); % get IDS
dataSetFinal = [ObsNames dataSetFinal]; 
params = ['NeuriteID' ; varName ;'Net Velocity'];
dataSetFinal = [params' ;dataSetFinal];


dataSetFinal = cell2dataset(dataSetFinal,'ReadObsNames',true);

save('dataSet.mat','dataSetFinal'); 


% collect p
% pMat = cellfun(@(x) x(1,2),p);
% %toTestC = toTest;
% %paramHits = paramsToCorrelate(p<0.05);
% outputHit = [paramsToCorrelate(pMat<0.05)  toTest(pMat<0.05)'  r(pMat<0.05)'   p(pMat<0.05)' ];
% resultsCorrScreen.Hit = outputHit;
%
%
% outputNonHit= [paramsToCorrelate(pMat>0.05)  toTest(pMat>0.05)' r(pMat>0.05)'  p(pMat>0.05)'];
% resultsCorrScreen.nonHit = outputNonHit;