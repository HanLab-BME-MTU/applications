function [output,dataSetFinal] = GCAAnalysisScreenCorrelations(outgrowthDeltas,toPlot)
%  currently only assumes 1 group in toPlot but will loop through all the
%  params

% Get all the parameters to correlate to global function
paramsToCorrelate = fieldnames(toPlot);
% take out the info field
paramsToCorrelate(cellfun(@(x) strcmpi(x,'info'),paramsToCorrelate)) = [];
outgrowthDeltas = outgrowthDeltas./10; % convert to um/min 
varName = paramsToCorrelate; 

% initiate a predictor/response Mat r equals number of observations (ie
% neurites in this case) and c is the number of variables - the last column
% will be the response variable which is in this case the neurite
% outgrowth.
forDataSetArray = nan(length(outgrowthDeltas),numel(paramsToCorrelate)+1);
forDataSetArray(:,end) = outgrowthDeltas; 
%%
% *BOLD TEXT* tArray(:,end) = outgrowthDeltas;

nGroups = numel(toPlot.info.names); 
for iParam = 1:numel(paramsToCorrelate)
     output(iParam).name = varName{iParam};
for iGroup = 1:nGroups

    
   
    % extract the current dataMat : it is format of rows = individual observations for
    % whole movie (if parameter has a calculation per frame (for example:
    % filopodia density) the rows will be equal to one.
    dataMat = toPlot.(paramsToCorrelate{iParam}).dataMat{iGroup};
    dataMat = real(dataMat);
   
    % to simplify take the median value of the parameter per movie.
    valuesPerGroup{iGroup,1} = nanmedian(dataMat,1); % take the median over the rows- will get one value per movie
    % (column) in list.
end % for iGroup 
    valuesC = horzcat(valuesPerGroup{:}); 
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
    output(iParam).v = valuesPerGroup; 
    
    
end
 
dataSetFinal = num2cell(forDataSetArray);
ObsNames = vertcat(toPlot.info.projList{:})  ; 
ObsNames = ObsNames(:,2); % get IDS
dataSetFinal = [ObsNames dataSetFinal]; 
params = ['NeuriteID' ; varName ;'Net Velocity'];
dataSetFinal = [params' ;dataSetFinal];


dataSetFinal = cell2dataset(dataSetFinal,'ReadObsNames',true);

save('dataSet.mat','dataSetFinal'); 


