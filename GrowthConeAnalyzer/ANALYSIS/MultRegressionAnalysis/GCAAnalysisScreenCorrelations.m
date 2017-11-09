function [output,dataSetFinal,forDataSetArray,rAll,pAll] = GCAAnalysisScreenCorrelations(outgrowthDeltas,toPlot,varargin)
%  currently only assumes 1 group in toPlot but will loop through all the
%  params


%% Check input 
ip = inputParser;
ip.CaseSensitive = false;
% Check input
ip.addRequired('toPlot');
% Feature Selection Options
ip.addRequired('outgrowthDeltas');

% 
% Function for movie stat
ip.addParameter('perNeuriteStatistic','nanmean'); % default is to take the
% median value of the measurement calculated over the entire movie.
ip.addParameter('matrixPlot',true);
ip.addParameter('type','Pearson'); 
ip.parse(outgrowthDeltas,toPlot,varargin{:});

%%
perNeuriteStat = str2func(ip.Results.perNeuriteStatistic); 

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
    valuesPerGroup{iGroup,1} = perNeuriteStat(dataMat,1); % take the median over the rows- will get one value per movie
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
    
    
   % [r,p] = corrcoef(output(iParam).values);
    [r,p] = corr(output(iParam).values,'type',ip.Results.type); 
    
    output(iParam).r = r;
    output(iParam).p = p;
    output(iParam).v = valuesPerGroup; 
    
    
end
 


dataSetFinal = num2cell(forDataSetArray);
[rAll,pAll] = corr(forDataSetArray,'type',ip.Results.type); 


ObsNames = vertcat(toPlot.info.projList{:})  ; 
ObsNames = ObsNames(:,2); % get IDS
dataSetFinal = [ObsNames dataSetFinal]; 
params = ['NeuriteID' ; varName ;'Net Elongation Velocity'];
dataSetFinal = [params' ;dataSetFinal];

 %[rAll,pValuesAll] = corrplotMine(forDataSetArray,'testR','on','varNames',[varName ; 'Net Elongation Velocity']); 

dataSetFinal = cell2dataset(dataSetFinal,'ReadObsNames',true);

% save(['dataSet_' ip.Results.perNeuriteStatistic '.mat'],'dataSetFinal'); 


