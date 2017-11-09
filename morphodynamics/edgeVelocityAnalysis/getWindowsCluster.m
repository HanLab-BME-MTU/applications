function cellData = getWindowsCluster(movieObj,varargin)
%This function cluster windows based on edge velocity characteristics
%
%USAGE:
%       cellData = getWindowsCluster(movieObj,varargin)
%
%Input:
%       movieObj      - movie list or movie data object
%       includeWin    - cell array with indexes of the windows to be included for the analysis. Each element of the array corresponds to a cell from the movie list
%       outLevel      - scalar for outlier detection(see detectOutlier)
%       trendType     - scalar - algorihtm used to detrend the data(See getTimeSeriesTrend)
%       minLength     - scalar - mininum time serie length
%       maxLag        - scalar - maximum cross/auto-correlation lag
%       edgeState     - vector - 1 for protrusion;2- for retraction. Default [1 2]
%
%       stateFeature  - vector - 1 - persistence time
%                                2 - average velocity  
%                                3 - maximum velocity      
%                                4 - median velocity    
%                                All these features come from getEdgeMotionPersistence.m            
%                       
%       featureStats  - statistical operation applied to the stateFeature time series(see findTimeSeriesCluster
%                       Default - [2 3]. [mean variance]
%
%       clusterSet   - cellArray - input options for the clustering algorithm
%
%Output:
%       cellData - structure array with all parameters for each cell
%                          cellData(iCell).cluster - cellArray where each element is a vector with indexes for each cluster
%
%Marco Vilela, 2013

ip = inputParser;
ip.addRequired('movieObj',@(x) isa(x,'MovieList') || isa(x,'MovieData'));
ip.addParamValue('nBoot',1e3,@isscalar);
ip.addParamValue('alpha',.05,@isscalar);
ip.addParamValue('cluster',false,@isscalar);
ip.addParamValue('nCluster',2,@isscalar);


if isa(movieObj,'MovieData')
    
    ML = movieData2movieList(movieObj);
    
else
    
    ML = movieObj;
    
end

nCell = numel(ML.movies_);

%% Time Series Pre-Processing operations
ip.addParamValue('includeWin', cell(1,nCell),@iscell);
ip.addParamValue('winInterval',num2cell(cell(1,nCell)),@iscell);
ip.addParamValue('outLevel',  zeros(1,nCell),@isvector);
ip.addParamValue('trendType',-ones(1,nCell),@isvector);
ip.addParamValue('minLength', 30*ones(1,nCell),@isvector);
ip.addParamValue('scale',false,@islogical);
ip.addParamValue('edgeState',1:2,@isvector);
ip.addParamValue('stateFeature',1:4,@isvector);
ip.addParamValue('featureStats',@isvector );
ip.addParamValue('clusterMethod',1,@isscalar);
ip.addParamValue('clusterSet', {2,'Distance','sqEuclidean','Replicates',10}, @iscell);                 
ip.addParamValue('interval',num2cell(cell(1,nCell)),@iscell);

ip.parse(movieObj,varargin{:});
scale       = ip.Results.scale;
includeWin  = ip.Results.includeWin;
winInterval = ip.Results.winInterval;
interval    = ip.Results.interval;
outLevel    = ip.Results.outLevel;
minLen      = ip.Results.minLength;
trend       = ip.Results.trendType;
edgeState   = ip.Results.edgeState;
edgeFeat    = ip.Results.stateFeature;
fVector     = ip.Results.featureStats;
clusterM    = ip.Results.clusterMethod;
clusterSet  = ip.Results.clusterSet;

%%
edgeInputParam = {'interval',interval,'winInterval',winInterval,'outLevel',outLevel,'minLength',minLen,'trendType',trend,'includeWin',includeWin,'gapSize',ones(1,nCell),'scale',scale};
cellData       = edgeVelocityQuantification(ML,edgeInputParam{:});

% State Features

feature{1}   = 'persTime';
feature{2}   = 'Veloc';
feature{3}   = 'maxVeloc';
feature{4}   = 'mednVeloc';

%Edge States
measures{1}  = 'protrusion';
measures{2}  = 'retraction';

aux1         = arrayfun(@(x) repmat(x,1,numel(edgeFeat)),measures(edgeState),'Unif',0);
aux2         = repmat(feature(edgeFeat),1,numel(edgeState));
aux3         = cellfun(@(z1,z2) cellfun(@(y) arrayfun(@(x) x.(z1),y.([z2 'Analysis']).windows,'Unif',0),cellData,'Unif',0),aux2,cat(2,aux1{:}),'Unif',0);
featureSpace = cellfun(@(x) cat(2,x{:}),aux3,'Unif',0);

statsVector(1:numel(edgeFeat)*numel(edgeState)) = {fVector};

%Getting clusters
[out,dataPoints] = findTimeSeriesCluster(featureSpace,statsVector,clusterM,'clusterSet',clusterSet);
cellData         = getCellIndex(cellData,out,dataPoints,edgeFeat,edgeState,feature,measures);

%recounting the window's index
for iCell = 1:nCell

    cellData{iCell}.clusterWin = arrayfun(@(x) cellData{iCell}.data.includedWin{1}(cellData{iCell}.clusterWin{x}),1:numel(out),'Unif',0);
         
end

%Save results
savingMovieResultsPerCell(ML,cellData,'clusterAnalysis','cluster')
savingMovieDataSetResults(ML,cellData,'clusterAnalysis','cluster')

end

function cellData = getCellIndex(cellData,out,dataPoints,edgeFeat,edgeState,feature,measures)

totalWin = cell2mat( cellfun(@(x) numel(x.data.procExcEdgeMotion),cellData,'Unif',0) );

testM  = cellfun(@(x) sum(x'*(1./cumsum(totalWin)) > 1,2)+1, out,'Unif',0 );
fixIdx = [0 totalWin(1:end-1)];

for iCell = 1:numel(cellData)
    
    cellData{iCell}.clusterWin = cellfun(@(x,y) x(y == iCell)-sum(fixIdx(1:iCell)),out,testM,'Unif',0);
    
    for iM = 1:numel(edgeState)
        
        for iF = 1:numel(edgeFeat)
            
            cellData{iCell}.clusterPoint.(measures{iM}).(feature{iF}) = cellfun(@(x) dataPoints{iM+iF-1}(x,:),cellData{iCell}.clusterWin,'Unif',0);
    
        end
        
    end
    
end

end