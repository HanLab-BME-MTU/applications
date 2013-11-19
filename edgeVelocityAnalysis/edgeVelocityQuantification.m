function [cellData,dataSet] = edgeVelocityQuantification(movieObj,varargin)
% This function quantifies:
%                           mean protrusion/retraction instantaneous values and confidence interval;
%                           mean protrusion/retraction persistence time and confidence interval;
%                           cluster velocities and persistence time
%                           for each window: persistence time;protrusion/retraction block;Upper and lower noise threshold
%
% All the above quantifications are done for each cell of the movieObj (movieList or movieData)
% If the input 'interval' is set, all the described quantifications are calculated for each time interval
%
%
%Usage: [cellData,dataSet] = edgeVelocityQuantification(movieObj,varargin)
%
% Input:
%       nBoot         - # of boostrap samples to be used (default value 1000)
%
%       alpha         - alpha used to generate the bootstrap confidence intervals
%                       (default value 0.05)
%
%       includeWin -
%
%       scale         - convert velocity from pixel/frame into nm/sec.Default:false
%
%       includeWin    - cell array. Each element is a vector with the window indexes to be used in the analysis. Enter [] to analize all windows.
%
%       outLevel      - scalar for outlier detection(see detectOutlier)
%
%       trendType     - scalar - algorihtm used to detrend the data(See getTimeSeriesTrend)
%
%       minLength     - scalar - mininum time serie length
%
%       gapSize       - scalar - length of the nan gap to be closed. Default:0
%
%
%       outputPath    - string. Name of the output folder. Default: 'EdgeVelocityQuantification'
%
%       fileName      - string. Name of the output file. Default: 'edgeVelocity'
%
% Output:
%       cellData - this is a long structure array cellData(for each cell) with the following fields:
%                .data.excludeWin - indexes of windows that were excluded from analysis . Ex: border windows
%                     .rawEdgeMotion - raw edge velocity time series. Protrusion map. This is redundant.
%                     .procEdgeMotion - pre-processed edge velocity time series. Ex: trend, mean and NaN removed.
%
%                .protrusionAnalysis(for each interval).meanValue.persTime
%                                                                .maxVeloc
%                                                                .minVeloc
%                                                                .meanVeloc
%                                                                .mednVeloc
%                                                                .meanValue
%                                                                .cluster.persTime
%                                                                        .maxVeloc
%                                                                        .minVeloc
%                                                                        .meanVeloc
%                                                                        .mednVeloc
%                                                                        .meanValue
%
%           The structure continues with the confidence interval for each of the measurement above
%
%                                                      CI.persTimeCI
%                                                        .maxVelocCI
%                                                        .minVelocCI
%                                                        .meanVelocCI
%                                                        .mednVelocCI
%                                                        .meanValueCI
%                                                        .cluster.persTime
%                                                                .maxVeloc
%                                                                .minVeloc
%                                                                .meanVeloc
%                                                                .mednVeloc
%                                                                .meanValue
%
%                                                       Analysis at the individual time series level
%
%                                                      .windows(for each window).limit
%                                                                               .PersTime
%                                                                               .BlockOut
%                                                                               .MaxVeloc
%                                                                               .MeanVeloc
%                                                                               .MinVeloc
%                                                                               .MednVeloc
%
%Marco Vilela, 2012


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
ip.addParamValue('includeWin',cell(1,nCell),@iscell);
ip.addParamValue('winInterval',num2cell(cell(1,nCell)),@iscell);
ip.addParamValue('outLevel',  zeros(1,nCell),@isvector);
ip.addParamValue('trendType',-ones(1,nCell),@isvector);
ip.addParamValue('minLength', 10*ones(1,nCell),@isvector);
ip.addParamValue('gapSize',   zeros(1,nCell),@isvector);
ip.addParamValue('scale',     false,@islogical);
ip.addParamValue('outputPath','EdgeVelocityQuantification',@isstr);
ip.addParamValue('fileName','EdgeMotion',@isstr);
ip.addParamValue('interval',num2cell(cell(1,nCell)),@iscell);
ip.addParamValue('lwPerc',2.5,@isscalar);
ip.addParamValue('upPerc',97.5,@isscalar);
ip.addParamValue('selectMotion',{[]},@iscell);

ip.parse(movieObj,varargin{:});
nBoot       = ip.Results.nBoot;
alpha       = ip.Results.alpha;
cluster     = ip.Results.cluster;
nCluster    = ip.Results.nCluster;
interval    = ip.Results.interval;
scale       = ip.Results.scale;
outputPath  = ip.Results.outputPath;
fileName    = ip.Results.fileName;

includeWin  = ip.Results.includeWin;
outLevel    = ip.Results.outLevel;
minLen      = ip.Results.minLength;
trend       = ip.Results.trendType;
gapSize     = ip.Results.gapSize;
winInterval = ip.Results.winInterval;
upPerc      = ip.Results.upPerc;
lwPerc      = ip.Results.lwPerc;
selection   = ip.Results.selectMotion;

if numel(interval) == 1
    interval = repmat(interval,1,nCell);
end

%% Formatting Time Series
operations = {'interval',interval,'outLevel',outLevel,'minLength',minLen,'trendType',trend,'gapSize',gapSize,'saveOn',false,'outputPath',outputPath,'fileName',fileName};
cellData   = formatMovieListTimeSeriesProcess(ML,'ProtrusionSamplingProcess',operations{:});


for iCell = 1:nCell
    
    sameWinInterval = false;
    sameIncludedWin = false;
    
    %% Setting up includedWin input if it's []
    
    
    if isempty( includeWin{iCell} )
       nWin              = size(cellData{iCell}.data.rawTimeSeries,1);
       includeWin{iCell} =  1:nWin;
    else
        nWin = numel(includeWin{iCell});
    end
    
    %% Setting up the winInterval input if it's [] for each cell
    
    if isempty(winInterval{iCell}{1})%If winInterval is []
        
        winInterval{iCell} = num2cell(repmat(1:cellData{iCell}.data.nFrames,nWin,1),2);
        
    else
        
        if numel(winInterval{iCell}) == 1 %If it's just one cell, repeat for all windows
            
            winInterval{iCell} = repmat(winInterval{iCell},1,nWin);
            
        end
        
    end
    
    %% Checking if the winInterval is the same as previously set
    
    if isfield(cellData{iCell}.data,'winInterval');
        
        sameWinInterval = all(cell2mat(cellfun(@(x,y) isequaln(x,y),winInterval{iCell},cellData{iCell}.data.winInterval,'Unif',0)));
        
        if ~sameWinInterval
            cellData{iCell}.data.winInterval = winInterval{iCell};
        end
    else
        
        cellData{iCell}.data.winInterval = winInterval{iCell};
        
    end
    
    %% Checking if the included windows are the same
    
    [currWin,winIdx] = setdiff(includeWin{iCell},cellData{iCell}.data.excludedWin{1});
    
    if isfield(cellData{iCell}.data,'includedWin');
        
        
        sameIncludedWin = isequaln(cellData{iCell}.data.includedWin{1},currWin);
        
        if ~sameIncludedWin
            cellData{iCell}.data.includedWin{1} = currWin;
        end
    else
        
        cellData{iCell}.data.includedWin{1} = currWin;
        
    end
    
    
    cellData{iCell}.data.analyzedLastRun = ~sameIncludedWin || ~sameWinInterval || cellData{iCell}.data.processedLastRun;
    
    %% If a different processing or different winInterval, then 
    if cellData{iCell}.data.analyzedLastRun
                
        scaling = 1;
        if scale
           scaling = ( cellData{iCell}.data.pixelSize/ cellData{iCell}.data.timeInterval);
        end
                
        cellData{iCell}.data.scaling           = scaling;
        cellData{iCell}.data.procEdgeMotion    = cellData{iCell}.data.procTimeSeries.*scaling;
        cellData{iCell}.data.procExcEdgeMotion = cellfun(@(win,time) cellData{iCell}.data.procEdgeMotion(win,time),...
                                                 num2cell(cellData{iCell}.data.includedWin{1}(:)),cellData{iCell}.data.winInterval(winIdx),'Unif',0);
    
    end
    
end

%% Getting Average Velocities and Persistence Time per Cell

runEdgeAnalysis = cellfun(@(x) x.data.analyzedLastRun,cellData);
commonGround    = @(x,z,y) mergingEdgeResults(x,'cluster',cluster,'nCluster',nCluster,'alpha',alpha,'nBoot',nBoot,'deltaT',z,'winInterval',y,'selection',selection);

if sum(runEdgeAnalysis) ~= 0
    
    %cellData(runEdgeAnalysis) = cellfun(@(x) rmfield(x,{'protrusionAnalysis','retractionAnalysis'}),cellData(runEdgeAnalysis),'Unif',0);
    firstLevel  = @(x,y,z) commonGround( cellfun(@(w) w(x),y,'Unif',0), z, {[]});
    secondLevel = @(x,y,z) cellfun(@(w) firstLevel(w,y,z),x,'Unif',0);
    
    [protrusion,retraction] ...
        = cellfun(@(x) secondLevel(x.data.interval,x.data.procExcEdgeMotion,x.data.timeInterval),cellData(runEdgeAnalysis),'Unif',0);

    [cellData,dataSet] = getDataSetAverage(cellData,protrusion,retraction,interval,lwPerc,upPerc,runEdgeAnalysis);    
    
    % Saving results
    savingMovieResultsPerCell(ML,cellData,outputPath,fileName)
    savingMovieDataSetResults(ML,dataSet,outputPath,fileName)
    
end

end%End of main function

function [cellData,dataSet] = getDataSetAverage(cellData,protrusion,retraction,interval,lwPerc,upPerc,idx)
%This function pull all the data from individual cells and calculates the dataSet mean value and CI
%It also formats the data structure for plotting
%Input:
%       cellData   - structure created by the function "formatEdgeVelocity.m".
%       interval   - Cell array containing the time intervals where the analysis is performed
%       alpha      - confidence interval Ex: - 0.05 = 95 percent
%       nBoot      - number of bootstrap samples
%
%Output:
%       dataSet    - averages and CI for the data set
%
nCell = numel(cellData);
cc    = 1;

for iCell = find(idx)
    
    
    for iiInt = 1:numel(interval{iCell})
        cellData{iCell}.protrusionAnalysis(iiInt) = protrusion{cc}{iiInt};
        cellData{iCell}.retractionAnalysis(iiInt) = retraction{cc}{iiInt};
    end
    cc = cc + 1;
    
end

dataSet = [];
nInter  = cellfun(@(x) numel(x),interval);

if sum(rem(nInter,nInter(1))) == 0
    
    for iInt = 1:numel(interval{1})
        
        total = struct('ProtPersTime',[],'ProtMaxVeloc',[],'ProtMinVeloc',[],'ProtMeanVeloc',[],'ProtMednVeloc',[],'RetrPersTime',[],'RetrMaxVeloc',[],'RetrMinVeloc',[],'RetrMeanVeloc',[],'RetrMednVeloc',[]);
        
        for iCell = 1:nCell
            
            total.ProtPersTime    = [total.ProtPersTime;cellData{iCell}.protrusionAnalysis(iInt).total.persTime];
            total.ProtMaxVeloc    = [total.ProtMaxVeloc;cellData{iCell}.protrusionAnalysis(iInt).total.maxVeloc];
            total.ProtMinVeloc    = [total.ProtMinVeloc;cellData{iCell}.protrusionAnalysis(iInt).total.minVeloc];
            total.ProtMeanVeloc   = [total.ProtMeanVeloc;cellData{iCell}.protrusionAnalysis(iInt).total.Veloc];
            total.ProtMednVeloc   = [total.ProtMednVeloc;cellData{iCell}.protrusionAnalysis(iInt).total.mednVeloc];
            
            total.RetrPersTime    = [total.RetrPersTime;cellData{iCell}.retractionAnalysis(iInt).total.persTime];
            total.RetrMaxVeloc    = [total.RetrMaxVeloc;cellData{iCell}.retractionAnalysis(iInt).total.maxVeloc];
            total.RetrMinVeloc    = [total.RetrMinVeloc;cellData{iCell}.retractionAnalysis(iInt).total.minVeloc];
            total.RetrMeanVeloc   = [total.RetrMeanVeloc;cellData{iCell}.retractionAnalysis(iInt).total.Veloc];
            total.RetrMednVeloc   = [total.RetrMednVeloc;cellData{iCell}.retractionAnalysis(iInt).total.mednVeloc];
            
        end
        
        dataSet.CI.interval(iInt  )        = structfun(@(x) prctile(x,[lwPerc upPerc]),total,'Unif',0);
        dataSet.medianValue.interval(iInt) = structfun(@(x) nanmedian(x),total,'Unif',0);
        dataSet.total.interval(iInt)       = total;
        
    end
    
end

end