function cellData = sampledSignalQuantification(movieObj,channel,varargin)
% This function quantifies:
%                           mean protrusion/retraction instantaneous values and confidence interval;
%                           mean protrusion/retraction persistence time and confidence interval;
%                           cluster velocities and persistence time
%                           for each window: persistence time;protrusion/retraction block;Upper and lower noise threshold
%
% All the above quantifications are do,1,signalInput{:}ne for each cell of the movieObj (movieList or movieData)
% If the input 'interval' is set, all the described quantifications are calculated for each time interval
%
% IMPORTANT: This function requires to format the time series first. To do that, run the function "formatEdgeVelocity".
%
%Usage: [cellData,dataSet] = edgeVelocityQuantification(movieObj,varargin)
%
% Input:
%       nBoot      - # of boostrap samples to be used (default value 1000)
%
%       alpha      - alpha used to generate the bootstrap confidence intervals
%                    (default value 0.05)
%
%       cluster    - scalar "1" to perform cluster analysis; "0" otherwise
%
%       nCluster   - number of cluster
%
%       interval   - cell array with the time intervals in frames
%                     Ex : interval = {[1:30],[20:40]} - Quantification will be done with velocities calculated at each element of the cell array
%
%       scale      - convert velocity from pixel/frame into nm/sec
%
%       includeWin -
%
%       outLevel   -
%
%       trendType  -
%
%       minLength  -
%
%       scale      -
%
% Output:
%       cellData - this is a long structure array cellData(for each cell) with the following fields:
%                .data.excludeWin - indexes of windows that were excluded from analysis . Ex: border windows
%                     .rawEdgeMotion - ,1,signalInput{:}raw edge velocity time series. Protrusion map. This is redundant.
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
%                                      ,1,signalInput{:}                CI.persTimeCI
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
%                                                                               .MednV,'channel',1eloc
%
%Marco Vilela, 2012


ip = inputParser;
ip.addRequired('movieObj',@(x) isa(x,'MovieList') || isa(x,'MovieData'));
ip.addRequired('channel',@isscalar);
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
ip.addParamValue('minLength', 10*ones(1,nCell),@isvector);
ip.addParamValue('gapSize',   zeros(1,nCell),@isvector);
ip.addParamValue('outputPath','sampledSignalQuantification',@isstr);
ip.addParamValue('fileName','Signal',@isstr);
ip.addParamValue('interval',num2cell(cell(1,nCell)),@iscell);

ip.parse(movieObj,channel,varargin{:});

includeWin  = ip.Results.includeWin;
outLevel    = ip.Results.outLevel;
minLen      = ip.Results.minLength;
trend       = ip.Results.trendType;
outputPath  = ip.Results.outputPath;
fileName    = ip.Results.fileName;
gapSize     = ip.Results.gapSize;
winInterval = ip.Results.winInterval;

%% Formatting Time Series
operations = {'channel',channel,'outLevel',outLevel,'minLength',minLen,'trendType',trend,'gapSize',gapSize,'outputPath',outputPath,'fileName',fileName};
cellData   = formatMovieListTimeSeriesProcess(ML,'WindowSamplingProcess',operations{:});

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
            
            winInterval{iCell} = repmat(winInterval{iCell},nWin,1);
            
        end
        
    end
    
    
    %% Checking if the included windows are the same
    
    [currWin] = setdiff(includeWin{iCell},cellData{iCell}.data.excludedWin{1});
    
    if isfield(cellData{iCell}.data,'includedWin');
        
        
        sameIncludedWin = isequaln(cellData{iCell}.data.includedWin{1},currWin);
        
        if ~sameIncludedWin
            cellData{iCell}.data.includedWin{1} = currWin;
        end
    else
        
        cellData{iCell}.data.includedWin{1} = currWin;
        
    end

    %% Checking if the winInterval is the same as previously set
    
     if isfield(cellData{iCell}.data,'winInterval') && sameIncludedWin
        
        if numel(winInterval{iCell}) == numel(cellData{iCell}.data.winInterval)
            
            sameWinInterval = all(cell2mat(cellfun(@(x,y) isequaln(x,y),winInterval{iCell},cellData{iCell}.data.winInterval,'Unif',0)));
            
            if ~sameWinInterval
                cellData{iCell}.data.winInterval = winInterval{iCell};
            end
        else
            cellData{iCell}.data.winInterval = winInterval{iCell};
        end
        
    else
        
        cellData{iCell}.data.winInterval = winInterval{iCell};
        
    end
    
    cellData{iCell}.data.analyzedLastRun   = ~sameIncludedWin || ~sameWinInterval || cellData{iCell}.data.processedLastRun;
    
    %%
    
    if cellData{iCell}.data.analyzedLastRun
        
        cellData{iCell}.data.rawSignal     = cellData{iCell}.data.rawTimeSeries;
        cellData{iCell}.data.procSignal    = cellData{iCell}.data.procTimeSeries;
        
       
        
        nLayer = size(cellData{iCell}.data.procSignal,3);
        for iLayer = 1:nLayer
            
            [currWin,winIdx] = setdiff(includeWin{iCell},cellData{iCell}.data.excludedWin{iLayer});
            
            if isfield(cellData{iCell}.data,'includedWin');
                
                pLayer = numel(cellData{iCell}.data.includedWin);
                
                if pLayer >= iLayer
                    
                    sameIncludedWin = isequaln(cellData{iCell}.data.includedWin{iLayer},currWin);
                    
                    if ~sameIncludedWin
                        cellData{iCell}.data.includedWin{iLayer} = currWin;
                    end
                else
                    cellData{iCell}.data.includedWin{iLayer} = currWin;
                end
                
            else
                
                cellData{iCell}.data.includedWin{iLayer} = currWin;
                
            end
            
            
            cellData{iCell}.data.procExcSignal{iLayer} = cellfun(@(win,time) cellData{iCell}.data.procSignal(win,time,iLayer),...
                                                         num2cell(cellData{iCell}.data.includedWin{iLayer}(:)),cellData{iCell}.data.winInterval(winIdx),'Unif',0);
                                         
            windows = cellData{iCell}.data.includedWin{iLayer};
            signal  = cellData{iCell}.data.procExcSignal{iLayer};
            
            cellData{iCell}.intensityOverTime(1:numel(windows),iLayer) = cellfun(@(x) nanmean(x),signal);
            cellData{iCell}.intensityOverTimeSpace(iLayer)             = nanmean( cellData{iCell}.intensityOverTime(:,iLayer) );
            
        end
        
    end
    
end

%% Saving results

savingMovieResultsPerCell(ML,cellData,outputPath,fileName)

end