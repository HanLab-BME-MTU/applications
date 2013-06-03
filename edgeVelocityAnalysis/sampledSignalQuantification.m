function cellData = sampledSignalQuantification(movieObj,channel,varargin)
% This function quantifies:
%                           mean protrusion/retraction instantaneous values and confidence interval;
%                           mean protrusion/retraction persistence time and confidence interval;
%                           cluster velocities and persistence time
%                           for each window: persistence time;protrusion/retraction block;Upper and lower noise threshold
%
% All the above quantifications are done for each cell of the movieObj (movieList or movieData)
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
ip.addRequired('channel',@isscalar);
ip.addParamValue('nBoot',1e3,@isscalar);
ip.addParamValue('alpha',.05,@isscalar);
ip.addParamValue('cluster',false,@isscalar);
ip.addParamValue('nCluster',2,@isscalar);
ip.addParamValue('interval',{[]},@iscell);

if isa(movieObj,'MovieData')
    
    ML = movieData2movieList(movieObj);
    
else
    
    ML = movieObj;
    
end

nCell = numel(ML.movies_);

%% Time Series Pre-Processing operations
ip.addParamValue('includeWin', cell(1,nCell),@iscell);
ip.addParamValue('outLevel',0,@isscalar);
ip.addParamValue('trendType',   -1,@isscalar);
ip.addParamValue('minLength',  10,@isscalar);


ip.parse(movieObj,channel,varargin{:});
nBoot    = ip.Results.nBoot;
alpha    = ip.Results.alpha;
interval = ip.Results.interval;


includeWin = ip.Results.includeWin;
outLevel   = ip.Results.outLevel;
minLen     = ip.Results.minLength;
trend      = ip.Results.trendType;

%% Formatting Time Series
operations = {'channel',channel,'includeWin',includeWin,'outLevel',outLevel,'minLength',minLen,'trendType',trend};
cellData   = formatMovieListTimeSeriesProcess(ML,'WindowSamplingProcess',operations{:});


for iCell = 1:nCell
    
    currMD                          = ML.movies_{iCell};
    cellData(iCell).data.pixelSize  = currMD.pixelSize_;
    cellData(iCell).data.frameRate  = currMD.timeInterval_;
    cellData(iCell).data.rawSignal  = cellData(iCell).data.rawTimeSeries;
    cellData(iCell).data.procSignal = cellData(iCell).data.procTimeSeries;
    cellData(iCell).data            = rmfield(cellData(iCell).data,{'rawTimeSeries','procTimeSeries'});
    
    [nWin,~,nLayer] = size(cellData(iCell).data.procSignal);
    for iLayer = 1:nLayer
        
        windows = cellData(iCell).data.includedWin{iLayer};
        signal  = squeeze( cellData(iCell).data.procSignal(:,:,iLayer) );
        
        cellData(iCell).intensityOverTime(iLayer,:)       = nan(1,nWin);
        cellData(iCell).intensityOverTime(iLayer,windows) = nanmean( signal(windows,:),2 );
        cellData(iCell).intensityOverTimeSpace(iLayer)    = nanmean( cellData(iCell).intensityOverTime(iLayer,windows) );
    
    end
    
end


%% Saving results

savingMovieResultsPerCell(ML,cellData,'sampledSignalQuantification')

end