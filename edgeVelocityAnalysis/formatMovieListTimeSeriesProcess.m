function [cellData,dataSet] = formatMovieListTimeSeriesProcess(movieObj,processType,varargin)
%This function takes the output of the protrusion sampling process and formats each edge velocity time series
%Format actually means TS pre-processing. It removes: outliers, mean, trend, NaN and close gaps
%
% Usage: cellData = formatMovieListTimeSeriesProcess(ML,varargin)
%
% INPUTS:
%       ML - movie list or movie data object
%
%       includeWin - cell array with the same size as the ML. Each element has the indexes of the windows(variables) to be included for analysis;
%                   All windows that are not in this array will be excluded. The default value is {[]}, which includes all windows.
%
%       outLevel  - # of sigmas considered for outlier removal (see detectOutliers)
%
%       trendType - optional: a scalar giving the type of trend to remove. If no value is inputed, nothing is done.
%                    0 - remove only the mean value
%                    1 - remove linear trend
%                    2 - remove exponential trend
%                    3 - remove double exponential trend
%                    4 - remove nonlinear local trend (trendFilteringEMD)
%                    5 - remove all determinitic component and spits out a stationary signal
%                        the trend in this case is a smoothed version of the input signal
%
%       minLength  - minimal length accepted. Any window that has a TS with less than
%                    minLength points will be discarded. (Default = 30)
%
%
%
%Output:
%       cellData - structure for each cell with the TS operation results
%
%       This function creates a folder (EdgeVelocityAnalysis) and writes a file edgeVelocity.mat
%
%See also: excludeWindowsFromAnalysis, edgeVelocityQuantification
%Marco Vilela, 2012

formattableProc = {'ProtrusionSamplingProcess','WindowSamplingProcess'};
testInput       = @(y) any(cell2mat(cellfun(@(x) strcmp(y,x),formattableProc,'Unif',0)));
ip = inputParser;
ip.addRequired('movieObj',@(x) isa(x,'MovieList') || isa(x,'MovieData'));
ip.addRequired('processType',@(x) testInput(x));


if isa(movieObj,'MovieData')
    
    ML = movieData2movieList(movieObj);
    
else
    
    ML = movieObj;
    
end

nCell = numel(ML.movies_);

ip.addParamValue('channel',   1,@isscalar);
ip.addParamValue('interval',num2cell(cell(1,nCell)),@iscell);
ip.addParamValue('outLevel',  zeros(1,nCell),@isvector);
ip.addParamValue('trendType',-ones(1,nCell),@isvector);
ip.addParamValue('minLength', 30*ones(1,nCell),@isvector);
ip.addParamValue('gapSize',   zeros(1,nCell),@isvector);
ip.addParamValue('saveOn',    false,@islogical);
ip.addParamValue('outputPath','EdgeVelocityQuantification',@isstr);
ip.addParamValue('fileName','EdgeMotion',@isstr);


ip.parse(movieObj,processType,varargin{:});
interval    = ip.Results.interval;
outLevel    = ip.Results.outLevel;
minLen      = ip.Results.minLength;
trend       = ip.Results.trendType;
gapSize     = ip.Results.gapSize;
channel     = ip.Results.channel;
saveOn      = ip.Results.saveOn;
outputPath  = ip.Results.outputPath;
fileName    = ip.Results.fileName;


dataSet                       = [];
timeSeriesOperations{1,nCell} = [];
cellData                      = loadingMovieResultsPerCell(ML,outputPath,fileName);

for iCell = 1:nCell
    
    %Building interval input in case it's []
    timeSeriesOperations{iCell} = {'outLevel',outLevel(iCell),'minLength',minLen(iCell),'trendType',trend(iCell),'gapSize',gapSize(iCell)};
    
    %Logical check variables
    hasNotBeenProc = isempty(cellData{iCell});
    sameTSoperat   = false;
    sameInterval   = false;
    
    %Checking interval input
    intervalEmpty  = all(cellfun(@(x) isempty(x),interval{iCell}));    
    
    if ~hasNotBeenProc
        
        if cellData{iCell}.data.channel ~= channel
            hasNotBeenProc = true;
        end
    end
    
    if hasNotBeenProc % If data has not been processed
       
        currMD     = ML.movies_{iCell};
        timeSeries = readingTimeSeries(currMD,formattableProc,processType,channel);
        nDim       = ndims(timeSeries);
        
        if nDim == 3
            cellData{iCell}.data.rawTimeSeries = permute(timeSeries,[1 3 2]);%Reverting Hunter's retarded decision
        else
            cellData{iCell}.data.rawTimeSeries = timeSeries;
        end
        if isempty(currMD.pixelSize_) || isempty(currMD.timeInterval_)
            error(['Movie ' num2str(iCell) 'does not have the pixel size and/or time interval setup'])
        end
        
        cellData{iCell}.data.nFrames           = currMD.nFrames_;
        cellData{iCell}.data.pixelSize         = currMD.pixelSize_;
        cellData{iCell}.data.timeInterval      = currMD.timeInterval_;
        cellData{iCell}.data.channel           = channel;
        
        [~,nObs,nLayer] = size(cellData{iCell}.data.rawTimeSeries);
        
        %Building interval input in case it's []
        if intervalEmpty
            cellData{iCell}.data.interval = {1:nObs};
        else
            cellData{iCell}.data.interval = interval{iCell};
        end
        
    else
        
        [~,nObs,nLayer] = size(cellData{iCell}.data.(['proc' fileName]));
        % Comparing the input TS operations with the previous TS processing 
        sameTSoperat = isequal(cellData{iCell}.data.timeSeriesOperations,timeSeriesOperations{iCell});
        
        %Building interval input in case it's []
        if intervalEmpty
            interval{iCell} = {1:nObs};
        end
        
        %Comparing the intervals
        if numel(cellData{iCell}.data.interval) == numel(interval{iCell})
            sameInterval = all( cellfun(@(x,y) isequaln(x,y),cellData{iCell}.data.interval(:),interval{iCell}(:)) );
            if ~sameInterval
                cellData{iCell}.data.interval = interval{iCell};
            end
        end
        
    end
    
    %%  If TS operations are different, process data with new settings
    
    excludeVar = cell(1,nLayer);
      
    cellData{iCell}.data.processedLastRun = hasNotBeenProc || ~sameTSoperat || ~sameInterval;
    
    
    if cellData{iCell}.data.processedLastRun
        
        %Applying Time Series Operations - applied for each interval
        cellData{iCell}.data.timeSeriesOperations = timeSeriesOperations{iCell};
        preProcessingInput                        = [timeSeriesOperations{iCell} {'interval',cellData{iCell}.data.interval}];
        
        for iLayer = 1:nLayer
            
            [cellData{iCell}.data.procTimeSeries(:,:,iLayer),excludeVar{iLayer}] = timeSeriesPreProcessing(squeeze(cellData{iCell}.data.rawTimeSeries(:,:,iLayer)),preProcessingInput{:});
            cellData{iCell}.data.excludedWin{iLayer}                             = excludeVar{iLayer};
            
        end
        
    end
    
    
end


%% Saving results per cell
if saveOn
    
    savingMovieResultsPerCell(ML,cellData,outputPath,fileName)
    
end

end

function samples = readingTimeSeries(currMD,forProc,processType,channel)

procIdx = currMD.getProcessIndex(processType);

if strcmp(processType,forProc{1})
    
    samples = currMD.processes_{procIdx}.loadChannelOutput.avgNormal;
    
elseif strcmp(processType,forProc{2})
    
    samples = currMD.processes_{procIdx}.loadChannelOutput(1,channel).avg;
    
end

end