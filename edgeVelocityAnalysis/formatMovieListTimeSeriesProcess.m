function cellData = formatMovieListTimeSeriesProcess(movieObj,processType,varargin)
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

ip.addParamValue('channel',   0,@isscalar);
ip.addParamValue('includeWin',cell(1,nCell),@iscell);
ip.addParamValue('outLevel',  0,@isscalar);
ip.addParamValue('trendType',-1,@isscalar);
ip.addParamValue('minLength', 30,@isscalar);
ip.addParamValue('gapSize',   0,@isscalar);
ip.addParamValue('saveOn',    false,@islogical);

ip.parse(movieObj,processType,varargin{:});
includeWin = ip.Results.includeWin;
outLevel   = ip.Results.outLevel;
minLen     = ip.Results.minLength;
trend      = ip.Results.trendType;
gapSize    = ip.Results.gapSize;
channel    = ip.Results.channel;
saveOn     = ip.Results.saveOn;

timeSeriesOperations = {'outLevel',outLevel,'minLength',minLen,'trendType',trend,'gapSize',gapSize};


dataS    = struct('includedWin',[],'excludedWin',[],'pixelSize',[],'frameRate',1,'rawTimeSeries',[],'procTimeSeries',[]);
cellData = struct('data',repmat({dataS},1,nCell)) ;

for iCell = 1:nCell
    
    currMD     = ML.movies_{iCell};
    timeSeries = readingTimeSeries(currMD,formattableProc,processType,channel);
    nDim       = ndims(timeSeries);
    nWin       = size(timeSeries,1);
    
    if isempty(includeWin{iCell})
        
        includeWin{iCell} = 1:nWin;
        
    end
    
    %Applying Time Series Operations
    cellData(iCell).data.timeSeriesOperations = timeSeriesOperations;
    if nDim == 3
        
        for iLayer = 1:size(timeSeries,2)
            
            [cellData(iCell).data.procTimeSeries(:,:,iLayer),excludeVar] = timeSeriesPreProcessing(cellData(iCell).data.rawTimeSeries,timeSeriesOperations{:});
            cellData(iCell).data.excludedWin{iLayer}                     = unique([setdiff(1:nWin,includeWin{iCell}) excludeVar]);
            cellData(iCell).data.includedWin{iLayer}                     = setdiff(includeWin{iCell},excludeVar);
            
        end
        
    else
        cellData(iCell).data.rawTimeSeries               = timeSeries; 
        [cellData(iCell).data.procTimeSeries,excludeVar] = timeSeriesPreProcessing(cellData(iCell).data.rawTimeSeries,timeSeriesOperations{:});
        cellData(iCell).data.excludedWin                 = unique([setdiff(1:nWin,includeWin{iCell}) excludeVar]);
        cellData(iCell).data.includedWin                 = setdiff(includeWin{iCell},excludeVar);
        
    end
    
end
% Performing windowing exclusion
%cellData = excludeWindowsFromAnalysis(ML,'excBorder',border,'cellData',cellData);

%% Saving results per cell
if saveOn
    
    savingMovieResultsPerCell(ML,cellData)
    
end

end

function samples = readingTimeSeries(currMD,forProc,processType,channel)

procIdx = currMD.getProcessIndex(processType);

if strcmp(processType,forProc{1})
    
    samples = currMD.processes_{procIdx}.loadChannelOutput.avgNormal;
    
elseif strcmp(processType,forProc{2})
    
    samples = currMD.processes_{procIdx}.loadChannelOutput(channel).avg;
    
end

end