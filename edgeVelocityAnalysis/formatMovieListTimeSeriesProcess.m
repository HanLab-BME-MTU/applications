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

ip.addParamValue('channel',   0,@isscalar);
ip.addParamValue('includeWin',cell(1,nCell), @iscell);
ip.addParamValue('winInterval',cell(1,nCell),@iscell);
ip.addParamValue('outLevel',  zeros(1,nCell),@isvector);
ip.addParamValue('trendType',-ones(1,nCell),@isvector);
ip.addParamValue('minLength', 30*ones(1,nCell),@isvector);
ip.addParamValue('gapSize',   zeros(1,nCell),@isvector);
ip.addParamValue('saveOn',    false,@islogical);
ip.addParamValue('outputPath','edgeVelocityQuantification',@isstr);
ip.addParamValue('fileName','edgeVelocity',@isstr);


ip.parse(movieObj,processType,varargin{:});
includeWin  = ip.Results.includeWin;
winInterval = ip.Results.winInterval;
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
    
    
    timeSeriesOperations{iCell} = {'outLevel',outLevel(iCell),'minLength',minLen(iCell),'trendType',trend(iCell),'gapSize',gapSize(iCell)};
    
    inp1 = isempty(cellData{iCell});
    inp3 = false;
    
    if inp1
       
        currMD     = ML.movies_{iCell};
        timeSeries = readingTimeSeries(currMD,formattableProc,processType,channel);
        nDim       = ndims(timeSeries);
        
        if nDim == 3
            cellData{iCell}.data.rawTimeSeries = permute(timeSeries,[1 3 2]);%Reverting Hunter's retarded decision
        else
            cellData{iCell}.data.rawTimeSeries = timeSeries;
        end

    else
        
        inp3 = isequal(cellData{iCell}.data.timeSeriesOperations,timeSeriesOperations{iCell});
        inp2 = false;
        inp5 = false;
        if ~isempty(winInterval{iCell}{1})
            
            if numel(cellData{iCell}.data.winInterval) == numel(winInterval{iCell})%Same number of intervals
              
                inp2 = sum( cellfun(@(x,y) isequaln(x,y),cellData{iCell}.data.winInterval(:),winInterval{iCell}(:)) ) == numel(winInterval{iCell});
                
            end
            
        end
        
        if ~isempty(includeWin{iCell})
        
            inp5 = isequal(cellData{iCell}.data.includedWin{1},includeWin{iCell});
        
        end
        
        if isfield(cellData{iCell}.data,'rawEdgeMotion')
            
            cellData{iCell}.data.rawTimeSeries  = cellData{iCell}.data.rawEdgeMotion;
            cellData{iCell}.data.procTimeSeries = cellData{iCell}.data.procEdgeMotion;
            
            if ~(inp3 && inp2 && inp5) && isfield(cellData{iCell},'protrusionAnalysis')
                cellData{iCell} = rmfield(cellData{iCell},{'protrusionAnalysis','retractionAnalysis'});
            end
            
        elseif isfield(cellData{iCell}.data,'rawSignal')
            
            cellData{iCell}.data.rawTimeSeries  = cellData{iCell}.data.rawSignal;
            cellData{iCell}.data.procTimeSeries = cellData{iCell}.data.procSignal;
            
            if ~( inp3 && inp2 && inp5)
                cellData{iCell} = rmfield(cellData{iCell},{'intensityOverTime','intensityOverTimeSpace'});
            end
            
        end
        
    end
    
    
    [nWin,nObs,nLayer] = size(cellData{iCell}.data.rawTimeSeries);
    excludeVar         = cell(1,nLayer);
    %%  If TS operations are different, process data with new settings
    if inp1 || ~inp3
        
        %Applying Time Series Operations
        cellData{iCell}.data.timeSeriesOperations = timeSeriesOperations{iCell};
        
        for iLayer = 1:nLayer
            
            [cellData{iCell}.data.procTimeSeries(:,:,iLayer),excludeVar{iLayer}] = timeSeriesPreProcessing(squeeze(cellData{iCell}.data.rawTimeSeries(:,:,iLayer)),timeSeriesOperations{iCell}{:});
            
        end
        
    end
    
    
    if isempty(includeWin{iCell})%If includeWin is [], include all windows
        
        includeWin{iCell} = 1:nWin;
        
    end
    
    nWin = numel(includeWin{iCell});
    if isempty(winInterval{iCell}{1})%If winInterval is [], include all time points for all windows
        
        cellData{iCell}.data.winInterval = num2cell(repmat(1:nObs,nWin,1),2);
        
    else
        
        if sum(cellfun(@(x,y) numel(x)-numel(y),winInterval,includeWin)) ~= 0
            error('Number of windows does not match number of intervals')
        end

        cellData{iCell}.data.winInterval = winInterval{iCell};
        
    end
    
    
    for iLayer = 1:nLayer
        
        cellData{iCell}.data.excludedWin{iLayer}       = unique([setdiff(1:size(cellData{iCell}.data.rawTimeSeries,2),includeWin{iCell}) excludeVar{iLayer}]);
        cellData{iCell}.data.includedWin{iLayer}       = setdiff(includeWin{iCell},excludeVar{iLayer});
        
        if numel(cellData{iCell}.data.includedWin{iLayer}) ~= numel(cellData{iCell}.data.winInterval)
            
            cellData{iCell}.data.procExcTimeSeries{iLayer} = {[]};
            
        else
            
            cellData{iCell}.data.procExcTimeSeries{iLayer} = cellfun(@(win,time) cellData{iCell}.data.procTimeSeries(win,time,iLayer),...
                                                             num2cell(cellData{iCell}.data.includedWin{iLayer}(:)),cellData{iCell}.data.winInterval(:),'Unif',0);
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
    
    samples = currMD.processes_{procIdx}.loadChannelOutput(channel).avg;
    
end

end