function cellData = formatEdgeVelocity(movieObj,varargin)
%This function takes the output of the protrusion sampling process and formats each edge velocity time series
%Format actually means TS pre-processing. It removes: outliers, mean, trend, NaN and close gaps 
%
% Usage: cellData = formatEdgeVelocity(ML,varargin)
%
% INPUTS:
%       ML - movie list or movie data object  
%
%       includeWin - cell array with the same size as the ML. Each element has the indexes of the windows(variables) to be included for analysis;
%                   All windows that are not in this array will be excluded. The default value is {inf}, which includes all windows.
%
%       outLevel  - # of sigmas considered for outlier removal (see detectOutliers)
%
%       trendType - optional: a scalar giving the type of trend to remove
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
%       scale      - convert velocity from pixel/frame into nm/sec
%
%       excBorder  - exclude border windows. This parameter is actually the length of the border
%
%Output:
%       cellData - structure for each cell with the TS operation results
%                   
%       This function creates a folder (EdgeVelocityAnalysis) and writes a file edgeVelocity.mat
%
%See also: excludeWindowsFromAnalysis, edgeVelocityQuantification
%Marco Vilela, 2012


ip = inputParser;
ip.addRequired('movieObj',@(x) isa(x,'MovieList') || isa(x,'MovieData'));

if isa(movieObj,'MovieData')
    
    ML = movieData2movieList(movieObj);
    
else
    
    ML = movieObj;
    
end
nCell = numel(ML.movies_);

ip.addParamValue('includeWin', num2cell(inf(1,nCell)),@iscell);
ip.addParamValue('excBorder',0,@isscalar);
ip.addParamValue('outLevel',0,@isscalar);
ip.addParamValue('trendType',   -1,@isscalar);
ip.addParamValue('minLength',  30,@isscalar);
ip.addParamValue('scale', false,@islogical);
ip.addParamValue('gapSize',0,@isscalar);

ip.parse(movieObj,varargin{:});
includeWin = ip.Results.includeWin;
outLevel   = ip.Results.outLevel;
minLen     = ip.Results.minLength;
scale      = ip.Results.scale;
trend      = ip.Results.trendType;
gapSize    = ip.Results.gapSize;
excBorder    = ip.Results.excBorder;

timeSeriesOperations = {'outLevel',outLevel,'minLength',minLen,'trendType',trend,'gapSize',gapSize};


dataS    = struct('includedWin',[],'excludedWin',[],'pixelSize',[],'frameRate',[],'rawEdgeMotion',[],'procEdgeMotion',[]);
cellData = struct('data',repmat({dataS},1,nCell)) ;

for iCell = 1:nCell
    
    curMD                             = ML.movies_{iCell};
    cellData(iCell).data.includedWin   = includeWin{iCell};
    % Get the windowing package
    windPack = curMD.getPackage(curMD.getPackageIndex('WindowingPackage'));
    % Protrusion sampling process
    protProc = windPack.processes_{3};
    protSamples = protProc.loadChannelOutput;
    cellData(iCell).data.rawEdgeMotion = protSamples.avgNormal;
    
    %Extracting outliers
    %Removing NaN and closing 1 frame gaps
    cellData(iCell).data.timeSeriesOperations        = timeSeriesOperations;
    [cellData(iCell).data.procEdgeMotion,excludeVar] = timeSeriesPreProcessing(cellData(iCell).data.rawEdgeMotion,timeSeriesOperations{:});    
    cellData(iCell).data.excludedWin                  = unique([cellData(iCell).data.excludedWin excludeVar]);
    
end
%% Performing windowing exclusion
cellData = excludeWindowsFromAnalysis(ML,'excBorder',excBorder,'cellData',cellData);

%% Saving results per cell
analysis = 'FormatEdgeVelocity'; [~,fileName]=fileparts(movieObj.getPath);
savingMovieResultsPerCell(ML,cellData,analysis,fileName)

end