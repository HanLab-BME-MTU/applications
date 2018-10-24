function cellData = excludeWindowsFromAnalysis(movieObj,varargin)
% This function exclude selected windows from analysis.
%
%Input:
%       movieObj - movie data or movie list object
%       excBorder - if only a region of the cell is windowed, this option
%                 allows the user to exclude the border windows
%                 - scalar that defines the number of windows in both
%                 borders (default = 0)
%
%       excludeW - cell array that contains the index of windows to be
%                  excluded for each movie.
%                  Ex: excludeW{1} = [4 6]       - excluding windows 4 and
%                      excludeW{2} = [1 9 8 4]
%
%Output:
%       cellData - structure with the updated processed data
%                  This function will write a file for each movie - see formatEdgeVelocity
%
% See also: formatEdgeVelocity, edgeVelocityQuantification
%
%Marco Vilela, 2012

ip = inputParser;
ip.addRequired('movieObj',@(x) isa(x,'MovieList') || isa(x,'MovieData'));
ip.addParamValue('excBorder',0,@isscalar);
ip.addParamValue('excludeW',cell(1),@iscell);
ip.addParamValue('cellData',[],@isstruct);

ip.parse(movieObj,varargin{:});
excBorder = ip.Results.excBorder;
excludeW  = ip.Results.excludeW;
cellData  = ip.Results.cellData;

if isa(movieObj,'MovieData')
    
    ML = movieData2movieList(movieObj);
    
else
    
    ML = movieObj;
    
end

if isempty(cellData)
    %Loading formatted data - see formatEdgeVelocity
    analysis = 'FormatEdgeVelocity'; [~,fileName]=fileparts(movieObj.getPath);
    cellData = loadingMovieResultsPerCell(ML,analysis,fileName);
end

nCell   = numel(ML.movies_);
nExc    = numel(excludeW);
exclude = cell(nCell,1);
exclude(1:nExc) = excludeW;


for iCell = 1:nCell
    curData = cellData(iCell);
    if iscell(curData)
        curData=curData{1};
    end
        
    if excBorder
        nWin   = size( curData.data.rawEdgeMotion,1 );
        border = [1:excBorder nWin-(excBorder-1):nWin];
        curData.data.excludeWin = unique([curData.data.excludeWin border]);
    end
    
    %Excluding pre-selected windows
    curData.data.excludedWin = unique([curData.data.excludedWin exclude{iCell}]);
    
    curData.data.excProcEdgeMotion                                  = num2cell(curData.data.procEdgeMotion,2);
    curData.data.excProcEdgeMotion(curData.data.excludedWin) = [];
    try
        cellData(iCell) = curData;
    catch
        cellData{iCell} = curData;
    end
end

analysis = 'ExcludedWindows'; [~,fileName]=fileparts(movieObj.getPath);
savingMovieResultsPerCell(ML,cellData,analysis,fileName)