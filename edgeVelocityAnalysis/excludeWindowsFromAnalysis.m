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

ip.parse(movieObj,varargin{:});
excBorder = ip.Results.excBorder;
excludeW  = ip.Results.excludeW;

if isa(movieObj,'MovieData')
    
    ML = movieData2movieList(movieObj);
    
else
    
    ML = movieObj;
    
end

nCell   = numel(ML.movies_);
nExc    = numel(excludeW);
exclude = cell(nCell,1);
exclude(1:nExc) = excludeW;

%Loading formatted data - see formatEdgeVelocity
cellData = loadingMovieResultsPerCell(ML);

for iCell = 1:nCell
    
    if excBorder
        nWin   = size( cellData(iCell).data.rawEdgeMotion,1 );
        border = [1:excBorder nWin-(excBorder-1):nWin];
        cellData(iCell).data.excludeWin = unique([cellData(iCell).data.excludeWin border]);
    end
    
     %Excluding pre-selected windows
    cellData(iCell).data.excludeWin = unique([cellData(iCell).data.excludeWin exclude{iCell}]);  
    
    cellData(iCell).data.excProcEdgeMotion                                  = cellData(iCell).data.procEdgeMotion;
    cellData(iCell).data.excProcEdgeMotion(cellData(iCell).data.excludeWin) = [];

end

savingMovieResultsPerCell(ML,cellData)