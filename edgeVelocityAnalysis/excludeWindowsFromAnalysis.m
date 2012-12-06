function cellData = excludeWindowsFromAnalysis(movieObj)


if isa(movieObj,'MovieData')
    
    ML = movieData2movieList(movieObj);
    
else
    
    ML = movieObj;
    
end

nCell    = numel(ML.movies_);
exclude  = cell(nCell,1);
cellData = loadingMovieResultsPerCell(ML);


%exclude{1} = [1:35 95:112];

for iCell = 1:nCell
    
    nWin   = size( cellData(iCell).data.rawEdgeMotion,1 );
    border = [1:3 nWin-2:nWin];
    %Excluding border
    cellData(iCell).data.excludeWin = unique([cellData(iCell).data.excludeWin border]);
    
    %Excluding pre-selected windows
    cellData(iCell).data.excludeWin = unique([cellData(iCell).data.excludeWin exclude{iCell}]);
    
end

savingMovieResultsPerCell(ML,cellData)