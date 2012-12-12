function cellData = excludeWindowsFromAnalysis(movieObj)

excBorder = false;

if isa(movieObj,'MovieData')
    
    ML = movieData2movieList(movieObj);
    
else
    
    ML = movieObj;
    
end

nCell    = numel(ML.movies_);
exclude  = cell(nCell,1);
cellData = loadingMovieResultsPerCell(ML);


%exclude{2} = [1:35 95:112];

for iCell = 1:nCell
    
    if excBorder
        nWin   = size( cellData(iCell).data.rawEdgeMotion,1 );
        border = [1:3 nWin-2:nWin];
        cellData(iCell).data.excludeWin = unique([cellData(iCell).data.excludeWin border]);
        %Excluding pre-selected windows
        cellData(iCell).data.excludeWin = unique([cellData(iCell).data.excludeWin exclude{iCell}]);        
    end
    
    cellData(iCell).data.procEdgeMotion(cellData(iCell).data.excludeWin) = [];

    
end


savingMovieResultsPerCell(ML,cellData)