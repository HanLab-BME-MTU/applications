function  out = loadingMovieResultsPerCell(ML,analysisDir,fileName)
%Loading edge velocity analysis for each cell in the MovieList

nCell = numel(ML.movies_);
cF    = 1;

for iCell = 1:nCell

    currMD = ML.movies_{iCell};
     %Saving results for each cell
    edgePath = [currMD.outputDirectory_ filesep analysisDir];
    filePath = [edgePath filesep fileName '.mat'];
    if exist(filePath,'file')
        
        aux      = load(filePath);
        out(cF)  = aux.analysisResults;
        cF       = cF + 1;
        
    else
        
        out  = [];
        break;
        
    end
        
end