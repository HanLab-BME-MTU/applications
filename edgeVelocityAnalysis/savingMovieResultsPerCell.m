function  savingMovieResultsPerCell(ML,results)
%Saving edge velocity analysis for each cell in the MovieList

nCell = numel(ML.movies_);

for iCell = 1:nCell

    currMD = ML.movies_{iCell};
     %Saving results for each cell
    edgePath = [currMD.outputDirectory_ filesep 'EdgeVelocityAnalysis'];
    if ~isdir(edgePath)
        mkdir(edgePath)
    end
    
    edgeVelocity = results(iCell);
    filePath     = [edgePath filesep 'edgeVelocity.mat'];
    
    if exist(filePath,'file')
        copyfile(filePath,[edgePath filesep 'edgeVelocity.old'])
    end
    
    save(filePath,'edgeVelocity');
        
end

    