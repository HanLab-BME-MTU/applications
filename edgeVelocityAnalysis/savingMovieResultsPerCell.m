function  savingMovieResultsPerCell(ML,results,analysis,fileName)
%Saving edge velocity analysis for each cell in the MovieList
%
%USAGE:
%       savingMovieResultsPerCell(ML,results,analysis,fileName)
%Input:
%       ML       - movieList object
%       results  - strucutre array resulting from analysis: EdgeVelocityQuantification, sampledSignalQuantification or AssociationEstimation
%       analysis - string with analysis type: edgeVelocity, sampledSignal or Associations
%       fileName - guess!
%
%Marco Vilela, 2013

nCell = numel(ML.movies_);

for iCell = 1:nCell

    currMD = ML.movies_{iCell};
     %Saving results for each cell
    edgePath = [currMD.outputDirectory_ filesep analysis];
    if ~isdir(edgePath)
        mkdir(edgePath)
    end
    
    analysisResults = results{iCell};
    filePath        = [edgePath filesep fileName '.mat'];
    
    if exist(filePath,'file')
        copyfile(filePath,[edgePath filesep fileName '.old'])
    end
    
    save(filePath,'analysisResults');
        
end

    