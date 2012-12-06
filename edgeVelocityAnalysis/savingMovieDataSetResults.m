function  savingMovieDataSetResults(ML,results)
%Saving edge velocity analysis for each cell in the MovieList


dataSetPath = [ML.outputDirectory_ filesep ML.movieListFileName_(1:end-4) 'EdgeVelocityAnalysis'];
if ~isdir(dataSetPath)
    mkdir(dataSetPath)
end

filePath = [dataSetPath filesep 'dataSetEdgeVelocity.mat'];

if exist(filePath,'file')
    copyfile(filePath,[dataSetPath filesep 'dataSetEdgeVelocity.old'])
end
save(filePath,'results');
        
