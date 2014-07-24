function  savingMovieDataSetResults(ML,results,analysisDir,fileName)
%Saving edge velocity analysis for each cell in the MovieList


dataSetPath = [ML.outputDirectory_ filesep ML.movieListFileName_(1:end-4) analysisDir];
if ~isdir(dataSetPath)
    mkdir(dataSetPath)
end

filePath = [dataSetPath filesep fileName '.mat'];

if exist(filePath,'file')
    copyfile(filePath,[dataSetPath filesep fileName '.old'])
end
save(filePath,'results');
        
