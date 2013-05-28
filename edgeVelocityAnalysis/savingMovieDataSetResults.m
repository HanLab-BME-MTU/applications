function  savingMovieDataSetResults(ML,results,analysis)
%Saving edge velocity analysis for each cell in the MovieList


dataSetPath = [ML.outputDirectory_ filesep ML.movieListFileName_(1:end-4) analysis];
if ~isdir(dataSetPath)
    mkdir(dataSetPath)
end

filePath = [dataSetPath filesep analysis '.mat'];

if exist(filePath,'file')
    copyfile(filePath,[dataSetPath filesep analysis '.old'])
end
save(filePath,'results');
        
