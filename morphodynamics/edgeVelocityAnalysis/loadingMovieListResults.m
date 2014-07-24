function  out = loadingMovieListResults(ML,analysisDir,fileName)

out = [];
dataSetPath = [ML.outputDirectory_ filesep ML.movieListFileName_(1:end-4) analysisDir];

filePath = [dataSetPath filesep fileName '.mat'];

if exist(filePath,'file')
    aux = load(filePath);
    out = aux.results;
end
