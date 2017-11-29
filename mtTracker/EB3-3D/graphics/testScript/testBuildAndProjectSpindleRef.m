%% Building test dataset
% testing on metaphase data
allMovieToAnalyse=readtable('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/phaseProgression/analysis/movieTables/allMovieToAnalyse.xlsx');
allMovieToAnalyse=allMovieToAnalyse(~(allMovieToAnalyse.blurred|allMovieToAnalyse.doubleCell),:);
outputPath='/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/proudot/u-track-3D/dynROI/buildRef';

movieIndex=find(ismember(allMovieToAnalyse.Setup_min_,'16'))';
ML=MovieList(allMovieToAnalyse.rawPath(movieIndex(3)),outputPath,'movieListFileName_','metaphaseTest');
ML=ML.addAnalysisFolder('/project/bioinformatics/Danuser_lab/externBetzig/raw/adavid/lattice/',outputPath)
%%
for mIdx=1:ML.getSize()
    buildAndProjectSpindleRef(ML.getMovie(mIdx));
end