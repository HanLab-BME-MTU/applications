%% load the data
load('project/movieList','ML');
MD = ML.getMovies();
% channel permutations
cpermute = [3 1 2 4];

%% make thumbnails
makeThumbs(MD{1},cpermute);
makeThumbs(MD{2});

%% plot correlation matrix
load('corrmatrix.mat');
load('corrmatrix2.mat');
corrmatrix2r = corrmatrix2(cpermute,:,cpermute,:);
plotcorrmatrix(corrmatrix2r,MD{1}.getFilename);
plotcorrmatrix(corrmatrix,MD{2}.getFilename);