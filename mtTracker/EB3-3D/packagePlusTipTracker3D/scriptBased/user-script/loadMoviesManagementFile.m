%%% Description
% Once the movie management files (MovieData and MovieList) are created. One
% can pick the movies to be loaded in a flexible manner

dataPath='/work/gdanuser/proudot/project/EB3-3D-track/packaging/alpha/plusTipTracker3D-alpha-1/data/A1_HeLa_Cells_EB1';

%% Loading all the sequence in the folder "anaphase": 
%  Load a MovieList file that links to each movie in the folder. 
anaphaseFolderPath=[dataPath filesep 'anaphase'];
anaphaseCells=MovieList.load([anaphaseFolderPath filesep 'movieList.mat']);

%% Loading specific cells in the "anaphase" folfer 
% Load the MovieData file associated to each selected movie and build a MovieList 
anaphaseFolderPath=[dataPath filesep 'anaphase'];

% Loading first Cell
anaphaseCell1Path=[anaphaseFolderPath filesep 'Cell1'];
MD1=MovieData.load([anaphaseCell1Path filesep 'analysis' filesep 'movieData.mat']);

% Loading second Cell
anaphaseCell1BisPath=[anaphaseFolderPath filesep 'Cell1bis'];
MD2=MovieData.load([anaphaseCell1BisPath filesep 'analysis' filesep 'movieData.mat']);

% Build movieList (using [MD1] or [MD1,MD2]lo)
anaphaseSelectedCells=MovieList([MD1],anaphaseFolderPath,'movieListFileName_','selectedCells.mat','movieListPath_',anaphaseFolderPath);
anaphaseSelectedCells.save();

%% Loading all the list of movies created for each condition 
% Load the file generated under <dataPath> yields all the movieList
% generated for that folder
tmp=load([dataPath filesep 'A1_HeLa_Cells_EB1.mat'])            % Load the file containing all the MovieList
arrayfun(@(x) x.sanityCheck,tmp.movieListArray,'unif',0);       % Validate each MovieList 
movieListArray=tmp.movieListArray;
