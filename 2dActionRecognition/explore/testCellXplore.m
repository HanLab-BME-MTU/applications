% script/steps to test/run cellXploreDR


load('mockData.mat');
genFakeMovies;

% playMovie(movies{1})
cellXploreDR(data, 'extra', [], 'movies', movies);

%% how to close out figures -- 
ff = findall(0,'Type', 'Figure');
close(ff);

% Example how to add new DR view choices
data.DR.rand = rand([150, 2])
