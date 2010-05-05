function status = checkMovieFigures(movieData)
%
% status = checkMovieFigures(movieData)
% Returns true if the movie specified by the input movieData has
% successfully created final figures using getMovieFigures.m.
%
% Sylvain Berlemont, 2010

%Ensure correct OS for filenames
movieData = setupMovieData(movieData);

status = isfield(movieData,'figures') && ...
    isfield(movieData.figures,'status') && ...
    movieData.figures.status == 1 && ...
    isfield(movieData.figures,'directory') && ...
    exist([movieData.figures.directory],'dir');
