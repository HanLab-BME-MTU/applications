function status = checkMovieDensity(movieData)
%
% status = checkMovieDensity(movieData)
% Returns true if the movie specified by the input movieData has been
% successfully compute density maps using getMovieDensity.m and false
% otherwise.
%
% Sylvain Berlemont, 12/2009
%

%Ensure correct OS for filenames
movieData = setupMovieData(movieData);

status = isfield(movieData,'density') && ...
    isfield(movieData.density,'status') && ...
    movieData.density.status == 1 && ...
    isfield(movieData.density,'directory') && ...
    exist([movieData.density.directory],'dir');
