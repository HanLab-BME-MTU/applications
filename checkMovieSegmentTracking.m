function status = checkMovieTracking(movieData)
%
% status = checkMovieTracking(movieData)
% Returns true if the movie specified by the input movieData has
% successfully run the tracking using getMovieTracking.m.
%
% Sylvain Berlemont, 2010

%Ensure correct OS for filenames
movieData = setupMovieData(movieData);

status = isfield(movieData,'tracking') && ...
    isfield(movieData.tracking,'status') && ...
    movieData.tracking.status == 1 && ...
    isfield(movieData.tracking,'directory') && ...
    exist([movieData.tracking.directory],'dir');
