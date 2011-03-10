function status = checkMovieDistanceTransform(movieData)
%
% status = checkMovieDistanceTransform(movieData)
% Returns true if the movie specified by the input movieData has
% successfully compute the distance transform using
% getMovieDistanceTransform.m
%
% Sylvain Berlemont, 12/2009
%

%Ensure correct OS for filenames
movieData = setupMovieData(movieData);

status = isfield(movieData,'distanceTransform') && ...
    isfield(movieData.distanceTransform,'status') && ...
    movieData.distanceTransform.status == 1 && ...
    isfield(movieData.distanceTransform,'directory') && ...
    exist([movieData.distanceTransform.directory],'dir');
