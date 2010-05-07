function status = checkMovieDetection(movieData)
%
% status = checkMovieDetection(movieData)
% Returns true if the movie specified by the input movieData has
% successfully run the detection using getMovieDetection.m.
%
% Sylvain Berlemont, 2010

%Ensure correct OS for filenames
movieData = setupMovieData(movieData);

status = isfield(movieData,'detection') && ...
    isfield(movieData.detection,'status') && ...
    movieData.detection.status == 1 && ...
    isfield(movieData.detection,'directory') && ...
    exist([movieData.detection.directory],'dir');
