function status = checkMovieSegmentTracking(movieData)
%
% status = checkMovieSegmentTracking(movieData)
% Returns true if the movie specified by the input movieData has
% successfully run the tracking using getMovieSegmentTracking.m.
%
% Sylvain Berlemont, 2010

%Ensure correct OS for filenames
movieData = setupMovieData(movieData);

status = isfield(movieData,'segmentTracking') && ...
    isfield(movieData.segmentTracking,'status') && ...
    movieData.segmentTracking.status == 1 && ...
    isfield(movieData.segmentTracking,'directory') && ...
    exist([movieData.segmentTracking.directory],'dir');
