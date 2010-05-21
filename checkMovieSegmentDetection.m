function status = checkMovieSegmentDetection(movieData)
%
% status = checkMovieSegmentDetection(movieData)
% Returns true if the movie specified by the input movieData has
% successfully run the detection using getMovieSegmentDetection.m.
%
% Sylvain Berlemont, 2010

%Ensure correct OS for filenames
movieData = setupMovieData(movieData);

status = isfield(movieData,'segmentDetection') && ...
    isfield(movieData.segmentDetection,'status') && ...
    movieData.segmentDetection.status == 1 && ...
    isfield(movieData.segmentDetection,'directory') && ...
    exist([movieData.segmentDetection.directory],'dir');
