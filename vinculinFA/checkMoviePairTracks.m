function status = checkMoviePairTracks(movieData)
%
% status = checkMoviePairTracks(movieData)
% Returns true if the movie specified by the input movieData has
% successfully run the pairTracks using getMoviePairTracks.m.
%
% Sylvain Berlemont, 2010

%Ensure correct OS for filenames
movieData = setupMovieData(movieData);

status = isfield(movieData,'pairTracks') && ...
    isfield(movieData.pairTracks,'status') && ...
    movieData.pairTracks.status == 1 && ...
    isfield(movieData.pairTracks,'directory') && ...
    exist([movieData.pairTracks.directory],'dir');
