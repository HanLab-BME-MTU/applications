function status = checkMovieBWDist(movieData)
%
% status = checkMovieBWDist(movieData)
% Returns true if the movie specified by the input movieData has
% successfully compute the distance transform using getMovieLabels.m.
%
% Sylvain Berlemont, 12/2009
%

%Ensure correct OS for filenames
movieData = setupMovieData(movieData);

status = isfield(movieData,'bwdist') && ...
    isfield(movieData.bwdist,'status') && ...
    movieData.bwdist.status == 1 && ...
    isfield(movieData.bwdist,'directory') && ...
    exist([movieData.bwdist.directory],'dir');
