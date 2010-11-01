function status = checkMovieParticleTracking(movieData)
%
% status = checkMovieParticleTracking(movieData)
% Returns true if the movie specified by the input movieData has
% successfully run the tracking using getMovieParticleTracking.m.
%
% Sylvain Berlemont, 2010

%Ensure correct OS for filenames
movieData = setupMovieData(movieData);

status = isfield(movieData,'particleTracking') && ...
    isfield(movieData.particleTracking,'status') && ...
    movieData.particleTracking.status == 1 && ...
    isfield(movieData.particleTracking,'directory') && ...
    exist([movieData.particleTracking.directory],'dir');
