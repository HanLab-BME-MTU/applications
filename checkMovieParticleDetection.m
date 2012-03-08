function status = checkMovieParticleDetection(movieData)
%
% status = checkMovieParticleDetection(movieData)
% Returns true if the movie specified by the input movieData has
% successfully run the detection using getMovieParticleDetection.m.
%
% Sylvain Berlemont, 2010

%Ensure correct OS for filenames
movieData = setupMovieData(movieData);

status = isfield(movieData,'particleDetection') && ...
    isfield(movieData.particleDetection,'status') && ...
    movieData.particleDetection.status == 1 && ...
    isfield(movieData.particleDetection,'directory') && ...
    exist([movieData.particleDetection.directory],'dir');
