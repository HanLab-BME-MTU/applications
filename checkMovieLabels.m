function status = checkMovieLabels(movieData)
%
% status = checkMovieLabels(movieData)
% Returns true if the movie specified by the input movieData has been
% successfully labeled using getMovieLabels.m and false otherwise
%
% Sylvain Berlemont, 12/2009
%

%Ensure correct OS for filenames
movieData = setupMovieData(movieData);

status = isfield(movieData,'labels') && ...
    isfield(movieData.labels,'status') && ...
    movieData.labels.status == 1 && ...
    isfield(movieData.labels,'directory') && ...
    exist([movieData.labels.directory],'dir');
