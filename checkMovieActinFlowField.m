function status = checkMovieActinFlowField(movieData)
%
% status = checkMovieActinFlowField(movieData)
% Returns true if the movie specified by the input movieData has
% successfully run fsmCenter post process flow field
%
% Sylvain Berlemont, 2011

%Ensure correct OS for filenames
movieData = setupMovieData(movieData);

status = isfield(movieData,'actinFlowField') && ...
    isfield(movieData.actinFlowField,'status') && ...
    movieData.actinFlowField.status == 1 && ...
    isfield(movieData.actinFlowField,'directory') && ...
    exist([movieData.actinFlowField.directory],'dir');
