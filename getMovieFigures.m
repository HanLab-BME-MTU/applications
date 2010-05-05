function movieData = getMovieFigures(movieData,batchMode)

%Indicate that figure creation was started
movieData.figures.status = 0;

% Check that tracking has been performed
assert(checkMovieTracking(movieData) == true);

movieData.figures.directory = [movieData.channels(1).analysisDirectory ...
    filesep 'figures'];

if ~exist(movieData.figures.directory, 'dir')
    mkdir(movieData.figures.directory);
end



% Load tracks
filename = [movieData.tracks.directory filesep movieData.tracks.filename];
load(filename);

% 

% Remove tracks smaller than 2
trackSEL = getTrackSEL(tracksFinal); %#ok<NODEF>
tracksFinal = tracksFinal(trackSEL(:,3) >= 2);

movieData.figures.dateTime = datestr(now);
movieData.figures.status = 1;

updateMovieData(movieData);
