function movieData = getMovieFigures(movieData,batchMode)

assert(checkMoviePairTracks(movieData));

movieData.figures.status = 0;

movieData.figures.directory = fullfile(movieData.analysisDirectory, 'figures');

% Create output directory
if ~exist(movieData.figures.directory, 'dir')
    mkdir(movieData.figures.directory);
end

nFrames = movieData.nImages(1);

%% Load segments

load(fullfile(movieData.pairTracks.directory, ['CCParams_iter=' ...
    num2str(movieData.pairTracks.params.maxIter-1) '_.mat']));

%% END
movieData.figures.status = 1;
