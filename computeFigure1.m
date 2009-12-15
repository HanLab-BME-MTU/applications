function movieData = computeFigure1(movieData, batchMode)

%Indicate that Figure 1 computing was started
movieData.output.fig1.status = 0;

%Verify that the labeling has been performed
if ~checkMovieLabels(movieData)
    error('Must label movie before computing figure 1.');
end

movieData.output.fig1.dateTime = datestr(now);
movieData.output.fig1.status = 1;

% Load labels
L = cell(currMovie.labels.nFrames, 1);
filenames = dir([currMovie.labels.directory filesep '*.tif']);
for iFrame = 1:currMovie.labels.nFrames
    L{iFrame} = imread([currMovie.labels.directory filesep filenames(iFrame).name]);
end

% Load speckles channel 1
% Load speckles channel 2