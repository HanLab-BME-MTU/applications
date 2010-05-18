function movieData = getMovieLabels(movieData, batchMode)

%Indicate that labeling was started
movieData.labels.status = 0;

imSize = movieData.imSize';

%Verify that the windowing has been performed
if ~checkMovieWindows(movieData)
    error('Must window movie before labeling windows.');
end

if ~checkMovieProtrusionSamples(movieData)
    error('Must sample protrusion before labeling windows.');
end

% Step 1: Classify protrusion speed into 3 classes:

% Add this classification into protrusionSamples variable. This should not
% be done here. Instead, add this lines into getProtrusionSamples.m (Talk
% to Hunter).

%Load protrusion sample
load(fullfile(movieData.protrusion.directory, ...
    movieData.protrusion.samples.fileName));

%Check that the protrusion sample map loaded
if ~exist('protrusionSamples','var')
    error(['Problem loading protrusionSamples from ' movieData.windows.directory ...
        filesep movieData.windows.fileName],mfilename);
end

if ~isfield(protrusionSamples,'classNames') || ~isfield(protrusionSamples,'classes') %#ok<NODEF>
    % Store class names
    protrusionSamples.classNames = {'Pause', 'Protrusion', 'Retraction'};
    
    % classify speeds
    protrusionSamples.classes = zeros(size(protrusionSamples.averageNormalComponent));
    protrusionSamples.classes(protrusionSamples.averageNormalComponent >= 1) = 1;
    protrusionSamples.classes(protrusionSamples.averageNormalComponent <= -1) = 2;

    % Update protrusionSamples file
    save(fullfile(movieData.protrusion.directory, ...
        movieData.protrusion.samples.fileName), 'protrusionSamples');
end

% Step 2: Go through each frame and save the windows to a file

%Load the windows
load([movieData.windows.directory filesep movieData.windows.fileName]);

%Check that they loaded
if ~exist('allWinPoly','var') || isempty(allWinPoly) %#ok<NODEF>
    error(['Problem loading windows from ' movieData.windows.directory ...
        filesep movieData.windows.fileName],mfilename);
end

movieData.labels.directory = [movieData.analysisDirectory filesep 'labels'];

if ~exist(movieData.labels.directory, 'dir')
    mkdir(movieData.labels.directory);
end

%Determine number of windows/bands
[nBands,nWindows,nFrames] = size(allWinPoly);

movieData.labels.nWindows = nWindows;
movieData.labels.nBands = nBands;
movieData.labels.nFrames = nFrames;

%Make the string for formatting
fString = strcat('%0',num2str(ceil(log10(nFrames)+1)),'.f');

if ~batchMode
    h = waitbar(0,'Please wait, labeling window frames...');
end

for iFrame = 1:nFrames-1
    winPoly = allWinPoly(:,:,iFrame);
    
    labels = createLabelsFromWindows(winPoly, imSize);

    imwrite(uint16(labels), [movieData.labels.directory filesep 'labels_' ...
        num2str(iFrame,fString) '.tif'], 'Compression', 'lzw');
    
    if ~batchMode && ishandle(h)
        waitbar(iFrame/nFrames,h)
    end
end

if ~batchMode && ishandle(h)
    close(h);
end

movieData.labels.dateTime = datestr(now);
movieData.labels.status = 1;

updateMovieData(movieData);

end