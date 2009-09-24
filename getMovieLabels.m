function movieData = getMovieLabels(movieData)

%Indicate that labeling was started
movieData.labels.status = 0;

%Verify that the windowing has been performed
if ~checkMovieWindows(movieData)
    error('Must window movie before splitting windows!!!')
end

%Load the windows
disp('Please wait, loading windows....')
load([movieData.windows.directory filesep movieData.windows.fileName]);

%Check that they loaded
if ~exist('allWinPoly','var') || isempty(allWinPoly)
    errordlg(['Problem loading windows from ' movieData.windows.directory filesep movieData.windows.fileName],mfilename);
    return
end

% Determine image dimensions according to the first image mask
if ~isfield(movieData, 'masks') || ~isfield(movieData.masks, 'directory') ||...
    ~exist(movieData.masks.directory, 'dir')
    error(['Problem loading image mask from ' movieData.analysisDirectory]);
end

maskFilenames = dir([movieData.masks.directory filesep '*.tif']);
img = imread([movieData.masks.directory filesep maskFilenames(1).name]);
[N M] = size(img);
clear img;

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

%Go through each frame and save the windows to a file
h = waitbar(0,'Please wait, labeling window frames....');

for iFrame = 1:nFrames    
    winPoly = allWinPoly(:,:,iFrame); %#ok<COLND>
    
    labels = createLabelsFromWindows(winPoly, N, M);

    imwrite(uint16(labels), [movieData.labels.directory filesep 'labels_' num2str(iFrame,fString) '.tif']);
    
    waitbar(iFrame/nFrames,h)    
end

close(h);

movieData.labels.dateTime = datestr(now);
movieData.labels.status = 1;

updateMovieData(movieData);

end