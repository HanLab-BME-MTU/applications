function movieData = getMovieFADetection(movieData,batchMode)

%Indicate that labeling was started
movieData.detection.status = 0;

if ~batchMode
    h = waitbar(0,'Please wait, focal adhesion detection....');
end

imagePath = movieData.channels(1).roiDirectory;
imageFiles = dir([imagePath filesep '*.tif']);

movieData.detection.directory = [movieData.channels(1).analysisDirectory ...
    filesep 'detection'];

if ~exist(movieData.detection.directory, 'dir')
    mkdir(movieData.detection.directory);
end

nFrames = numel(imageFiles);

% !!!! LINK THIS WITH THE REAL PARAMS:
% sigmaPSF = vectorialPSFSigma(1.4, 509, 67)
sigmaPSF = 1.6255;
minSize = 2;

%Make the string for formatting
fString = strcat('%0',num2str(ceil(log10(nFrames)+1)),'.f');

for i = 1:nFrames
    % Read images
    I = imread([imagePath filesep imageFiles(i).name]);
    
    % Make sure I is type double
    I = double(I);

    % Get initial segment parameters
    FA = getInitialSegmentParams(I,sigmaPSF,minSize); %#ok<NASGU>

    save([movieData.detection.directory filesep 'FA_' num2str(i,fString) '.mat'], 'FA');
    
    if ~batchMode && ishandle(h)
        waitbar(i/nFrames, h)
    end
end

movieData.detection.dateTime = datestr(now);
movieData.detection.status = 1;

if ~batchMode && ishandle(h)
    close(h);
end

updateMovieData(movieData);