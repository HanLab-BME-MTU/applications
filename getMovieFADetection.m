function movieData = getMovieFADetection(movieData,batchMode)

%Indicate that detection was started
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

movieData.detection.filename = 'segmentParams.mat';

nFrames = numel(imageFiles);

% !!!! LINK THIS WITH THE REAL PARAMS:
% sigmaPSF = vectorialPSFSigma(1.4, 509, 67)
sigmaPSF = 1.6255;
minSize = 2;

segmentParams = cell(nFrames, 1);

for i = 1:nFrames
    % Read images
    I = imread([imagePath filesep imageFiles(i).name]);
    
    % Make sure I is type double
    I = double(I);

    % Get initial segment parameters
    segmentParams{i} = getInitialSegmentParams(I,sigmaPSF,minSize);
    
    if ~batchMode && ishandle(h)
        waitbar(i/nFrames, h)
    end
end

save([movieData.detection.directory filesep movieData.detection.filename], ...
    'segmentParams');

movieData.detection.dateTime = datestr(now);
movieData.detection.status = 1;

if ~batchMode && ishandle(h)
    close(h);
end

updateMovieData(movieData);