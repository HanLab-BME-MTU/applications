function movieData = getMovieFADetection(movieData, batchMode)

%Indicate that labeling was started
movieData.detection.status = 0;

% Go through each frame and save the windows to a file
if ~batchMode
    h = waitbar(0,'Please wait, image detection....');
end

imagePath = movieData.channels(1).roiDirectory;
imageFiles = dir([imagePath filesep '*.tif']);

maskPath = movieData.masks.directory;
maskFiles = dir([maskPath filesep '*.tif']);

detectionPath = [movieData.channels(1).analysisDirectory filesep ...
    'detection'];

if ~exist(detectionPath, 'dir')
    mkdir(movieData.channels(1).analysisDirectory, 'detection');
end

nFrames = numel(imageFiles);

% sigmaPSF = vectorialPSFSigma(1.4, 509, 67)
sigmaPSF = 1.6255;

for i = 1:nFrames
    I = imread([imagePath filesep imageFiles(i).name]);
    BW = imread([maskPah filesep maskFiles(i).name]);
    
    [FA, Im] = focalAdhesionDetector(I,BW,sigmaPSF);
    
    save();
    
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