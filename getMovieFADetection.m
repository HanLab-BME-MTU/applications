function movieData = getMovieSegmentation(movieData, batchMode)

%Indicate that labeling was started
movieData.segmentation.status = 0;

%Go through each frame and save the windows to a file
if ~batchMode
    h = waitbar(0,'Please wait, image segmentation....');
end

imagePath = movieData.channels(1).roiDirectory;
imageFiles = dir([imagePath filesep '*.tif']);

segmentationPath = [movieData.channels(1).analysisDirectory filesep ...
    'segmentation'];

if ~exist(segmentationPath, 'dir')
    mkdir(movieData.channels(1).analysisDirectory, 'segmentation');
end

nFrames = numel(imageFiles);

for i = 1:nFrames
    I = imread([imagePath filesep imageFiles(i).name]);

    BW = blobSegmentThreshold(I, 1, 0);
    
    imwrite(logical(BW), [segmentationPath filesep imageFiles(i).name], ...
        'Compression', 'none');
    
    if ~batchMode && ishandle(h)
        waitbar(i/nFrames, h)
    end
end

movieData.segmentation.dateTime = datestr(now);
movieData.segmentation.status = 1;

if ~batchMode && ishandle(h)
    close(h);
end

updateMovieData(movieData);