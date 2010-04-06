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

movieData.detection.directory = [movieData.channels(1).analysisDirectory ...
    filesep 'detection'];

if ~exist(movieData.detection.directory, 'dir')
    mkdir(movieData.detection.directory);
end

nFrames = numel(imageFiles);

% sigmaPSF = vectorialPSFSigma(1.4, 509, 67)
sigmaPSF = 1.6255;

%Make the string for formatting
fString = strcat('%0',num2str(ceil(log10(nFrames)+1)),'.f');

for i = 1:nFrames
    I = imread([imagePath filesep imageFiles(i).name]);
    BW = imread([maskPath filesep maskFiles(i).name]);
    
    [FA, Im] = focalAdhesionDetector(I,BW,sigmaPSF); %#ok<NASGU>
    %load([movieData.detection.directory filesep 'FA_' num2str(i,fString) '.mat']);
    %load([movieData.detection.directory filesep 'Im_' num2str(i,fString) '.mat']);
    
    save([movieData.detection.directory filesep 'FA_' num2str(i,fString) '.mat'], 'FA');
    save([movieData.detection.directory filesep 'Im_' num2str(i,fString) '.mat'], 'Im');
    
    % Save image overlaid by FA
    n = size(FA,1);
    xMin = round(FA(:,1) - (FA(:,4) / 2) .* cos(FA(:,5)));
    xMax = round(FA(:,1) + (FA(:,4) / 2) .* cos(FA(:,5)));
    yMin = round(FA(:,2) - (FA(:,4) / 2) .* sin(FA(:,5)));
    yMax = round(FA(:,2) + (FA(:,4) / 2) .* sin(FA(:,5)));
    
    pts = cell2mat(arrayfun(@(i) bresenham([yMin(i) xMin(i)],[yMax(i) xMax(i)]), ...
        (1:size(FA,1))', 'UniformOutput', false));
    
    s = pts(:,1) >= 1 & pts(:,1) <= size(I,1) & ...
        pts(:,2) >= 1 & pts(:,2) <= size(I,2);
    pts = pts(s == 1,:);
    
    indPts = sub2ind(size(I), pts(:,1), pts(:,2));
    Z = zeros(size(I));
    Z(indPts) = 1;

    I = double(I);
    I = (I - min(I(:))) / (max(I(:)) - min(I(:)));
    I(indPts) = 0;    
    I = repmat(I, [1 1 3]);
    I(:,:,1) = I(:,:,1) + Z;
 
    imwrite(I, [movieData.detection.directory filesep 'overlayFA_' ...
        num2str(i,fString) '.tif'], 'Compression', 'none');
    
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