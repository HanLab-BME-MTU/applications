function [movieInfo]=plusTipGaussianCometDetector(projData,psfSigma,varargin)
% plusTipGaussianCometDetector locates plus tip comets (or other blobs) in a movie stack
%
%SYNOPSIS [movieInfo]=plusTipCometDetector(projData,timeRange,bitDepth,savePlots)%
%INPUT  projData          : structure containing fields .anDir, which gives
%                           the full path to the roi_x directory
%                           and .imDir, which gives the full path to the
%                           folder containing the images for overlay.
%                           if given as [], program will query user for
%                           roi_x directory.
%       psfSigma          : value of the psfSigma of the point-spread function
%                           will be used in the anisotropic gaussian fit as
%                           the lowest value of the gaussian psfSigma.
%       timeRange         : row vector of the form [startFrame endFrame]
%                           indicating time range to plot. if not given or
%                           given as [], tracks from the whole movie will
%                           be displayed.
%       bitDepth          : bit depth of the images - should be 12, 14, or 16
%       savePlots         : 1 to save overlay plots of detection results,
%                           0 if not
%
%OUTPUT movieInfo         : nFrames-structure containing x/y coordinates
%
% Sebastien Besson 5/2011

warningState = warning;

%%%%%% OPTIONS TO CHANGE MANUALLY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
removeSatPixels = 0; % put one if you want to turn on this option

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CHECK INPUT AND SET UP DIRECTORIES
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('projData', @(x) isstruct(x) || isempty(x));
ip.addRequired('psfSigma', @(x) isnumeric(x));
ip.addOptional('timeRange',[],@(x) isequal(sort(size(x)),[1 2]) || isempty(x));
ip.addOptional('bitDepth',[], @(x) isnumeric(x) || isempty(x));
ip.addOptional('savePlots',1,@isscalar);
ip.addParamValue('minDist',.5, @(x) isnumeric(x));
ip.addParamValue('alpha',.01, @(x) isnumeric(x));
ip.addParamValue('displayFirstImage',0, @isscalar);
ip.parse(projData,psfSigma,varargin{:});

timeRange = ip.Results.timeRange;
bitDepth = ip.Results.bitDepth;
savePlots = ip.Results.savePlots;
minDist = ip.Results.minDist;
alpha = ip.Results.alpha;
displayFirstImage = ip.Results.displayFirstImage;

% get projData in correct format
if isempty(projData)
    % if not given as input, ask user for ROI directory
    % assume images directory is at same level
    projData.anDir=uigetdir(pwd,'Please select ROI directory');
    homeDir=pwd;
    cd(projData.anDir);
    cd('..');
    projData.imDir=[pwd filesep 'images'];
    cd(homeDir)
else
    % adjust for OS
    if ~isfield(projData,'imDir') || ~isfield(projData,'anDir')
        error('--plusTipGaussianCometDetector: first argument should be a structure with fields imDir and anDir');
    else
        [projData.anDir] = formatPath(projData.anDir);
        [projData.imDir] = formatPath(projData.imDir);
    end
end

% Get list of Images in Image Directory and count them.
listOfImages = imDir(projData.imDir);
nImTot = size(listOfImages,1);

% check timeRange input, assign start and end frame
if isempty(timeRange)
    startFrame = 1;
    endFrame = nImTot;
else
    if timeRange(1)<=timeRange(2) && timeRange(2)<=nImTot
        startFrame = timeRange(1);
        endFrame = timeRange(2);
    else
        startFrame = 1;
        endFrame = nImTot;
    end
end

nFrames = endFrame-startFrame+1;

fileNameIm = [projData.imDir filesep listOfImages(1).name];

img = double(imread(fileNameIm));
[imL,imW] = size(img);
maxIntensity = max(img(:));

% get bit depth if not given
if isempty(bitDepth)
    imgData = imfinfo(fileNameIm);
    bitDepth = imgData.BitDepth;
    disp(['bitDepth estimated to be ' num2str(bitDepth)])
end

% check bit depth to make sure it is 12, 14, or 16 and that its dynamic
% range is not greater than the provided bitDepth
if sum(bitDepth==[12 14 16])~=1 || maxIntensity > 2^bitDepth-1
    error('--plusTipCometDetector: bit depth should be 12, 14, or 16');
end

% make feat directory if it doesn't exist from batch
featDir = [projData.anDir filesep 'feat'];
if isdir(featDir)
    rmdir(featDir,'s')
end
mkdir(featDir)

overlayDir=[featDir filesep 'overlayImages'];
tifOverlayDir=[overlayDir filesep 'tifs'];

if savePlots==1
    mkdir(overlayDir);
    mkdir(tifOverlayDir);
end

% look for region of interest info from project setup step
if ~exist([projData.anDir filesep 'roiMask.tif'],'file')...
        && ~exist([projData.anDir filesep 'masks'],'dir');
    % not roi selected; use the whole image
    roiMask = ones(imL,imW);
    roiYX=[1 1; imL 1; imL imW; 1 imW; 1 1];
    multMasks = 0;
else
    %if there exists a folder called masks in the current roi dir
    if exist([projData.anDir filesep 'masks'],'dir');
        % set the mask to the mask of image 1 (first mask in list)
        roiMask = double(imread([projData.anDir filesep 'masks' filesep 'roiMask1.tif']));
        % tell the program that there is a mask for each image
        % so it will know to update it when going through the detection
        multMasks = 1;
        
    else % only one mask to load, load it here
        roiMask = double(imread([projData.anDir filesep 'roiMask.tif']));
        multMasks = 0;
    end % if exist
    
    % load roi edge pixels
    roiYX=load([projData.anDir filesep 'roiYX']);
    roiYX=roiYX.roiYX;
    
    % for some reason the donut mask will open not as 1 and 0s but as
    % 255 (an 8 bit image) and 0s, Check for this and correct so the donut mask will
    % run through detection (added MB 09/2010)
    if isempty(find(roiMask==1,1))
        roiMask(roiMask==255) = 1;
    else % keep the same
    end % isempty
    
end % if ~exist

% string for number of files
s1 = length(num2str(endFrame));
strg1 = sprintf('%%.%dd',s1);


%% START DETECTION

% initialize structure to store info for tracking
movieInfo(nImTot,1) = ...
    struct('xCoord',[],'yCoord',[],'amp',[],'sigmaX',[],'sigmaY',[],'theta',[],'bkg',[]);

% % loop thru frames and detect
count=1;
progressText(0,'Detecting comets');
if savePlots==1, saveFig = figure('Visible','off'); end
for iFrame = startFrame:endFrame
    
    progressText(count/nFrames,'Detecting comets');
    fileNameIm = [projData.imDir filesep listOfImages(iFrame).name];
    img = double(imread(fileNameIm))./((2^bitDepth)-1);
    
    if removeSatPixels == 1, img(img==1)= 0; end
    
    % if there is a mask for each image file
    if multMasks == 1 && iFrame > 1
        maskFilename = ['roiMask' num2str(iFrame) '.tif'];
        roiMask = double(imread([projData.anDir filesep 'masks' filesep maskFilename]));
    end
    
    % If mask encoded in 8 bits instead of binary
    if isempty(find(roiMask==1,1))
        roiMask(roiMask==255) = 1;
    end
    
    % Core anisotropic detection function. Returns a structure compatible
    % with Khuloud's tracker
    movieInfo(iFrame,1) = cometDetection(img,logical(roiMask),psfSigma,...
        'minDist',minDist,'alpha',alpha);
    
    indxStr1 = sprintf(strg1,iFrame); % frame
    
    %plot feat outlines and centroid on image
    if savePlots==1
        if displayFirstImage && iFrame == startFrame, 
            visibleState='on'; 
        else
            visibleState='off';
        end
        saveFig = figure('Visible',visibleState);
        imagesc(img)
        hold on
        scatter(movieInfo(iFrame,1).xCoord(:,1),movieInfo(iFrame,1).yCoord(:,1),'c.'); % plot centroid in cyan
        colormap gray
        plot(roiYX(2),roiYX(1),'w')
        axis equal
        print(saveFig,'-dtiff',[tifOverlayDir filesep 'overlay' indxStr1 '.tif']);
        saveas(saveFig,[overlayDir filesep 'overlay' indxStr1 '.fig']);
        close(saveFig)
    end
    count=count+1;


end
save([featDir filesep 'movieInfo'],'movieInfo');

warning(warningState);
