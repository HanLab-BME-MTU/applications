function [runInfo]=normImgSeries(runInfo,nFrames)
% NORMIMGSERIES subtracts avg bg and normalizes images to 0-1
%
% DESCRIPTION: normImgSeries subtracts the average background of the movie
% and uses the max calculated from nFrames to normalize to 0-1. 
%
% SYNOPSIS: normImgSeries(runInfo,nFrames)
%
% INPUT: runInfo : structure containing path to image and analysis
%                  directories. use runInfo=[] to query user for
%                  directories and normalize all frames
%        nFrames : number of frames to normalize 
%                  use nFrames = 0 (default) to normalize all frames.
%
% OUTPUT: /analysis/norm/normMats contains .mat's of the normalized images,
%                                 which are between 0 and 1.
%                  /norm/normTifs contains .tiff's of the normalized images
%                  /norm/negEdgeTifs contains .tiffs of the image negative
%                                 with the cell edge shown in black.  This
%                                 can be used to see how well the cell mask
%                                 (used to create edge) segments the image.
%                  /edge/cell_mask/bgMasks contains tiff images showing the
%                                 regions used to calculate the mean and
%                                 std (using bgStats, called here) in
%                                 white, and the excluded regions in black.
%         The scaled bg standard deviation (for whole movie) and fg mean
%         (for each frame) are saved in runInfo.
%
% This function will delete contents of norm from a previous run.
%
% MATLAB VERSION (originally written on): 7.2.0.232 (R2006a) Windows_NT
%
%
%
% USERNAME: kathomps
% DATE: 14-Jun-2007
%
%

warningState=warning;
warning('off','MATLAB:intConvertNonIntVal')

% this is the number of frames used to calculate the local average
% intensity to subtract - should be odd.  see bgStats for details.
frameRange=5; % this is defunct now, don't worry about changing it

[anDir,imDir,cmDir]=checkInputParameters(runInfo,nFrames);

% get list of images and total number of images
[listOfImages] = searchFiles('.tif',[],imDir);
nImTot=size(listOfImages,1);

% give error if nFrames is larger than the total number of images
if nFrames>nImTot
    disp('normImgSeries: nFrames must be <= the total number of images');
    disp('               nFrames will now be the total number of images');
    nFrames=nImTot;
end

% create string for naming files with correct number of digits
% we want to match the naming of the image files, so use nImTot for this
s=length(num2str(nImTot));
strg=sprintf('%%.%dd',s);

% output the normalized images to subdirectory of project directory
% delete contents if it already exists
normDir=[anDir filesep 'norm'];
if ~isdir(normDir)
    mkdir(normDir);
end
if ~isdir([normDir filesep 'normTifs'])
    mkdir([normDir filesep 'normTifs']);
else
    delete([normDir filesep 'normTifs' filesep '*tif'])
end
if ~isdir([normDir filesep 'normMats'])
    mkdir([normDir filesep 'normMats']);
else
    delete([normDir filesep 'normMats' filesep '*mat'])
end
if ~isdir([normDir filesep 'negEdgeTifs'])
    mkdir([normDir filesep 'negEdgeTifs']);
else
    delete([normDir filesep 'negEdgeTifs' filesep '*tif'])
end
if isdir([cmDir filesep 'bgMasks']) % mkdir done by bgStats.m
   delete([cmDir filesep 'bgMasks' filesep '*tif'])
end

% nImTot is given the user-set nFrames if not 0
if nFrames~=0
    nImTot=nFrames; 
end

% we take the global min to be the minimum gray level from the background
% pixels. (if we take the overall global min, it will always be zero
% because of the zero regions from the translation to the cell frame of
% reference.)  we take the global max to be the actual global max from the
% whole image because we believe all these values are accurate. we subtract
% off the background mean and scale the std accordingly.
[bgStd, bgMin, globalMax]=bgStats(imDir,cmDir,nImTot,frameRange);
bgStd=bgStd/(globalMax-bgMin); 

[listOfCellMasks] = searchFiles('.tif',[],cmDir,0);

% now we normalize and save images
fgMean=zeros(nImTot,1);
for i=1:nImTot
    fileNameIm=[char(listOfImages(i,2)) filesep char(listOfImages(i,1))];
    im=double(imread(fileNameIm)); %im is raw image
    normImg=(im-bgMin)./(globalMax-bgMin); % normalize using global max/min
    
    % cut off values under 0 or over 1
    normImg(normImg<0)=0;
    normImg(normImg>1)=1; 
    
    % save normalized image as .mat and .tif
    indxStr=sprintf(strg,i);
    save([normDir filesep 'normMats' filesep 'norm_image' indxStr],'normImg');
    imwrite(normImg,[normDir filesep 'normTifs' filesep 'norm_image',indxStr,'.tif']);

    % create image negative and show cell edge from cell mask
    fileNameMask=[char(listOfCellMasks(i,2)) filesep char(listOfCellMasks(i,1))];
    cm=double(imread(fileNameMask)); %cm is cell mask
    pixelEdgeMask=bwmorph(cm,'remove');
    pixelEdgeMask=swapMaskValues(pixelEdgeMask);
    negEdgeImg=abs(normImg-1).*pixelEdgeMask;
    imwrite(negEdgeImg,[normDir filesep 'negEdgeTifs' filesep 'neg_edge_image',indxStr,'.tif']);
    
    cm=swapMaskValues(cm,0,nan);
    maskedIm=normImg.*cm;
    fgMean(i)=nanmean(maskedIm(:));
end

% plot the foreground mean per frame
scatter(1:nImTot,fgMean,'.')
ylim([0 1])
saveas(gcf,[normDir filesep 'cellBodyMeanIntensityPerFrame.fig']);
save([normDir filesep 'bgStd'],'bgStd');
save([normDir filesep 'fgMean'],'fgMean');

runInfo.bgStd=bgStd; 
runInfo.fgMean=fgMean;

warning(warningState);

% --------------------------
% subfunction to check input
function [anDir,imDir,cmDir]=checkInputParameters(runInfo,nFrames)

% check nFrames
if nargin<2 || isempty(nFrames) || ~isnumeric(nFrames) 
        error('normImgSeries: nFrames must be a positive integer, or use 0 to normalize over all frames');
else
    if nFrames<0 || (nFrames>0 && mod(nFrames,1)~=0)
        % error if not zero and either negative or not an integer
        error('normImgSeries: nFrames must be a positive integer, or use 0 to normalize over all frames');
    end
end

% if no arguments given, query user for directories and use nFrames=0
if nargin<1 || ~isstruct(runInfo)
    anDir=uigetdir(pwd ,'Please select project analysis directory');
    imDir=uigetdir(anDir,'Choose directory containing tifs to be normalized');
end

% check that runInfo contains imDir and anDir and that they are valid paths
if isstruct(runInfo)
    if isfield(runInfo,'imDir') && isfield(runInfo,'anDir') && isdir(runInfo.anDir) && isdir(runInfo.imDir)
        anDir=runInfo.anDir;
        imDir=runInfo.imDir;
    else
        error('normImgSeries: runInfo must be a structure containing imDir and anDir paths, or you may use [] to query user');
    end
end

cmDir=[anDir filesep 'edge' filesep 'cell_mask'];
if ~isdir(cmDir)
    error('normImgSeries: Please run edge tracking')
end