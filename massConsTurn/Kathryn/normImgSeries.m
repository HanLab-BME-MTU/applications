function normImgSeries(runInfo,nFrames)
% NORMIMGSERIES subtracts avg bg and normalizes images to 0-1
%
% DESCRIPTION: normImgSeries subtracts the average background of the movie
% and uses the max calculated from nFrames to normalize to 0-1. 
%
% SYNOPSIS: [normInfo]=normImgSeries(runInfo,nFrames)
%
% INPUT: runInfo : structure containing path to image and analysis
%                  directories.
%        nFrames : number of frames to normalize 
%                  use nFrames = 0 (default) to normalize all frames.
%
% OUTPUT: /images/norm/normMats contains .mat's of the normalized images,
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
%         (for each frame) are saved in normInfo.
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

% USER can change this: should be 1 if the max intensity value should be
% found from within the cell across the whole movie, or 0 if it should be
% found from the whole image.
insideCellMaxFlag=1;

% --- CHECK USER INPUT ---
if nargin<1
    error('normImgSeries: Not enough input parameters')
end
if ~isstruct(runInfo)
    runInfo=struct;
end

if ~isfield(runInfo,'imDir') || ~isfield(runInfo,'anDir')
    error('normImgSeries: runInfo should contain fields imDir and anDir');
else
    [anDir] = formatPath(runInfo.anDir);
    [imDir] = formatPath(runInfo.imDir);
end

% cell mask directory
cmDir=[anDir filesep 'edge' filesep 'cell_mask'];
if ~isdir(cmDir)
    error('normImgSeries: Please run edge tracking')
end

[listOfCellMasks] = searchFiles('.tif',[],cmDir,0);
nImTot=size(listOfCellMasks,1);

% check nFrames
if nargin<2 || isempty(nFrames) || ~isnumeric(nFrames) || nFrames==0
    nFrames=nImTot;
elseif nFrames>nImTot
    error('normImgSeries: nFrames must be <= the total number of images');
end

% output the normalized images to subdirectory of image directory
% delete contents if it already exists
normDir=[imDir filesep 'norm'];
if isdir(normDir)
    rmdir(normDir,'s');
end
mkdir(normDir);
mkdir([normDir filesep 'normTifs']);
mkdir([normDir filesep 'normMats']);
mkdir([normDir filesep 'negEdgeTifs']);

% create string for naming files with correct number of digits
% we want to match the naming of the image files, so use nImTot for this
s=length(num2str(nImTot));
strg=sprintf('%%.%dd',s);

% --- END CHECK USER INPUT ---

warningState=warning;
warning('off','MATLAB:intConvertNonIntVal')

% we take the global min to be the minimum gray level from the background
% pixels. (if we take the overall global min, it will always be zero
% because of the zero regions from the translation to the cell frame of
% reference.)  we take the global max to be the actual global max from the
% whole image or from the pixels inside the cell mask. we subtract
% off the background mean and scale the std accordingly.
[bgStd, bgMin, globalMax]=bgStats(imDir,cmDir,nFrames,insideCellMaxFlag);
bgStd=bgStd/(globalMax-bgMin); 

% now we normalize and save images
fgMean=zeros(nFrames,1);
edgePix=cell(1,nFrames);
[listOfImages] = searchFiles('.tif',[],imDir,0);

for i=1:nFrames
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

    edgePix{i}=find(pixelEdgeMask);

    pixelEdgeMask=swapMaskValues(pixelEdgeMask,[0 1],[1 0]);
    negEdgeImg=abs(normImg-1).*pixelEdgeMask;
    imwrite(negEdgeImg,[normDir filesep 'negEdgeTifs' filesep 'neg_edge_image',indxStr,'.tif']);
    
    cm=swapMaskValues(cm,0,nan);
    maskedIm=normImg.*cm;
    fgMean(i)=nanmean(maskedIm(:));
end

% plot the foreground mean per frame
h=get(0,'CurrentFigure');
if h==1
    close(h)
end
figure(1);
scatter(1:nFrames,fgMean,'.')
ylim([0 1])
saveas(gcf,[normDir filesep 'cellBodyMeanIntensityPerFrame.tif']);
close

save([anDir filesep 'edgePix'],'edgePix');

normInfo.bgStd=bgStd; 
normInfo.fgMean=fgMean;
save([normDir filesep 'normInfo'],'normInfo');

warning(warningState);





