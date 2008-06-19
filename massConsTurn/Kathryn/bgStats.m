function [bgStd, bgMin, globalMax]=bgStats(imDir,cmDir,nFrames,nFms2avg)
%BGSTATS gives background mean and std for image series
%
% DESCRIPTION: bgStats finds the intensity max and background min,std for a
% non-normalized image series.  it finds the background region by ignoring
% pixels outside a dilated cell mask (ie, in case image segmentation wasn't
% very good, we expand the mask to make sure the values used are from the
% background). also, bgStats avoids regions where image intensity is 0 (eg,
% in keratocyte cell-frame-of-reference images, the rotation procedure
% creates artificial boundaries).
%
% bgStats accepts .tif or .mat images and masks, though if .mat, images
% need to load to 'normImg' and masks need to load to 'mask'.  be aware
% that because of the dilation step around 0-intensity regions, this
% function SHOULD NOT BE USED ON NORMALIZED IMAGES where the background
% has been subtracted.
%
% SYNOPSIS: [bgStd, bgMin, globalMax]=bgStats(imDir,cmDir,nFrames,nFms2avg)
%
% INPUT: imDir     : path to image directory
%        cmDir     : path to cell mask directory
%        nFrames   : number of frames over which to find stats
%                    use 0 (default) to find stats for each frame in imDir
%        nFms2avg  : the number of frames used to calculate local average
%                    (e.g. if nFms2avg=5, bgStats will calculate the avg bg
%                    from frames i-2, i-1, i, i+1, i+2)
%                    use 0 (default) to return stats calculated over
%                    nFrames THIS PARAMETER IS NO LONGER FUNCTIONAL...
%
% OUTPUT: bgStd    : standard deviation from background pixels
%         bgMin    : minimum intensity from background pixels
%         globalMax: max intensity from any pixel in any frame
%         
%         cmDir/bgMasks: a directory is created with masks of the bg
%                        regions used for each frame
%
% MATLAB VERSION (originally written on): 7.2.0.232 (R2006a) Windows_NT
%
% USERNAME: kathomps
% DATE: 13-Jun-2007
%
%



% assume images are tifs, but check for .mat
matformatim=0;
matformatcm=0;

% need at least imDir and cmDir
if nargin<2
    error('BGSTATS: not enough input arguments');
end

% check for valid input of imDir
if ~isstr(imDir)
    error('BGSTATS: imDir must be a string');
elseif ~isdir(imDir)
    error('BGSTATS: imDir does not exist')
else % if valid directory, check for tifs
    [listOfImages] = searchFiles('.tif',[],imDir,0);
    if isempty(listOfImages) % if no tifs, check for .mat files
        [listOfImages] = searchFiles('.mat',[],imDir,0);
        if isempty(listOfImages) % if none, report error
            error('BGSTATS: imDir is empty')
        end
        fileNameIm=[char(listOfImages(1,2)) filesep char(listOfImages(1,1))];
        im=load(fileNameIm);
        if ~isfield(im,'normImg') % if .mat files, make sure they load to normImg
            error('BGSTATS: .mat images must load to normImg')
        end
        matformatim=1;
    end
end

% check for valid input of cmDir
if ~isstr(cmDir)
    error('BGSTATS: cmDir must be a string');
elseif ~isdir(cmDir)
    error('BGSTATS: cmDir does not exist')
else % if valid directory, check for tifs
    [listOfCellMasks] = searchFiles('.tif',[],cmDir,0);
    if isempty(listOfCellMasks) % if no tifs, check for .mat files
        [listOfCellMasks] = searchFiles('.mat',[],cmDir,0);
        if isempty(listOfCellMasks) % if none, report error
            error('BGSTATS: cmDir is empty')
        end
        fileNameMask=[char(listOfCellMasks(1,2)) filesep char(listOfCellMasks(1,1))];
        cm=load(fileNameMask);
        if ~isfield(cm,'mask') % if .mat files, make sure they load to normImg
            error('BGSTATS: .mat masks must load to mask')
        end
        matformatcm=1;
    end
end

% output the bg masks to subdirectory of cmDir
bgMaskDir=[cmDir filesep 'bgMasks'];
if ~isdir(bgMaskDir)
    mkdir(bgMaskDir);
else
    delete([bgMaskDir,filesep,'bgMask' filesep '*tif'])
end

% default for nFrames = 0 (all the frames)
if nargin<3 || isempty(nFrames)
    nFrames=0;
else
    if nFrames<0 || (nFrames>0 && mod(nFrames,1)~=0)
        % error if negative or positive but not an integer
        error('BGSTATS: nFrames must be a whole number');
    end
end

% default for nFms2avg = 0 (get stats over nFrames)
if nargin<4 || isempty(nFms2avg)
    nFms2avg=0;
else
    if nFms2avg<0 || (nFms2avg>0 && (mod(nFms2avg,1)~=0 || mod(nFms2avg,2)==0))
        % error if not zero and either negative, not an integer, or not odd
        error('BGSTATS: nFms2avg must be an odd whole number');
    end
end

% can't average over more the number of frames
if nFms2avg>nFrames
    error('BGSTATS: nFms2avg must be less than or equal to nFrames');
end

% find number of images and masks
nImTot=size(listOfImages,1);
nCmTot=size(listOfCellMasks,1);

s=length(num2str(nImTot));
strg=sprintf('%%.%dd',s);

% if nFrames = 0 (default), use the total number of frames in the directory
if nFrames==0
    nFrames=nImTot;
end

% in case the user wants to calculate stats on fewer than the total number
% of frames, or if some images don't have masks (or vice versa), we use the
% minimum number
nFrames2Use=min([nImTot nCmTot nFrames]);

% initialize vector to store data from all frames
fileNameIm=[char(listOfImages(1,2)) filesep char(listOfImages(1,1))];
if matformatim==0 % tiff
    im=double(imread(fileNameIm)); %use first image to get size
else % normImg.mat when loaded
    im=load(fileNameIm);
    im=im.normImg;
end
[imL imW]=size(im);
bgPixVal=nan*zeros(imL,imW,nFrames2Use);

% loop thru frames and put bg masks into bgPixVal
globalMax=0;
for i=1:nFrames2Use

    % read or load image
    fileNameIm=[char(listOfImages(i,2)) filesep char(listOfImages(i,1))];
    if matformatim==0 % tiff
        im=double(imread(fileNameIm)); %im is raw image
    else % normImg.mat when loaded
        im=load(fileNameIm);
        im=im.normImg;
    end

    % keep track of global max
    globalMax=max(globalMax,nanmax(im(:))); 
    
    % read or load mask
    fileNameMask=[char(listOfCellMasks(i,2)) filesep char(listOfCellMasks(i,1))];
    if matformatcm==0 % tiff
        cm=double(imread(fileNameMask)); %cm is cell mask
    else
        cm=load(fileNameMask);
        cm=double(cm.mask);
    end

    % dilate the cell mask in case segmentation wasn't perfect
    cmDil=bwmorph(cm,'dilate',40);

    % bgm is 1 wherever the raw image is 0
    bgm=zeros(size(im));
    bgm(im==0)=1;

    % in the keratocyte cell frame of reference, there is often a zone of
    % 0's from the alignment procedure.  we dilate this zone to remove edge
    % effects
    bgmDil=bwmorph(bgm,'dilate',10);
    bgMask=zeros(size(im));
    bgMask(bgmDil | cmDil)=1; % this is inverse of what we want

    % make it NaN in cell body and in former 0 zones, and 1 across bg
    % in other words, 1=bg, NaN=everywhere else
    bgMask=swapMaskValues(bgMask,[0 1],[1 nan]);

    % save BW mask in cmDir/bgMask as tiff: 1=bg, 0=everywhere else
    bgMaskBW=swapMaskValues(bgMask,[nan],[0]);
    indxStr=sprintf(strg,i);
    imwrite(logical(bgMaskBW),[bgMaskDir,filesep,'bgMask',indxStr,'.tif']);

    % create image with NaNs everywhere except background
    bgPixVal(:,:,i)=bgMask.*im;

    % this commented-out section was used to verify that the std/sqrt(n) is
    % essentially the same if you actually sample an image versus finding
    % the std over the whole movie.  this is an uninteresting result, but i
    % nevertheless wanted to see it for myself.
%     % for the first frame, go through the possible window sizes and get
%     % measures of the intensity std over the bg
%     if i==1
%         winRange=[15:-2:5]; nWinSizes=length(winRange);
%         c1=1; % counter for window sizes
%         for winL=winRange
%             [tiledArray]=tileSquaresWithIndex(im,winL); % tile first image (bg only) with boxes (side length = winL)
%             [arrayL,arrayW]=size(tiledArray);
%             imCropped=bgPixVal(1:arrayL,1:arrayW,1); % crop the image to same size as tiledArray
% 
%             nBoxes=max(tiledArray(:)); % number of boxes that fit in the tiled array
%             % initialize matrix to contain std/sqrt(N) values based on largest box size
%             if winL==winRange(1)
%                 perPixVar=zeros(nBoxes,nWinSizes);
%                 nValues=nBoxes;
%             end
%             c2=1; %counter for squares that have non-NaN values
%             for square=1:nBoxes % loop through boxes in image to get std/sqrt(N)                
%                 intensities=imCropped(tiledArray==square);
%                 varVal=std(intensities)/winL;
%                 % discard if NaN - occurs when part/all of box is outside
%                 % bg of image
%                 if ~isnan(varVal)
%                     perPixVar(c2,end-c1+1)=varVal; % record small-to-large
%                     c2=c2+1;
%                     if c2>nValues
%                         break
%                     end
%                     
%                 end
%             end
%             if winL==winRange(1)
%                 perPixVar(perPixVar(:,end)==0,:)=[];
%                 nValues=size(perPixVar,1);
%             end
%             c1=c1+1;
%         end
%     end
% 
end

% here we get the bg mean and std over the whole movie
bgMean=round(nanmean(bgPixVal(:)));
bgStd=round(nanstd(bgPixVal(:)));
bgMin=nanmin(bgPixVal(:));


% this commented section is vestigial from when i was calculating the mean
% and std of the bg pixels using a moving window around the frame. it
% essentially is meaningless for
% if nFms2avg==0 % average over all nFrames
%     bgMean=nanmean(bgPixVal(:));
%     bgStd=nanstd(bgPixVal(:));
% else % do calculation for each frame
%     for i=1:nFrames2Use
%         sF=i-(nFms2avg-1)/2; %starting frame
%         eF=i+(nFms2avg-1)/2; %ending frame
%
%         if sF<1
%             sF=1;
%         end
%         if eF>nFrames2Use
%             eF=nFrames2Use;
%         end
%
%         usefulPart=bgPixVal(:,:,sF:eF);
%         bgMean(i)=nanmean(usefulPart(:));
%         bgStd(i)=nanstd(usefulPart(:));
%     end
% end
