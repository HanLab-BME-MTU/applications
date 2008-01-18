function [bgMean, bgStd]=bgStats(imDir,cmDir,nFrames,nFms2avg)
%BGSTATS finds the mean and std of the background for an image series
%
% DESCRIPTION: bgStats finds the background mean and std for a
% non-normalized image series.  it finds the background region by ignoring
% pixels outside a dilated cell mask (ie, in case image segmentation wasn't
% very good, we expand the mask to make sure the values used are from the
% background). also, bgStats avoids regions where image intensity is 0 (eg,
% in keratocyte cell-frame-of-reference images, the rotation procedure
% creates artificial boundaries).  please note that if nFms2avg > 1, the
% stats for the first and last few frames will be based off of fewer
% images.  eg) if nFrames=10 and nFms2avg=5, then bgMean for frame 1 will
% be calcuated from frames 1-3 (3 frames), for frame 2 from 1-4 (4), 3 from
% 1-5 (5), 4 from 2-6 (5),..., 9 from 7-10 (4), 10 from 8-10 (3).
%
% bgStats accepts .tif or .mat images and masks, though if .mat, images
% need to load to 'normImg' and masks need to load to 'mask'.  be aware
% that because of the dilation step around 0-intensity regions, this
% function SHOULD NOT BE USED ON NORMALIZED IMAGES where the background
% has been subtracted.
%
% SYNOPSIS: [bgMean, bgStd]=bgStats(imDir,cmDir,nFrames,nFms2avg)
%
% INPUT: imDir     : path to image directory
%        cmDir     : path to cell mask directory
%        nFrames   : number of frames over which to find stats
%                    use 0 (default) to find stats for each frame in imDir
%        nFms2avg  : the number of frames used to calculate local average
%                    (e.g. if nFms2avg=5, bgStats will calculate the avg bg
%                    from frames i-2, i-1, i, i+1, i+2)
%                    use 0 (default) to return stats calculated over
%                    nFrames
%
% OUTPUT: bgMean   : average background intensity
%         bgStd    : standard deviation
%         these are calculated for every image unless nFms2avg=0 (see
%         input)
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

% initialize vectors to store stats for frames
if nFms2avg==0 % will average over whole series
    bgMean=0;
    bgStd=0;
else % will find value for each frame by averaging over nFms2avg
    bgMean=zeros(nFrames2Use,1);
    bgStd=zeros(nFrames2Use,1);
end

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
for i=1:nFrames2Use

    % read or load image
    fileNameIm=[char(listOfImages(i,2)) filesep char(listOfImages(i,1))];
    if matformatim==0 % tiff
        im=double(imread(fileNameIm)); %im is raw image
    else % normImg.mat when loaded
        im=load(fileNameIm);
        im=im.normImg;
    end

    % read or load mask
    fileNameMask=[char(listOfCellMasks(i,2)) filesep char(listOfCellMasks(i,1))];
    if matformatcm==0 % tiff
        cm=imread(fileNameMask); %cm is cell mask
    else
        cm=load(fileNameMask);
        cm=cm.mask;
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
end

if nFms2avg==0 % average over all nFrames
    bgMean=nanmean(bgPixVal(:));
    bgStd=nanstd(bgPixVal(:));
else % do calculation for each frame
    for i=1:nFrames2Use
        sF=i-(nFms2avg-1)/2; %starting frame
        eF=i+(nFms2avg-1)/2; %ending frame

        if sF<1
            sF=1;
        end
        if eF>nFrames2Use
            eF=nFrames2Use;
        end

        usefulPart=bgPixVal(:,:,sF:eF);
        bgMean(i)=nanmean(usefulPart(:));
        bgStd(i)=nanstd(usefulPart(:));
    end
end
