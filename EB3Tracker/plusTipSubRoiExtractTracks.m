function [projData,M]=plusTipSubRoiExtractTracks(subRoiDir,excludeMask)
% this fn is called by plusTipSubdivideRoi or as a standalone fn after
% re-running tracking if you don't want to re-draw ROIs

% INPUT: subRoiDir: path of sub_x folder of interest
% OUTPUT: projData: project data similar to post-processing output, only
% for growth tracks (no fgaps or bgaps) that exist for 3 or more frames
% inside the region.  There are some new fields with info about the track
% mean speeds in/out of the ROI at the end.

homeDir=pwd;

subRoiDir=formatPath(subRoiDir);
outputDir=[subRoiDir filesep 'meta'];
if ~isdir(outputDir)
    mkdir(outputDir)
end


if nargin<2 || isempty(excludeMask)
    excludeMask=[];
end

cd(subRoiDir)
cd ..
cd ..
cd 'meta'
load 'projData'
sourceProjData=projData;
clear('projData')
mainProjDir=pwd;



% all the original tracks
aT=sourceProjData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix;

% get growth track indices only and coordinates
idx=find(aT(:,5)==1);
[xMat,yMat]=plusTipGetSubtrackCoords(sourceProjData,idx);

% get which tracks have their first point NOT in the exclude region
if ~isempty(excludeMask)
    [imL,imW]=size(excludeMask);
    x=zeros(length(idx),1);
    y=zeros(length(idx),1);
    for i=1:length(idx)
        x(i)=ceil(xMat(i,find(~isnan(xMat(i,:)),1,'first'))-.5);
        y(i)=ceil(yMat(i,find(~isnan(yMat(i,:)),1,'first'))-.5);
    end
    pixIdx=sub2ind([imL,imW],y,x); % pixel index of all the starts
    inIncludeRegion=find(~excludeMask(pixIdx));
else
    inIncludeRegion=1:length(idx);
end

% figure;
% imshow(~excludeMask)
% hold on
% scatter(x(inIncludeRegion),y(inIncludeRegion))

% start/end frames based on range stored in projData
sF=max([sourceProjData.detectionFrameRange(1); sourceProjData.trackingFrameRange(1); sourceProjData.postTrackFrameRange(1)]);
eF=min([sourceProjData.detectionFrameRange(2); sourceProjData.trackingFrameRange(2); sourceProjData.postTrackFrameRange(2)]);

% limit tracks to only these frames
xMat=xMat(:,sF:eF);
yMat=yMat(:,sF:eF);

% load sub-roi outline pixels
%roiYX=load([subRoiDir filesep 'roiYX']);
%yv=roiYX.roiYX(:,1);
%xv=roiYX.roiYX(:,2);

roiMask=imread([subRoiDir filesep 'roiMask.tif']);
pixIdx=sub2ind([imL,imW],ceil(yMat-.5),ceil(xMat-.5));
% assume that the first pixel in the image will be a zero
% but make sure this is so by making the first roiMask pixel = 0
pixIdx(isnan(pixIdx))=1; 
roiMask(1,1)=0;

% find which of the track indices are in the roiMask
IN=roiMask(pixIdx);


% calculate which growth track features are in the sub-ROI
%IN1=inpolygon(xMat,yMat,xv,yv);

% track indices where tracks exist for at least three time points
trckIdxIn=intersect(find(sum(IN,2)>=3),inIncludeRegion);

% limit data to these tracks
xMat=xMat(trckIdxIn,:);
yMat=yMat(trckIdxIn,:);
IN=IN(trckIdxIn,:);
IN=swapMaskValues(IN,0,NaN); % 1 for the sections of the tracks inside the sub-roi
OUT=swapMaskValues(IN); % 1 for the sections of the track outside the sub-roi


projData=sourceProjData;
projData.anDir=subRoiDir; % path to sub-roi
projData.imDir=formatPath(projData.imDir);
projData.numTracks=length(trckIdxIn); % number of growth tracks

% keep only the coordinates, speeds, etc. corresponding to tracks remaining
projData.xCoord=nan(size(sourceProjData.xCoord));
projData.yCoord=nan(size(sourceProjData.yCoord));
projData.featArea=nan(size(sourceProjData.featArea));
projData.featInt=nan(size(sourceProjData.featInt));
projData.frame2frameVel_micPerMin=nan(size(sourceProjData.frame2frameVel_micPerMin));
projData.segGapAvgVel_micPerMin=nan(size(sourceProjData.segGapAvgVel_micPerMin));
for iSub=1:length(trckIdxIn)
    projData.xCoord(aT(idx(trckIdxIn(iSub)),1),aT(idx(trckIdxIn(iSub)),2):aT(idx(trckIdxIn(iSub)),3))=sourceProjData.xCoord(aT(idx(trckIdxIn(iSub)),1),aT(idx(trckIdxIn(iSub)),2):aT(idx(trckIdxIn(iSub)),3));
    projData.yCoord(aT(idx(trckIdxIn(iSub)),1),aT(idx(trckIdxIn(iSub)),2):aT(idx(trckIdxIn(iSub)),3))=sourceProjData.yCoord(aT(idx(trckIdxIn(iSub)),1),aT(idx(trckIdxIn(iSub)),2):aT(idx(trckIdxIn(iSub)),3));

    projData.featArea(aT(idx(trckIdxIn(iSub)),1),aT(idx(trckIdxIn(iSub)),2):aT(idx(trckIdxIn(iSub)),3))=sourceProjData.featArea(aT(idx(trckIdxIn(iSub)),1),aT(idx(trckIdxIn(iSub)),2):aT(idx(trckIdxIn(iSub)),3));
    projData.featInt(aT(idx(trckIdxIn(iSub)),1),aT(idx(trckIdxIn(iSub)),2):aT(idx(trckIdxIn(iSub)),3))=sourceProjData.featInt(aT(idx(trckIdxIn(iSub)),1),aT(idx(trckIdxIn(iSub)),2):aT(idx(trckIdxIn(iSub)),3));

    projData.frame2frameVel_micPerMin(aT(idx(trckIdxIn(iSub)),1),aT(idx(trckIdxIn(iSub)),2):aT(idx(trckIdxIn(iSub)),3)-1)=sourceProjData.frame2frameVel_micPerMin(aT(idx(trckIdxIn(iSub)),1),aT(idx(trckIdxIn(iSub)),2):aT(idx(trckIdxIn(iSub)),3)-1);
    projData.segGapAvgVel_micPerMin(aT(idx(trckIdxIn(iSub)),1),aT(idx(trckIdxIn(iSub)),2):aT(idx(trckIdxIn(iSub)),3)-1)=sourceProjData.segGapAvgVel_micPerMin(aT(idx(trckIdxIn(iSub)),1),aT(idx(trckIdxIn(iSub)),2):aT(idx(trckIdxIn(iSub)),3)-1);
end

if projData.numTracks~=0
    % get frame-to-frame displacement for growth only (not forward/backward gaps)
    frame2frameDispPix=sqrt(diff(projData.xCoord,1,2).^2+diff(projData.yCoord,1,2).^2);
    % get rid of NaNs and linearize the vector
    projData.frame2frameDispPix=frame2frameDispPix(~isnan(frame2frameDispPix(:)));
else
    projData.frame2frameDispPix=NaN;
end

if projData.numTracks~=0
    % get change in velocity between frame *pairs* for segments only
    pair2pairDiffPix=diff(frame2frameDispPix,1,2);
    % get rid of NaNs and linearize the vector
    projData.pair2pairDiffPix=pair2pairDiffPix(~isnan(pair2pairDiffPix(:)));
else
    projData.pair2pairDiffPix=NaN;
end
% std (microns/min) of delta growthSpeed btw frames
projData.pair2pairDiffMicPerMinStd=std(pixPerFrame2umPerMin(projData.pair2pairDiffPix,projData.secPerFrame,projData.pixSizeNm));

projData.medNNdistWithinFramePix=NaN;
projData.meanDisp2medianNNDistRatio=NaN;

% there are no track numbers that contain an fgap or bgap
projData.percentFgapsReclass=NaN;
projData.percentBgapsReclass=NaN;
projData.tracksWithFgap = NaN;
projData.tracksWithBgap = NaN;

% calculate stats using the matrix where beginning/end data have NOT
% been removed. M records speeds (microns/min), lifetimes (sec), and
% displacements (microns) for growths, fgaps,and bgaps (of which the
% latter two do not exist here)
aT=aT(idx(trckIdxIn),:);
aTSecMic=aT;
aTSecMic(:,6)=aT(:,6).*projData.secPerFrame;
aTSecMic(:,7)=aT(:,7).*(projData.pixSizeNm/1000);

[projData.stats,M]=plusTipDynamParam(aTSecMic);

projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix=aT;



% NEW STUFF NOT IN SOURCE PROJDATA.MAT
if projData.numTracks~=0
    projData.sourceSubtrackIdx=idx(trckIdxIn); % original subtrack index from sourceProjData.nTrack...
    projData.trackLifeFrames=sum(~isnan(xMat),2); % total lifetime (frames)
    projData.framesInSubRoi=nansum(IN,2); % lifetime within sub-roi (frames)
    projData.percentLifeInside=100*(projData.framesInSubRoi./projData.trackLifeFrames); % percent time within sub-roi

    speedIn=nanmean(sqrt(diff(xMat.*IN,[],2).^2+diff(yMat.*IN,[],2).^2),2);
    speedOut=nanmean(sqrt(diff(xMat.*OUT,[],2).^2+diff(yMat.*OUT,[],2).^2),2);

    projData.speedInMicPerMin=pixPerFrame2umPerMin(speedIn,projData.secPerFrame,projData.pixSizeNm);
    projData.speedOutMicPerMin=pixPerFrame2umPerMin(speedOut,projData.secPerFrame,projData.pixSizeNm);
    projData.trackLifeSec=projData.trackLifeFrames.*projData.secPerFrame;

    projData.startOrEnd=~isnan(xMat(:,1)) | ~isnan(xMat(:,end));
    projData.percentAtStartOrEnd=sum(projData.startOrEnd)./projData.numTracks;
else
    projData.sourceSubtrackIdx=NaN;
    projData.trackLifeFrames=NaN;
    projData.framesInSubRoi=NaN;
    projData.percentLifeInside=NaN;

    projData.speedInMicPerMin=NaN;
    projData.speedOutMicPerMin=NaN;
    projData.trackLifeSec=NaN;

    projData.startOrEnd=NaN;
    projData.percentAtStartOrEnd=NaN;
end


% save projData in meta folder
save([outputDir filesep 'projData'],'projData')

cd(homeDir)