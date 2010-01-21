function [projData,M]=plusTipSubRoiExtractTracks(subRoiDir,timeUnits,timeVal)
% this fn is called by plusTipSubRoiTool
%
% NOTE: the extracted tracks are the merged growth tracks only - no fgaps or
% bgaps.

homeDir=pwd;

if isempty(strfind(subRoiDir,'sub'))
    return
end

if ispc
    fileExt='.emf';
else
    fileExt='.tif';
end

% check input
if ~ismember(lower(timeUnits),{'fraction','seconds'})
    error('plusTipSubRoiTool: timeUnits must be fraction or seconds')
end
if isempty(timeVal)
    error('plusTipSubRoiTool: time input missing')
end
if ~isempty(strmatch(lower(timeUnits),'fraction')) && ~(timeVal>0 && timeVal<=1)
    error('plusTipSubRoiTool: timeUnits is fraction, timeVal must be in 0-1')
end

subRoiDir=subRoiDir;

cd(subRoiDir) % sub_x directory
cd ..
subanDir=pwd; % subROIs directory
cd ..
cd 'meta' % roi_x/meta directory
load 'projData'
sourceProjData=projData; % projData from roi_x
clear('projData')

% wholeCellRoiMask is the roi_x mask
wholeCellRoiMask=imread([subanDir filesep 'roiMask.tif']);
% roiMask is the sub_x mask
roiMask=imread([subRoiDir filesep 'roiMask.tif']);
[imL,imW]=size(roiMask);

% if there is an exclude mask, use it; otherwise, use the inverse of the
% whole cell mask and write that into the sub_x folder
if exist([subRoiDir filesep 'excludeMask.tif'],'file')~=0
    excludeMask=imread([subRoiDir filesep 'excludeMask.tif']);
else
    excludeMask=~wholeCellRoiMask;
    imwrite(excludeMask,[subRoiDir filesep 'excludeMask.tif']);
end


% get all the original merged tracks and convert frames to seconds and
% pixels to microns
dataMatMerge=plusTipMergeSubtracks(sourceProjData); % merged data is first output
allData=abs(dataMatMerge);
allData(:,6)=allData(:,6).*sourceProjData.secPerFrame;
allData(:,7)=allData(:,7).*(sourceProjData.pixSizeNm./1000);

% get growth track indices only and coordinates for those tracks only
idx=find(allData(:,5)==1);
dataMatMerge=dataMatMerge(idx,:);
allData=allData(idx,:);
[xMat,yMat]=plusTipGetSubtrackCoords(sourceProjData,idx,1);

% get which tracks have their first point NOT in the exclude region
firstPtFrIdx=arrayfun(@(i) find(~isnan(xMat(i,:)),1,'first'),1:length(idx))';
x=ceil(arrayfun(@(i) xMat(i,firstPtFrIdx(i)),1:length(idx))'-.5);
y=ceil(arrayfun(@(i) yMat(i,firstPtFrIdx(i)),1:length(idx))'-.5);

pixIdx=sub2ind([imL,imW],y,x); % pixel index of all the starts
inIncludeRegion=find(~excludeMask(pixIdx));

% find which of the track indices are in the roiMask
% assume that the first pixel in the image will be a zero
% but make sure this is so by making the first roiMask pixel = 0
pixIdx=sub2ind([imL,imW],ceil(yMat-.5),ceil(xMat-.5));
pixIdx(isnan(pixIdx))=1;
roiMask(1,1)=0;
IN=roiMask(pixIdx);

lifeSec=allData(:,6); % track lifetimes, seconds
insideSec=sourceProjData.secPerFrame.*(sum(IN,2)-1);

% calculate which growth track features are in the sub-ROI
if ~isempty(strmatch(timeUnits,'fraction'))
    trckIdxIn=intersect(find(insideSec./lifeSec>=timeVal),inIncludeRegion);
else
    trckIdxIn=intersect(find(insideSec>=timeVal),inIncludeRegion);
end

% limit data to these tracks
xMat=xMat(trckIdxIn,:);
yMat=yMat(trckIdxIn,:);
IN=IN(trckIdxIn,:);
IN=swapMaskValues(IN,0,NaN); % 1 for the sections of the tracks inside the sub-roi
OUT=swapMaskValues(IN); % 1 for the sections of the track outside the sub-roi

% plot all member tracks in red on top of the mask
% figure; 
% imshow(roiMask); 
% hold on; 
% plot(xMat',yMat','r')
% saveas(gcf,[subRoiDir filesep 'tracksInSubRoi' fileExt])
% close(gcf)

projData=sourceProjData;
projData.anDir=subRoiDir; % path to sub-roi
projData.imDir=projData.imDir;
projData.nTracks=length(trckIdxIn); % number of growth tracks

% keep only the coordinates, speeds, etc. corresponding to tracks remaining
projData.xCoord=nan(size(sourceProjData.xCoord));
projData.yCoord=nan(size(sourceProjData.yCoord));
projData.featArea=nan(size(sourceProjData.featArea));
projData.featInt=nan(size(sourceProjData.featInt));
projData.frame2frameVel_micPerMin=nan(size(sourceProjData.frame2frameVel_micPerMin));
projData.segGapAvgVel_micPerMin=nan(size(sourceProjData.segGapAvgVel_micPerMin));

for iSub=1:length(trckIdxIn)
    k=trckIdxIn(iSub);

    projData.xCoord(allData(k,1),allData(k,2):allData(k,3))=sourceProjData.xCoord(allData(k,1),allData(k,2):allData(k,3));
    projData.yCoord(allData(k,1),allData(k,2):allData(k,3))=sourceProjData.yCoord(allData(k,1),allData(k,2):allData(k,3));

    projData.featArea(allData(k,1),allData(k,2):allData(k,3))=sourceProjData.featArea(allData(k,1),allData(k,2):allData(k,3));
    projData.featInt(allData(k,1),allData(k,2):allData(k,3))=sourceProjData.featInt(allData(k,1),allData(k,2):allData(k,3));

    projData.frame2frameVel_micPerMin(allData(k,1),allData(k,2):allData(k,3)-1)=sourceProjData.frame2frameVel_micPerMin(allData(k,1),allData(k,2):allData(k,3)-1);
    projData.segGapAvgVel_micPerMin(allData(k,1),allData(k,2):allData(k,3)-1)=sourceProjData.segGapAvgVel_micPerMin(allData(k,1),allData(k,2):allData(k,3)-1);
end

if projData.nTracks~=0
    % get frame-to-frame displacement for growth only (not forward/backward gaps)
    frame2frameDispPix=sqrt(diff(projData.xCoord,1,2).^2+diff(projData.yCoord,1,2).^2);
    % get rid of NaNs and linearize the vector
    projData.frame2frameDispPix=frame2frameDispPix(~isnan(frame2frameDispPix(:)));
else
    projData.frame2frameDispPix=NaN;
end

if projData.nTracks~=0
    % get change in velocity between frame *pairs* for segments only
    pair2pairDiffPix=diff(frame2frameDispPix,1,2);
    % get rid of NaNs and linearize the vector
    projData.pair2pairDiffPix=pair2pairDiffPix(~isnan(pair2pairDiffPix(:)));
else
    projData.pair2pairDiffPix=NaN;
end
% std (microns/min) of delta growthSpeed btw frames
projData.pair2pairDiffMicPerMinStd=std(pixPerFrame2umPerMin(projData.pair2pairDiffPix,...
    projData.secPerFrame,projData.pixSizeNm));

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
[projData.stats,M]=plusTipDynamParam(allData);

% put here the merged data in frames and pixels
projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix=dataMatMerge(trckIdxIn,:);

% NEW STUFF NOT IN SOURCE PROJDATA.MAT
if projData.nTracks~=0
    pixSizMic=projData.pixSizeNm/1000; % side of a pixel in microns
    pixAreaSqMic=pixSizMic^2; % area of a pixel in square microns
    cellAreaSqMic=sum(wholeCellRoiMask(:)).*pixAreaSqMic;
    % calculate area in square microns of this roi and its
    % percentage of the full cell's area
    subRoiAreaSqMic=sum(roiMask(:)).*pixAreaSqMic;
    percentRoiArea=100*(subRoiAreaSqMic/cellAreaSqMic);
    projData.subRoiAreaSqMic=subRoiAreaSqMic;
    projData.percentRoiArea=percentRoiArea;

    projData.lifeSec=lifeSec(trckIdxIn); % total lifetime (seconds)
    projData.insideSec=insideSec(trckIdxIn); % lifetime within sub-roi (seconds)
    projData.percentLifeInside=100*(projData.insideSec./projData.lifeSec); % percent time within sub-roi

    speedIn=nanmean(sqrt(diff(xMat.*IN,[],2).^2+diff(yMat.*IN,[],2).^2),2);
    speedOut=nanmean(sqrt(diff(xMat.*OUT,[],2).^2+diff(yMat.*OUT,[],2).^2),2);

    projData.speedInMicPerMin=pixPerFrame2umPerMin(speedIn,projData.secPerFrame,projData.pixSizeNm);
    projData.speedOutMicPerMin=pixPerFrame2umPerMin(speedOut,projData.secPerFrame,projData.pixSizeNm);

    projData.startOrEnd=~isnan(xMat(:,1)) | ~isnan(xMat(:,end));
    projData.percentAtStartOrEnd=sum(projData.startOrEnd)./projData.nTracks;
else
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
save([subRoiDir filesep 'meta' filesep 'projData'],'projData')
% write out speed/lifetime/displacement distributions into a text file
dlmwrite([subRoiDir filesep 'meta' filesep 'gs_fs_bs_gl_fl_bl_gd_fd_bd.txt'], M, 'precision', 3,'delimiter', '\t','newline', 'pc');

cd(homeDir)