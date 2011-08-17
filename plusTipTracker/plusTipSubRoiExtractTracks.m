function [projData,M]=plusTipSubRoiExtractTracks(subRoiDir,timeUnits,timeVal)
% this fn is called by plusTipSubRoiTool
%
% NOTE: This Function was MODIFIED from the original release of
% plusTipTracker.  Now partitions pause and shrinkage information from the
% given subRegions.  Also fixes small bug in original code where the 
% Fraction option was ALWAYS bi-passed, even if it was 
% selected in the GUI, due to a small capitalization error. 
% Maria Bagonis (MB) 03/27/11
%
% 
%

%% new input to be incorporated later
onlyTarget =  0; % If set to one will select only those tracks that are 
%targeted (ie the last point of their growth track is within the given sub
%region if set to zero will use the settings specified in the
%plusTipSeeTracks GUI
useCropped = 1;

%% Check Input
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

%% Load Necessary Files
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

%% Extract Data: Growth

% notes : tech could change this to just read in the projData converted
% which is already saved (projData.dataMatForStats should do the trick
% here) though not sure if the way i did the indexing would allow for this 
% might be safer to just use the original data- then crop if user asks 



% get all the original merged tracks and convert frames to seconds and
% pixels to microns % this merging step is redundant we could maybe save in
% original analysis and read in

[dataMatMergeAll] =plusTipMergeSubtracks(sourceProjData); % merged data is first output


allDataAll=dataMatMergeAll; % save a separate matrix for converting
%Convert Growth Lifetimes Frames to Seconds

allDataAll(:,6)=allDataAll(:,6).*sourceProjData.secPerFrame;
%Convert Growth Displacements From Pixels to Microns
allDataAll(:,7)=allDataAll(:,7).*(sourceProjData.pixSizeNm./1000);




% get growth track indices only and coordinates for those tracks only
idx=find(allDataAll(:,5)==1);
dataMatMerge=dataMatMergeAll(idx,:); % just growth subtracks
allData=allDataAll(idx,:);
[xMat,yMat]=plusTipGetSubtrackCoords(sourceProjData,idx,1);
%Note: xMat: x-coordinate of detected particle: row # = track ID, col.# = frame number

% get which tracks have their first point NOT in the exclude region
%firstPtFrIdx gives the start frame for each growth subtrack
firstPtFrIdx=arrayfun(@(i) find(~isnan(xMat(i,:)),1,'first'),1:length(idx))';
%Extract the Coordinates Corresponding to the first particle of 
% each subtrack
x=ceil(arrayfun(@(i) xMat(i,firstPtFrIdx(i)),1:length(idx))'-.5);
y=ceil(arrayfun(@(i) yMat(i,firstPtFrIdx(i)),1:length(idx))'-.5);

%Convert from xy coordinate to pixel index
pixIdx=sub2ind([imL,imW],y,x); % pixel index of all the starts
%Get List of those start sites in the larger region of interest
%(Makes sure if the user excluded any region that these tracks will 
% not be considered
inIncludeRegion=find(~excludeMask(pixIdx));

if onlyTarget == 1;
    
%Find tracks targeted to sub region of interest
lastPtFrIdx=arrayfun(@(i) find(~isnan(xMat(i,:)),1,'last'),1:length(idx))';
firstPtFrIdx=arrayfun(@(i) find(~isnan(xMat(i,:)),1,'first'),1:length(idx))';

%Extract the Coordinates Corresponding to the last particle 
% of each track
xLast=ceil(arrayfun(@(i) xMat(i,lastPtFrIdx(i)),1:length(idx))'-.5);
yLast=ceil(arrayfun(@(i) yMat(i,lastPtFrIdx(i)),1:length(idx))'-.5);

xFirst=ceil(arrayfun(@(i) xMat(i,firstPtFrIdx(i)),1:length(idx))'-.5);
yFirst=ceil(arrayfun(@(i) yMat(i,firstPtFrIdx(i)),1:length(idx))'-.5);



% Do same for last point
%Convert from xy coordinate to pixel index 
pixIdxLast=sub2ind([imL,imW],yLast,xLast);% pixel index of all last points 
% of track 

pixIdxFirst=sub2ind([imL,imW],yFirst,xFirst);

%find which tracks have their last point within the subregion of interest
inIncludeRegionLast=find(roiMask(pixIdxLast));
inIncludeRegionFirst=find(roiMask(pixIdxFirst));
else % No need to calculate these variables
end

%inIncludeRegionFirst=find(roiMask(pixIdx));


% find which of the TRACK indices (not just start sites) are in the roiMask
% assume that the first pixel in the image will be a zero
% but make sure this is so by making the first roiMask pixel = 0
%Convert from xy indices to pixel index
%pixIdx maintains structure of xMat/yMat ie row = subTrack ID 
% and column = frame number.  Instead of x/y coordinates now point 
% denoted by pixel number in a imL by imW image. 
pixIdx=sub2ind([imL,imW],ceil(yMat-.5),ceil(xMat-.5));
pixIdx(isnan(pixIdx))=1;
roiMask(1,1)=0;

% IN is a matrix of row = subtrack ID and col = frame number
% places a 1 in those frames where the point is inside the 
% roiMask and zero where it is not
IN=roiMask(pixIdx);

lifeSec=allData(:,6); % track lifetimes, seconds
%the lifetieme(ie number of frames) 
%inside the ROI is just the sum of each row (sum in the 2nd D)-1
insideSec=sourceProjData.secPerFrame.*(sum(IN,2)-1);

if onlyTarget == 1 
    
%find track coordinates cooresponding to those tracks that are targeted 
%to a given subregion

trckIdxInFirst = intersect(find(insideSec>0),inIncludeRegionFirst);
xMatFirst = xMat(trckIdxInFirst,:);
yMatFirst = yMat(trckIdxInFirst,:);

trckIdxIn= intersect(find(insideSec>0),inIncludeRegionLast);
xMatLast = xMat(trckIdxIn,:);
yMatLast = yMat(trckIdxIn,:);




    IN=IN(trckIdxIn,:); 
    % Just exchange zero values for NaN, again IN is a matrix of the same 
    %form as xMat and yMat only replaced by 1 where the tracked coordinate should be 
    %be maintained in the subROI calc
    IN=swapMaskValues(IN,0,NaN); % 1 for the sections of the tracks inside the sub-roi
    OUT=swapMaskValues(IN); % 1 for the sections of the track outside the sub-roi

% plot all member tracks in red on top of the mask
 figure; 
 imshow(roiMask); 
 hold on; 
 plot(xMatFirst',yMatFirst','g'); % tracks initiated in subregion, not included in analysis
 hold on;
 plot(xMatLast',yMatLast','r') % tracks that end in subregion, included in analysis
 
 saveas(gcf,[subRoiDir filesep 'tracksTargetedToSubRoiRed_plusLeavingTracksGreen' fileExt])
 close(gcf)
 
else % get those tracks located within the subregion based with a given
    % time criteria (but NO requirement for whether the track starts or 
    % ends in the ROI of interest)

    % Get Track IDs for those Growth Sub-Tracks in the Sub-ROI of Interest
    if ~isempty(strmatch(timeUnits,'Fraction'))
        % Track index of those tracks you are including
        trckIdxIn=intersect(find(insideSec./lifeSec>=timeVal),inIncludeRegion);
        % Track index of those tracks that start in the region of interest but 
        % are not analyzed because do not make lifetime criteria. 
        trckIdxOut=intersect(find(insideSec./lifeSec<timeVal & insideSec>0),inIncludeRegion);
    else
        trckIdxIn=intersect(find(insideSec>=timeVal),inIncludeRegion);
        trckIdxOut=intersect(find(insideSec<timeVal & insideSec>0),inIncludeRegion);
    end

%find all tracks that start within a given subregion
%trckIdxInStart=intersect(find(insideSec>0),inIncludeRegionFirst);
%xMatStart=xMat(trckIdxInStart,:);
%yMatStart=yMat(trckIdxInStart,:);

    %limit data to these tracks
    xMatIn = xMat(trckIdxIn,:); 
    yMatIn = yMat(trckIdxIn,:);

    xMatOut = xMat(trckIdxOut,:);
    yMatOut = yMat(trckIdxOut,:);


    IN=IN(trckIdxIn,:); 
    % Just exchange zero values for NaN, again IN is a matrix of the same 
    %form as xMat and yMat only replaced by 1 where the tracked coordinate should be 
    %be maintained in the subROI calc
    IN=swapMaskValues(IN,0,NaN); % 1 for the sections of the tracks inside the sub-roi
    OUT=swapMaskValues(IN); % 1 for the sections of the track outside the sub-roi


% plot all member tracks in red on top of the mask
 figure; 
 imshow(roiMask); 
 hold on; 
 plot(xMatIn',yMatIn','r')
 saveas(gcf,[subRoiDir filesep 'tracksInSubRoi' fileExt])
 close(gcf)
 
%plot all tracks excluded from analysis
  figure
  imshow(roiMask);
  hold on;
  plot(xMatOut',yMatOut','g');
  saveas(gcf,[subRoiDir filesep 'tracksExcludedFromSubRoi' fileExt]);
  close(gcf)
 
end
 
%% Extract Data Pause
%if onlyTarget == 1  % Don't calculate the pause data 
 %   else % continue with pause calculations
%Get Pause Data 
idxPause=find(allDataAll(:,5)==2);
dataMatMergePause=dataMatMergeAll(idxPause,:);
allDataPause=allDataAll(idxPause,:);
[xMatPause,yMatPause]=plusTipGetSubtrackCoords(sourceProjData,idxPause,1);

% get which tracks have their first point NOT in the exclude region
%firstPtFrIdx gives the start frame for each growth subtrack
firstPtFrIdxPause=arrayfun(@(i) find(~isnan(xMatPause(i,:)),1,'first'),1:length(idxPause))';
%Extract the Coordinates Corresponding to the first particle of 
% each subtrack
xPause=ceil(arrayfun(@(i) xMatPause(i,firstPtFrIdxPause(i)),1:length(idxPause))'-.5);
yPause=ceil(arrayfun(@(i) yMatPause(i,firstPtFrIdxPause(i)),1:length(idxPause))'-.5);

%Convert from xy coordinate to pixel index
pixIdxPause=sub2ind([imL,imW],yPause,xPause); % pixel index of all the starts
%Get List of those start sites in the larger region of interest
%(Makes sure if the user excluded any region that these tracks will 
% not be considered
inIncludeRegionPause=find(~excludeMask(pixIdxPause));

pixIdxPause=sub2ind([imL,imW],ceil(yMatPause-.5),ceil(xMatPause-.5));
pixIdxPause(isnan(pixIdxPause))=1;
roiMask(1,1)=0;
% IN is a matrix of row = subtrack ID and col = frame number
% places a 1 in those frames where the point is inside the 
% roiMask and zero where it is not
INPause=roiMask(pixIdxPause);

lifeSecPause=allDataPause(:,6); % track lifetimes, seconds
%the lifetieme(ie number of frames) 
%inside the ROI is just the sum of each row (sum in the 2nd D)-1
insideSecPause =sourceProjData.secPerFrame.*(sum(INPause,2)-1);

% Get Track IDs for those Growth Sub-Tracks in the Sub-ROI of Interest
if ~isempty(strmatch(timeUnits,'Fraction'))
    trckIdxInPause=intersect(find(insideSecPause./lifeSecPause>=timeVal),inIncludeRegionPause);
    trckIdxOutPause=intersect(find(insideSecPause./lifeSecPause<timeVal & insideSecPause > 0),inIncludeRegionPause);
else
    trckIdxInPause=intersect(find(insideSec>=timeVal),inIncludeRegion);
    trckIdxOutPause=intersect(find(insideSec<timeVal & insideSecPause > 0),inIncludeRegion);
end

% limit data to these tracks
xMatPause=xMatPause(trckIdxInPause,:);
yMatPause=yMatPause(trckIdxInPause,:);
INPause=INPause(trckIdxInPause,:);
% Just exchange zero values for NaN, again IN is a matrix of the same 
%form as xMat and yMat only replaced by 1 where the tracked coordinate should be 
%be maintained in the subROI calc
INPause=swapMaskValues(INPause,0,NaN); % 1 for the sections of the tracks inside the sub-roi
OUTPause=swapMaskValues(INPause); % 1 for the sections of the track outside the sub-roi



% plot all member tracks in red on top of the mask
 figure; 
 imshow(roiMask); 
 hold on; 
 plot(xMatPause',yMatPause','b')
 saveas(gcf,[subRoiDir filesep 'PausesInSubRoi' fileExt])
 close(gcf)

%end
%% Extract Data Shrinkage

%if onlyTarget == 1 % Don't calculate the shrinkage data
 %   else % Continue with the shrinkage calculations

%Get Shrinkage Data

idxShrink=find(allDataAll(:,5)==3);
dataMatMergeShrink=dataMatMergeAll(idxShrink,:);
allDataShrink=allDataAll(idxShrink,:);
[xMatShrink,yMatShrink]=plusTipGetSubtrackCoords(sourceProjData,idxShrink,1);

% get which tracks have their first point NOT in the exclude region
%firstPtFrIdx gives the start frame for each growth subtrack
firstPtFrIdxShrink=arrayfun(@(i) find(~isnan(xMatShrink(i,:)),1,'first'),1:length(idxShrink))';
%Extract the Coordinates Corresponding to the first particle of 
% each subtrack
xShrink=ceil(arrayfun(@(i) xMatShrink(i,firstPtFrIdxShrink(i)),1:length(idxShrink))'-.5);
yShrink=ceil(arrayfun(@(i) yMatShrink(i,firstPtFrIdxShrink(i)),1:length(idxShrink))'-.5);

%Convert from xy coordinate to pixel index
pixIdxShrink=sub2ind([imL,imW],yShrink,xShrink); % pixel index of all the starts
%Get List of those start sites in the larger region of interest
%(Makes sure if the user excluded any reg1.0ion that these tracks will 
% not be considered
inIncludeRegionShrink=find(~excludeMask(pixIdxShrink));

pixIdxShrink=sub2ind([imL,imW],ceil(yMatShrink-.5),ceil(xMatShrink-.5));
pixIdxShrink(isnan(pixIdxShrink))=1;
roiMask(1,1)=0;
% IN is a matrix of row = subtrack ID and col = frame number
% places a 1 in those frames where the point is inside the 
% roiMask and zero where it is not
INShrink=roiMask(pixIdxShrink);

lifeSecShrink=allDataShrink(:,6); % track lifetimes, seconds
%the lifetieme(ie number of frames) 
%inside the ROI is just the sum of each row (sum in the 2nd D)-1
insideSecShrink =sourceProjData.secPerFrame.*(sum(INShrink,2)-1);

% Get Track IDs for those Growth Sub-Tracks in the Sub-ROI of Interest
if ~isempty(strmatch(timeUnits,'Fraction'))
    trckIdxInShrink=intersect(find(insideSecShrink./lifeSecShrink>=timeVal),inIncludeRegionShrink);
    trckIdxOutShrink=intersect(find(insideSecShrink./lifeSecShrink<timeVal & insideSecShrink > 0),inIncludeRegionShrink);
else
    trckIdxInShrink=intersect(find(insideSec>=timeVal),inIncludeRegion); 
    trckIdxOutShrink = intersect(find(insideSecShrink<timeVal & insideSecShrink > 0),inIncludeRegion);
end

% limit data to these tracks
xMatShrink=xMatShrink(trckIdxInShrink,:);
yMatShrink=yMatShrink(trckIdxInShrink,:);
INShrink=INShrink(trckIdxInShrink,:);
% Just exchange zero values for NaN, again IN is a matrix of the same 
%form as xMat and yMat only replaced by 1 where the tracked coordinate should be 
%be maintained in the subROI calc
INShrink=swapMaskValues(INShrink,0,NaN); % 1 for the sections of the tracks inside the sub-roi
OUTShrink=swapMaskValues(INShrink); % 1 for the sections of the track outside the sub-roi

% plot all member tracks in red on top of the mask
 figure; 
 imshow(roiMask); 
 hold on; 
 plot(xMatShrink',yMatShrink','bl')
 saveas(gcf,[subRoiDir filesep 'ShrinkagesInSubRoi' fileExt])
 close(gcf)
 
%end
%% Write Stats to Output

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


% there are no track numbers that contain an fgap or bgap- NOW THERE IS!!
% FIX
projData.percentFgapsReclass=NaN;
projData.percentBgapsReclass=NaN;
projData.tracksWithFgap = NaN;
projData.tracksWithBgap = NaN;

% calculate stats using the matrix where beginning/end data have NOT
% been removed. M records speeds (microns/min), lifetimes (sec), and
% displacements (microns) for growths, fgaps,and bgaps (of which the
% latter two do not exist here)

%if onlyTarget == 1
%Collect Information Regarding all SubTracks Within 
% ROI
%dataTotROI = dataMatMerge(trckIdxIn,:); % data total 

%[projData.stats,M]=plusTipDynamParam(dataTotROI);

%else

dataGrowthROI = dataMatMerge(trckIdxIn,:);
dataPauseROI = dataMatMergePause(trckIdxInPause,:);
dataShrinkROI = dataMatMergeShrink(trckIdxInShrink,:);

dataTotROI = [dataGrowthROI ; dataPauseROI ; dataShrinkROI];
dataTotROI_Frames_Pix = sortrows(dataTotROI);

%Convert (allDataAll is likewise converted from frame to sec and pix to um)
% however lets just do it again directly for now so I make sure indexing is
% correct (CHECK LATER)
dataTotROI_Sec_Mic = dataTotROI_Frames_Pix;
dataTotROI_Sec_Mic(:,6) = dataTotROI_Frames_Pix(:,6).*sourceProjData.secPerFrame;
dataTotROI_Sec_Mic(:,7) = dataTotROI_Sec_Mic(:,7).*sourceProjData.pixSizeNm./1000;

%Those tracks with start sites in the region of interest 
%but were excluded based on user specified time requirement 
%dataGrowthExclude = dataMatMerge(trckIdxOut,:);
%dataPauseExclude = dataMatMerge(trckIdxOutPause,:);
%dataShrinkExclude = dataMatMergeShrink(trckIdxOutShrink,:);

%dataTotExclude = [dataGrowthExclude ; dataPauseExclude];
%dataTotExclude = [dataTotExclude ; dataShrinkExclude];


% index of one at the end tells the program it is calling this function
% from subRoiAnalysis and to not calculate certain params and these might
% be buggy

if useCropped == 1
[dataTotROI_Sec_Mic] = plusTipRemBegEnd(dataTotROI_Sec_Mic, projData);
 projData.dataMatCropped_MicPerMin_Sec_Mic = dataTotROI_Sec_Mic; % save a copy
 projData.removeBegEnd = 'yes'; 
else 
     projData =  rmfield(projData,'dataMatCropped_MicPerMin_Sec_Mic'); % remove the field
     projData.removeBegEnd = 'no'; 
end 
    [projData,M]=plusTipDynamParam(dataTotROI_Sec_Mic,projData,0,1);
   
    
%projData.tracksExcludedFromAnalysis = dataTotExclude;

%end


    

% put here the merged data in frames and pixels
projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix=dataTotROI_Frames_Pix;

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

    %  
    if onlyTarget ~= 1
   
    speedIn=nanmean(sqrt(diff(xMatIn.*IN,[],2).^2+diff(yMatIn.*IN,[],2).^2),2);
    speedOut=nanmean(sqrt(diff(xMatIn.*OUT,[],2).^2+diff(yMatIn.*OUT,[],2).^2),2);

    
    projData.speedInMicPerMin=pixPerFrame2umPerMin(speedIn,projData.secPerFrame,projData.pixSizeNm);
    projData.speedOutMicPerMin=pixPerFrame2umPerMin(speedOut,projData.secPerFrame,projData.pixSizeNm);

    projData.startOrEnd=~isnan(xMat(:,1)) | ~isnan(xMat(:,end));
    projData.percentAtStartOrEnd=sum(projData.startOrEnd)./projData.nTracks;
    
    else % only calculate at this point if we aren't doing targeting I should 
        % go back and look at this as this value might be helpful 
    end 
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