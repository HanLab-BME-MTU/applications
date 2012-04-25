function [projData,M]=plusTipSubRoiExtractTracks(subRoiDir,timeUnits,timeVal)
% this fn is called by plusTipSubRoiTool
%
% NOTE: This Function was MODIFIED from the original release of
% plusTipTracker.  Now partitions pause and shrinkage information from the
% given subRegions.  Also fixes small bug in original code where the 
% Fraction option was ALWAYS bi-passed, even if it was 
% selected in the GUI, due to a small capitalization error. 
% Furthermore, it has internal option to partition subtracks based on if 
% their start site (onlyInitiate option) or end site (onlyTarget)
% is located within the region of interest
% Maria Bagonis (MB) 03/27/11 
% 
%
% 
%
%% new input to be incorporated later
onlyTarget =  1; % 

% onlyTarget: Value of 1: flag to bypass the GUI options and 
% select only those growth SUBTRACKS that are 
% targeted to a subROI
% (ie the last point of their growth track is within the given sub
% region). Subdivision of bgap and
% fgaps among subRois will still be 
% based on the GUI specifications. 
% Value of 0: Use the settings specified in the
% plusTipGroupAnalysis GUI for subdividing all subtrack types

onlyInitiate = 0; % Value of 1: Flag to bypass the GUI options and select
% only those SUBTRACKS that are initiated (ie the first point of their 
% growth track is within the given subregion). Subdivision of bgap and
% fgaps among subRois will still be based on the GUI input specifications. 
% Value of 0 will use the settings specified in the
% plusTipGroupAnalysis GUI for subdividing all subtracks types

onlyNuc = 0;% value of 1: Flag to bypass the GUI options and select
% only those SUBTRACKS that are NUCLEATED (ie the first point of their 
% growth track is within the given subregion and NOT after a bgap or fgap). 
% Subdivision of bgap and
% fgaps among subRois will still be based on the GUI input specifications. 
% Value of 0 will use the settings specified in the
% plusTipGroupAnalysis GUI for subdividing all subtracks types


remBegEnd = 1; % Value of 1 will remove tracks from the beginning and end

collectPlots = 1; % collect all subRoi track plots in one folder (series of 
%tifs that can be read into imagej for easy viewing of each cell). 

showResultsProc = 1 ; % show results during processing...will still collect these plots just will not pop-up 
                  % for troubleshooting during processing. 

%% Check Input
homeDir=pwd;

if isempty(strfind(subRoiDir,'sub'))
    return
end

fileExt='.tif';


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
if onlyNuc ==1 
    onlyInitiate = 1; 
end 

if showResultsProc == 1 
   visible = 'on'; 
else 
   visible = 'off'; 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

up1 = getFilenameBody(subRoiDir);
up2 = getFilenameBody(up1); 
[up3 roi numSub] = getFilenameBody(up2); 
[collectedDataPath folder numProj] = getFilenameBody(up3); 
collectedDataPath = [collectedDataPath filesep 'collectedSubRoiPlots']; 
roi = [roi numSub]; 
name = [folder numProj]; 
projName = [name roi];
projNameTitle = regexprep(projName,'_',' '); 


%% Collect Data to Be Subdivided

% get all the original merged tracks and convert frames to seconds and
% pixels to microns

% Reads in merged output
dataMatMergeAll = sourceProjData.mergedDataMatAllSubTracksConverted; 


%% Extract Growth SubTracks

% get growth track indices only and coordinates for those tracks only
idx=find(dataMatMergeAll(:,5)==1);
dataMatMerge=dataMatMergeAll(idx,:); % just growth subtracks

%Note: xMat: x-coordinate of detected particle: row # = subtrack ID, col.# = frame number
% 1 indicates that it will use merged data
[xMat,yMat]=plusTipGetSubtrackCoords(sourceProjData,idx,1);

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

% Use this index to make sure if the user excluded any region 
% that these tracks will not be considered (in the case of no selection 
% the exclude mask will be outside of the ROI) 
inIncludeRegion=find(~excludeMask(pixIdx));


if (onlyTarget == 1 || onlyInitiate == 1);
    
%Find tracks targeted to sub region of interest
lastPtFrIdx=arrayfun(@(i) find(~isnan(xMat(i,:)),1,'last'),1:length(idx))';
firstPtFrIdx=arrayfun(@(i) find(~isnan(xMat(i,:)),1,'first'),1:length(idx))';

%Extract the Coordinates Corresponding to the Last/First particle 
% of each track
xLast=ceil(arrayfun(@(i) xMat(i,lastPtFrIdx(i)),1:length(idx))'-.5);
yLast=ceil(arrayfun(@(i) yMat(i,lastPtFrIdx(i)),1:length(idx))'-.5);

xFirst=ceil(arrayfun(@(i) xMat(i,firstPtFrIdx(i)),1:length(idx))'-.5);
yFirst=ceil(arrayfun(@(i) yMat(i,firstPtFrIdx(i)),1:length(idx))'-.5);

%Convert from xy coordinate to pixel index 
pixIdxLast=sub2ind([imL,imW],yLast,xLast);% pixel index of all last points 
% of track 

pixIdxFirst=sub2ind([imL,imW],yFirst,xFirst);

%find which tracks have their last point within the subregion of interest
% idx based on dataMatMerge
inIncludeRegionLast=find(roiMask(pixIdxLast));
inIncludeRegionFirst=find(roiMask(pixIdxFirst));
if onlyNuc == 1
    nucIdx = find(dataMatMerge(:,8) == 1);
    inIncludeRegionFirst = intersect(inIncludeRegionFirst,nucIdx); 
end 
else % No need to calculate these variables
end

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

lifeSec=dataMatMerge(:,6); % track lifetimes, seconds
%the lifetime(ie number of frames) 
%inside the ROI is just the sum of each row (sum in the 2nd D)-1
insideSec=sourceProjData.secPerFrame.*(sum(IN,2)-1);



if (onlyTarget == 1 || onlyInitiate ==1) 
    
%find track coordinates cooresponding to those tracks that are targeted 
%to a given subregion 

% Note inside sec have to be greater than 0 because can potentially 
% have only 1-coordinate in region: we would like to exclude these as they 
% are not useful for stats. 



trckIdxInFirst = intersect(find(insideSec>0),inIncludeRegionFirst);

xMatFirst = xMat(trckIdxInFirst,:);
yMatFirst = yMat(trckIdxInFirst,:);

trckIdxInLast= intersect(find(insideSec>0),inIncludeRegionLast);
xMatLast = xMat(trckIdxInLast,:);
yMatLast = yMat(trckIdxInLast,:);

% find those that start and end in the same region
trckIdxInRegion = intersect(trckIdxInFirst,trckIdxInLast); 
xMatInRegion = xMat(trckIdxInRegion,:);
yMatInRegion = yMat(trckIdxInRegion,:); 

xMat3 = xMatInRegion; 
yMat3 = yMatInRegion; 

% find those that start or end in the current region
trckIdxOrRegion = union(trckIdxInFirst,trckIdxInLast); 
xMatOrRegion = xMat(trckIdxOrRegion,:);
yMatOrRegion = yMat(trckIdxOrRegion,:); 

xMat4 = xMatOrRegion; 
yMat4 = yMatOrRegion; 

if onlyTarget == 1 % include those tracks with subTrack end sites in analysis 
    trckIdxIn = trckIdxInLast; 
    xMat1 = xMatFirst;
    yMat1 = yMatFirst;
    xMatIn = xMatLast;
    yMatIn = yMatLast;
    
    forTitle1 = 'subTracks Terminated In SubRoi: Included in Analyis (Red and Pink)';
    forTitle2 = 'subTracks Terminated In Other Regions: Not Included in Analysis (Green) '; 
    forTitle3 = 'subTracks Initiated AND Terminated In SubRoi: Included in Analysis (Pink)';
    filename = 'subTracksTargetedToSubRoiRedAndPink_plusSubTracksLeavingRegionGreen';
    
elseif onlyInitiate == 1 % include those tracks with subtrack initiation sites in analysis 
    trckIdxIn = trckIdxInFirst; 
    xMat1 = xMatLast; 
    yMat1 = yMatLast;
    xMatIn = xMatFirst;
    yMatIn = yMatFirst; 
    
    if onlyNuc ==1 
        x = ' Nucleated ';
        y = ' Nucleated in Other Regions OR Not a Nucleation Event ';
    else 
        x = ' Initiated ';
        y = ' Initiated in Other Regions ';
    end 
    
    forTitle1 = ['subTracks', x,  'In SubRoi: Included in Analyis (Red and Pink)'];
    forTitle2 = ['subTracks', y, ':Not Included in Analysis (Green) ']; 
    forTitle3 = ['subTracks', x,  'AND Terminated In SubRoi: Included in Analysis (Pink)'];
    filename = 'subTracksInitiatedInSubRoiRedAndPink_plusSubTracksInitiatedInOtherRegionsGreen';
    
end 

 

    IN=IN(trckIdxIn,:); 
    % Just exchange zero values for NaN, again IN is a matrix of the same 
    %form as xMat and yMat only replaced by 1 where the tracked coordinate should be 
    %be maintained in the subROI calc
    IN=swapMaskValues(IN,0,NaN); % 1 for the sections of the tracks inside the sub-roi
    OUT=swapMaskValues(IN); % 1 for the sections of the track outside the sub-roi

% plot all member tracks in red on top of the mask
 figure('Visible',visible); 
 imshow(roiMask); 
 hold on; 
 plot(xMat1',yMat1','g'); % plot all tracks targeted/initiated in region

 hold on;
 plot(xMatIn',yMatIn','r') %  plot all tracks targeted/initiate in region- those 
 % that are both will be overwritten in red.  
 hold on; 
 plot(xMat3',yMat3','m')
  
 title({projNameTitle; forTitle1; forTitle2; forTitle3}); 
 % saveas(gcf,[subRoiDir filesep title fileExt])
 saveas(gcf,[subRoiDir filesep filename '.eps'],'psc2');
 
 
 % if option on put all in one folder for easy viewing
 if collectPlots == 1
 [path sub num] = getFilenameBody(subRoiDir);
 sub = [sub num];
 
    if exist([collectedDataPath filesep sub],'dir') == 0
     mkdir([collectedDataPath filesep sub]); 
    end
    
    
    collectedDataPathTracks = [collectedDataPath filesep sub filesep 'subRoiGrowthTracks'];
    if exist(collectedDataPathTracks,'dir') == 0;
        mkdir(collectedDataPathTracks); 
    end 
    
    if exist([collectedDataPathTracks filesep 'tifs'],'dir') == 0
        mkdir([collectedDataPathTracks filesep 'tifs']);
        mkdir([collectedDataPathTracks filesep 'eps']); 
    end 
    
    saveas(gcf,[collectedDataPathTracks filesep 'tifs' filesep projName '.tif']);
    saveas(gcf,[collectedDataPathTracks filesep 'eps' filesep projName '.eps'],'psc2'); 
 end % if collectPlots
 close(gcf); 

%%
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
 figure('Visible',visible); 
 imshow(roiMask); 
 hold on; 
 plot(xMatIn',yMatIn','r')
 title({projNameTitle; 'Tracks Included in SubRoi Regional Analysis (Red)'})
 %saveas(gcf,[subRoiDir filesep 'tracksInSubRoi' fileExt])
 saveas(gcf,[subRoiDir filesep 'tracksInSubRoi.eps'],'psc2');
 close(gcf)
 
 if collectPlots == 1
 [path sub num] = getFilenameBody(subRoiDir);
 sub = [sub num];
 
    if exist([collectedDataPath filesep sub],'dir') == 0
     mkdir([collectedDataPath filesep sub]); 
    end
    
    
    collectedDataPathTracks = [collectedDataPath filesep sub filesep 'subRoiGrowthTracks'];
    if exist(collectedDataPathTracks,'dir') == 0;
        mkdir(collectedDataPathTracks); 
    end 
    
    if exist([collectedDataPath filesep 'tifs'],'dir') == 0
        mkdir([collectedDataPath filesep 'tifs']);
        mkdir([collectedDataPath filesep 'eps']); 
    end 
    projName = [name roi];
    saveas(gcf,[collectedDataPathTracks filesep 'tifs' filesep projName '.tif']);
    saveas(gcf,[collectedDataPathTracks filesep 'eps' filesep projName '.eps'],'psc2'); 
 end % if collectPlots
 
 
 
 
 
 
 
 
 
 
%plot all tracks excluded from analysis
  figure('Visible',visible)
  imshow(roiMask);
  hold on;
  plot(xMatOut',yMatOut','g');
  title({projNameTitle; 'Tracks In SubRoi Region But Excluded From SubRoi Regional Analysis (Green)'}); 
  
  %saveas(gcf,[subRoiDir filesep 'tracksExcludedFromSubRoi' fileExt]);
  saveas(gcf,[subRoiDir filesep 'tracksExcludedFromSubRoi.eps'],'psc2' ); 
  close(gcf)
  
  
  
 


end % if (onlyTarget == 1 || onlyInitiate ==1)


 
%% Extract Data Pause
%if onlyTarget == 1  % Don't calculate the pause data 
 %   else % continue with pause calculations
%Get Pause Data 
idxPause=find(dataMatMergeAll(:,5)==2);
dataMatMergePause=dataMatMergeAll(idxPause,:);

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

lifeSecPause=dataMatMergePause(:,6); % track lifetimes, seconds
%the lifetieme(ie number of frames) 
%inside the ROI is just the sum of each row (sum in the 2nd D)-1
insideSecPause =sourceProjData.secPerFrame.*(sum(INPause,2)-1);

% Get Track IDs for those Growth Sub-Tracks in the Sub-ROI of Interest
if ~isempty(strmatch(timeUnits,'Fraction'))
    trckIdxInPause=intersect(find(insideSecPause./lifeSecPause>=timeVal),inIncludeRegionPause);
    trckIdxOutPause=intersect(find(insideSecPause./lifeSecPause<timeVal & insideSecPause > 0),inIncludeRegionPause);
else
    trckIdxInPause=intersect(find(insideSec>=timeVal),inIncludeRegion);
    trckIdxOutPause=intersect(find(insideSecPause<timeVal & insideSecPause > 0),inIncludeRegion);
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
 figure('Visible',visible); 
 imshow(roiMask); 
 hold on; 
 plot(xMatPause',yMatPause','b')
 title({projNameTitle; 'Pauses In SubRoi'}); 
 %saveas(gcf,[subRoiDir filesep 'PausesInSubRoi' fileExt])
 saveas(gcf,[subRoiDir filesep 'PausesInSubRoi.eps'],'psc2')
 close(gcf)

%end
%% Extract Data Shrinkage

%if onlyTarget == 1 % Don't calculate the shrinkage data
 %   else % Continue with the shrinkage calculations

%Get Shrinkage Data

idxShrink=find(dataMatMergeAll(:,5)==3);
dataMatMergeShrink=dataMatMergeAll(idxShrink,:);

% output is index for merged coordinates
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

lifeSecShrink=dataMatMergeShrink(:,6); % track lifetimes, seconds
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
 figure('Visible',visible); 
 imshow(roiMask); 
 hold on; 
 plot(xMatShrink',yMatShrink','k')
 title({projNameTitle;'Shrinkages In SubRoi'})
 %saveas(gcf,[subRoiDir filesep 'ShrinkagesInSubRoi' fileExt])
 saveas(gcf,[subRoiDir filesep 'ShrinkagesInSubRoi.eps'],'psc2'); 
 close(gcf)
 
%end
%% Write Stats to Output

projData=sourceProjData;
projData.anDir=subRoiDir; % path to sub-roi
projData.imDir=projData.imDir;
 
  xCoordIn = xMatIn.*IN; % calc coordinates before removing beg and end for plotting
    yCoordIn = yMatIn.*IN; % calc coordinates before removing beg and end for plotting
    
    
    xCoordOut = xMatIn.*OUT; 
    yCoordOut = yMatIn.*OUT; 
    
    projData.xCoordInAllTracks = xCoordIn;  % save these coordinates before removing beg and end for plotting
    projData.yCoordInAllTracks = yCoordIn;

    projData.xCoordOutAllTracks= xCoordOut ; 
    projData.yCoordOutAllTracks = yCoordOut; 


if ~isempty(trckIdxIn)
projData.nTracks = length(unique(dataMatMerge(trckIdxIn,1))); % number of
% compound tracks represented (NOTE: depending on the user specifications the 
 %complete compound track may not be in the region)
%projData.stats.percentTracksStartAndEndInRegion = length(trckIdxInRegion)/length(trckIdxIn)*100;  

else 
    projData.nTracks = 0; 
    projData.stats.percentTracksStartAndEndInRegion = 0;    
end 


%Fix some of the old fields in projData
projData.percentFgapsReclass='See projData for whole cell';
projData.percentBgapsReclass='See projData for whole cell';

if isfield(projData,'tracksWithFgap'); 
    projData = rmfield(projData,'tracksWithFgap'); 
end 

if isfield(projData,'tracksWithBgap'); 
    projData = rmfield(projData,'tracksWithBgap'); 
end 

% specifically write those compound tracks with an Fgap or Bgap included in
% region 
if ~isempty(trckIdxInPause)
projData.tracksWithFgapSubRoi = dataMatMergePause(trckIdxInPause,1);
projData.tracksWithBgapSubRoi = dataMatMergeShrink(trckIdxInShrink,1);
else 
    projData.tracksWithFgapSubRoi = 0; 
    projData.tracksWithBgapSubRoi = 0; 
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Median Nearest Neighbor For SubRoi 

 [numSubTracks numFrames] = size(xMatIn);
 if (isempty(xMatIn) || numSubTracks == 1)
     projData.stats.medNNdistWithinFrameMic = NaN; 
 else  
  NNdist=nan(length(xMatIn(:,1))*numFrames,1);
    count=1;
   
    for iFrame=5:numFrames
        
        if ~isempty(xMatIn)
        xCoord = xMatIn(:,iFrame); 
        yCoord = yMatIn(:,iFrame);

            D=createDistanceMatrix([xCoord yCoord],[xCoord yCoord]);
            [sD,idx]=sort(D,2);

            NNdist(count:count+length(xCoord)-1)=sD(:,2);

        end
        count=count+length(xCoord);
    end

    projData.stats.medNNdistWithinFrameMic = nanmedian(NNdist)*projData.pixSizeNm/1000; 
 end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% keep only the coordinates, speeds, etc. corresponding to tracks remaining
% put in xCoord form (row = compoundtrackID, column = Frame)

% Initiate
projData.xCoord=nan(size(sourceProjData.xCoord));
projData.yCoord=nan(size(sourceProjData.yCoord));
% projData.featArea=nan(size(sourceProjData.featArea));
% projData.featInt=nan(size(sourceProjData.featInt));
projData.frame2frameVel_micPerMin=nan(size(sourceProjData.frame2frameVel_micPerMin));
projData.segGapAvgVel_micPerMin=nan(size(sourceProjData.segGapAvgVel_micPerMin));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write the the coords for growth that are located in the subRoi

for iSub=1:length(trckIdxIn)
    k=trckIdxIn(iSub);

    projData.xCoord(dataMatMerge(k,1),dataMatMerge(k,2):dataMatMerge(k,3))=sourceProjData.xCoord(dataMatMerge(k,1),dataMatMerge(k,2):dataMatMerge(k,3));
    projData.yCoord(dataMatMerge(k,1),dataMatMerge(k,2):dataMatMerge(k,3))=sourceProjData.yCoord(dataMatMerge(k,1),dataMatMerge(k,2):dataMatMerge(k,3));

%     projData.featArea(dataMatMerge(k,1),dataMatMerge(k,2):dataMatMerge(k,3))=sourceProjData.featArea(dataMatMerge(k,1),dataMatMerge(k,2):dataMatMerge(k,3));
%     projData.featInt(dataMatMerge(k,1),dataMatMerge(k,2):dataMatMerge(k,3))=sourceProjData.featInt(dataMatMerge(k,1),dataMatMerge(k,2):dataMatMerge(k,3));

    projData.frame2frameVel_micPerMin(dataMatMerge(k,1),dataMatMerge(k,2):dataMatMerge(k,3)-1)=sourceProjData.frame2frameVel_micPerMin(dataMatMerge(k,1),dataMatMerge(k,2):dataMatMerge(k,3)-1);
    projData.segGapAvgVel_micPerMin(dataMatMerge(k,1),dataMatMerge(k,2):dataMatMerge(k,3)-1)=sourceProjData.segGapAvgVel_micPerMin(dataMatMerge(k,1),dataMatMerge(k,2):dataMatMerge(k,3)-1);
end

%%
% Calculate some other variables based on growth coordinates for SubRoi

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


%%
% Create a DataMat of just SubRoi Tracks for Stats

dataGrowthROI = dataMatMerge(trckIdxIn,:);
dataPauseROI = dataMatMergePause(trckIdxInPause,:);
dataShrinkROI = dataMatMergeShrink(trckIdxInShrink,:);

dataTotROI = [dataGrowthROI ; dataPauseROI ; dataShrinkROI];
dataTotROI = sortrows(dataTotROI);

% put here the merged data in frames and pixels 
dataTotROI_frame_pix = dataTotROI; 
dataTotROI_frame_pix(:,6) = dataTotROI_frame_pix(:,6)./sourceProjData.secPerFrame;
dataTotROI_frame_pix(:,7) = dataTotROI_frame_pix(:,7)./sourceProjData.pixSizeNm; 
projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix= dataTotROI_frame_pix; % 

%put here the merged data converted 
projData.mergedDataMatAllSubTracksConverted = dataTotROI;

%Data regarding excluded tracks
%but were excluded based on user specified time requirement 
%dataGrowthExclude = dataMatMerge(trckIdxOut,:);
%dataPauseExclude = dataMatMerge(trckIdxOutPause,:);
%dataShrinkExclude = dataMatMergeShrink(trckIdxOutShrink,:);

%dataTotExclude = [dataGrowthExclude ; dataPauseExclude];
%dataTotExclude = [dataTotExclude ; dataShrinkExclude];
%projData.tracksExcludedFromAnalysis = dataTotExclude;
%%
% Get Lifetimes Specific to Growth SubTracks: Total and Just Inside Region
% to Collect Stats

lifeSec = lifeSec(trckIdxIn); 
insideSec= insideSec(trckIdxIn); 

projData.insideSecAllTracks = insideSec; % save these tracks before removing subtracks 
% at the beginning or end of movie (for plotting)

% calculate speed in before removing the beginning and the end tracks (for
% plots)
    speedIn=nanmean(sqrt(diff(xMatIn.*IN,[],2).^2+diff(yMatIn.*IN,[],2).^2),2);
    speedOut=nanmean(sqrt(diff(xMatIn.*OUT,[],2).^2+diff(yMatIn.*OUT,[],2).^2),2);
    
    projData.speedInMicPerMinAllTracks=pixPerFrame2umPerMin(speedIn,projData.secPerFrame,projData.pixSizeNm);
    projData.speedOutMicPerMinAllTracks=pixPerFrame2umPerMin(speedOut,projData.secPerFrame,projData.pixSizeNm); 

% these values will be used for plotting

%%
% Remove Growth SubTracks (and any preceding fgap/bgap) that begin in frame
% first frame or end in last frame if user specifies
% Adjust lifeSec and insideSec column to include only those sub-tracks included
% in stats

if remBegEnd== 1
[dataTotROI] = plusTipRemBegEnd(dataTotROI, projData);
 projData.removeBegEnd = 'yes'; 
 
 % note for the user the percentage of growths at the start or end that 
 % they are removing
 projData.startOrEnd=~isnan(xMatIn(:,1)) | ~isnan(xMatIn(:,end));
 %projData.percentGrowthAtStartOrEnd=sum(projData.startOrEnd)./projData.stats.nGrowths*100;
 
 % modify lifeSec and insideLifeSec to only include those included in stats

 
lifeSec = lifeSec(~projData.startOrEnd); % remove subtracks at the beg or end of  movie for stats
insideSec = insideSec(~projData.startOrEnd); 



 % modify the coordinate matrix and the lifeSec matrix so can calculate 
 % some other stats from cropped data structure 
 
 xMatIn = xMatIn(~projData.startOrEnd,:); 
 yMatIn = yMatIn(~projData.startOrEnd,:); 
 IN = IN(~projData.startOrEnd,:); 
 OUT = OUT(~projData.startOrEnd,:); 
%  
else
     if isfield(projData,'dataMatCropped_MicPerMin_Sec_Mic')== 1; 
     projData =  rmfield(projData,'dataMatCropped_MicPerMin_Sec_Mic'); % a field in an older version
     end 
     projData.removeBegEnd = 'no';     
end 

%%
% WRITE NEW FIELDS NOT IN SOURCE (WHOLE CELL) PROJDATA.MAT
if projData.nTracks~=0
    pixSizMic=projData.pixSizeNm/1000; % size of a pixel in microns
    pixAreaSqMic=pixSizMic^2; % area of a pixel in square microns
    cellAreaSqMic=sum(wholeCellRoiMask(:)).*pixAreaSqMic;
    % calculate area in square microns of this roi and its
    % percentage of the full cell's area
    subRoiAreaSqMic=sum(roiMask(:)).*pixAreaSqMic;
    percentRoiArea=100*(subRoiAreaSqMic/cellAreaSqMic);
    projData.subRoiAreaSqMic=subRoiAreaSqMic;
    projData.percentRoiArea=percentRoiArea;
   
    % record in the projData a list of the lifetime of each growth subtrack
    % ONLY within the subregion of interest. (the subroi partitioning 
    % allows you to place a long track to a specific subregion but this value 
    % will calculate two local parameters(speed and growth lifetime) 
    % within this region)
     
    projData.lifeSec=lifeSec; % total lifetime (seconds) % only those cacluated in stats
    projData.insideSec=insideSec; % lifetime within sub-roi (seconds)
    projData.percentLifeInside=100*(projData.insideSec./projData.lifeSec); % percent time within sub-roi
    
    % record the median and mean values of the lifetime of the 
    % the subdivided growth events
    
    projData.stats.growth_lifetime_median_INSIDE_REGION = median(projData.insideSec); 
    projData.stats.growth_lifetime_mean_INSIDE_REGION = mean(projData.insideSec); 
    
    speedIn=nanmean(sqrt(diff(xMatIn.*IN,[],2).^2+diff(yMatIn.*IN,[],2).^2),2);
    speedOut=nanmean(sqrt(diff(xMatIn.*OUT,[],2).^2+diff(yMatIn.*OUT,[],2).^2),2);
    
    projData.speedInMicPerMin=pixPerFrame2umPerMin(speedIn,projData.secPerFrame,projData.pixSizeNm);
    projData.speedOutMicPerMin=pixPerFrame2umPerMin(speedOut,projData.secPerFrame,projData.pixSizeNm);
    
    projData.stats.growth_speed_median_INSIDE_REGION = nanmedian(projData.speedInMicPerMin); 
    projData.stats.growth_speed_mean_INSIDE_REGION = nanmean(projData.speedInMicPerMin); 
    
else
    
    projData.lifeSec=NaN;
    projData.insideSec=NaN;
    projData.percentLifeInside=NaN;

    projData.speedInMicPerMin=NaN;
    projData.speedOutMicPerMin=NaN;
   

    projData.startOrEnd=NaN;
    projData.percentGrowthAtStartOrEnd=NaN;
    projData.stats.growth_speed_mean_INSIDE_REGION = NaN; 
    projData.stats.growth_lifetime_median_INSIDE_REGION = NaN; 
    projData.stats.growth_lifetime_mean_INSIDE_REGION = NaN; 
    
end


%%
%%
% Calculate the stats 
% Note the index of one at the end tells the program it is calling this function
% from subRoiAnalysis and to not calculate certain params that require
% intact compound tracks

    [projData,M]=plusTipDynamParam(dataTotROI,projData,0,1);
   



turnHistOn  = 0 ; 
if turnHistOn == 1
% Make Histograms % make option to turn these on
 if projData.nTracks ~= 0 
plusTipMakeHistograms(M,[subRoiDir filesep 'meta' filesep 'histograms']);

%%
%Calculate Polarity Plots for Growth Tracks In SubRoi (uses projData.xCoord/projData.yCoord)

% for all angles assigned to subRegion
xCoord = xMatIn; 
yCoord = yMatIn; 

[anglesFinal,subTrackAnglesIndMean, subTrackAnglesIndStd] = plusTipPlotTrackAngles2(xCoord,yCoord); 

projNameTit = regexprep(projName,'_',' '); 

saveFig1 =figure('Visible',visible);
rose(anglesFinal); % based on all frame to frame angles 
projData.anglesFrame2Frame = anglesFinal; 

title({projNameTit; 'Rose Plot of All Frame-to-Frame Growth Displacement (Converted to Theta)'});
saveas(saveFig1,[subRoiDir filesep 'meta' filesep 'histograms' filesep 'angles_histogram.eps'], 'psc2');


projData.subTrackAnglesIndMean = [dataTotROI(dataTotROI(:,5)==1,1), subTrackAnglesIndMean]; 
projData.subTrackAnglesIndStd = [dataTotROI(dataTotROI(:,5)==1,1),subTrackAnglesIndStd]; 


projData.stats.polarCoordMeanOfAllSubtracks = nanmean(projData.subTrackAnglesIndMean(:,2)); 
projData.stats.polarCoordMedianOfAllSubtracks = nanmedian(projData.subTrackAnglesIndMean(:,2)); 
projData.stats.polarCoordStdOfAllSubtracks = nanstd(projData.subTrackAnglesIndMean(:,2)); 


% for angles only locally in subRegion 
xCoordInside = xMatIn.*IN; 
yCoordInside = yMatIn.*IN;

xCoordOutside = xMatIn.*OUT; 
yCoordOutside = yMatIn.*OUT; 



[anglesFinalInside,subTrackAnglesIndMeanInside] = plusTipPlotTrackAngles2(xCoordInside,yCoordInside); 



projData.subTrackAnglesIndMean_INSIDE_REGION = [dataTotROI(dataTotROI(:,5)==1,1),subTrackAnglesIndMeanInside]; 

projData.stats.polarCoordMeanOfAllSubtracks_INSIDE_REGION= nanmean(projData.subTrackAnglesIndMean_INSIDE_REGION(:,2)); 
projData.stats.polarCoordMedianOfAllSubtracks_INSIDE_REGION = nanmedian(projData.subTrackAnglesIndMean_INSIDE_REGION(:,2)); 
projData.stats.polarCoordStdOfAllSubtracks = nanstd(projData.subTrackAnglesIndMean_INSIDE_REGION(:,2)); 

saveFig2 = figure('Visible',visible); 
rose(anglesFinalInside);
projData.anglesFrame2Frame_INSIDE_REGION = anglesFinalInside; 
title({projNameTit ; 'Rose Plot of All Frame-to-Frame Growth Displacement (Converted to Theta): INSIDE SubRoi Only'});
saveas(saveFig2,[subRoiDir filesep 'meta' filesep 'histograms' filesep 'angles_histogramInside.eps'], 'psc2');
 

collectedDataPathPolarity = [collectedDataPath filesep 'sub' filesep 'Polarity'];
if collectPlots == 1
    if ~exist(collectedDataPathPolarity,'dir')==1;
        mkdir(collectedDataPathPolarity);
    end    
    
    if ~exist([collectedDataPathPolarity filesep 'tifs'],'dir') == 1
        mkdir([collectedDataPathPolarity filesep 'tifs']); 
    end 
    if ~exist([collectedDataPathPolarity filesep 'eps'],'dir') == 1
       
        mkdir([collectedDataPathPolarity filesep 'eps']);
    end 
    saveas(saveFig1,[collectedDataPathPolarity filesep 'tifs' filesep projName 'angles_histogram.tif']); 
    saveas(saveFig2,[collectedDataPathPolarity filesep 'tifs' filesep projName 'angles_histogrameInside.tif']); 
    saveas(saveFig1,[collectedDataPathPolarity filesep 'eps' filesep projName 'angles_histogram.eps'], 'psc2');
    saveas(saveFig2,[collectedDataPathPolarity filesep 'eps' filesep projName 'angles_histogramInside.eps'], 'psc2');
end 
close(saveFig1)
close(saveFig2)



% % Calculate the angle that MT hit the cell boundary
% load([up2 filesep 'contour_normal.mat']);
% projData.contourYX = contourYX;
% projData.normalYX = normalYX;
% clear('contourYX','normalYX');
% 
% projData=plusTipIncidence(projData,xCoord,yCoord);
   


end 
%test = figure; 
 %imshow(roiMask); 
 %hold on; 
 %plot(xCoordInside',yCoordInside','r'); % plot all tracks targeted/initiated in region
 %saveas(test,[collectedDataPath filesep sub filesep 'Polarity' filesep projName 'Inside.eps'], 'psc2');
 %close(test)
%%
end 
% Now write the interpolated coordinates corresponding to pause (necessary
% for plotting functions)

for iSub=1:length(trckIdxInPause)
    k=trckIdxInPause(iSub);

    projData.xCoord(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3))=sourceProjData.xCoord(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3));
    projData.yCoord(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3))=sourceProjData.yCoord(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3));

%     projData.featArea(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3))=sourceProjData.featArea(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3));
%     projData.featInt(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3))=sourceProjData.featInt(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3));

    projData.frame2frameVel_micPerMin(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3)-1)=sourceProjData.frame2frameVel_micPerMin(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3)-1);
    projData.segGapAvgVel_micPerMin(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3)-1)=sourceProjData.segGapAvgVel_micPerMin(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3)-1);
end

%%
%Now write the interpolated coordinates corresponding to shrink (necessary
%for plotting functions)

for iSub=1:length(trckIdxInShrink)
    k=trckIdxInShrink(iSub);

    projData.xCoord(dataMatMergeShrink(k,1),dataMatMergeShrink(k,2):dataMatMergeShrink(k,3))=sourceProjData.xCoord(dataMatMergeShrink(k,1),dataMatMergeShrink(k,2):dataMatMergeShrink(k,3));
    projData.yCoord(dataMatMergeShrink(k,1),dataMatMergeShrink(k,2):dataMatMergeShrink(k,3))=sourceProjData.yCoord(dataMatMergeShrink(k,1),dataMatMergeShrink(k,2):dataMatMergeShrink(k,3));
% 
%     projData.featArea(dataMatMergeShrink(k,1),dataMatMergeShrink(k,2):dataMatMergeShrink(k,3))=sourceProjData.featArea(dataMatMergeShrink(k,1),dataMatMergeShrink(k,2):dataMatMergeShrink(k,3));
%     projData.featInt(dataMatMergeShrink(k,1),dataMatMergeShrink(k,2):dataMatMergeShrink(k,3))=sourceProjData.featInt(dataMatMergeShrink(k,1),dataMatMergeShrink(k,2):dataMatMergeShrink(k,3));
% 
    projData.frame2frameVel_micPerMin(dataMatMergeShrink(k,1),dataMatMergeShrink(k,2):dataMatMergeShrink(k,3)-1)=sourceProjData.frame2frameVel_micPerMin(dataMatMergeShrink(k,1),dataMatMergeShrink(k,2):dataMatMergeShrink(k,3)-1);
    projData.segGapAvgVel_micPerMin(dataMatMergeShrink(k,1),dataMatMergeShrink(k,2):dataMatMergeShrink(k,3)-1)=sourceProjData.segGapAvgVel_micPerMin(dataMatMergeShrink(k,1),dataMatMergeShrink(k,2):dataMatMergeShrink(k,3)-1);
end




% save projData in meta folder
save([subRoiDir filesep 'meta' filesep 'projData'],'projData')
% write out speed/lifetime/displacement distributions into a text file
dlmwrite([subRoiDir filesep 'meta' filesep 'gs_fs_bs_gl_fl_bl_gd_fd_bd.txt'], M, 'precision', 3,'delimiter', '\t','newline', 'pc');

% Write stats results into a text file
statsFile = [subRoiDir filesep 'meta' filesep 'Stats.txt'];
statsData= struct2cell(projData.stats);
statsName = fieldnames(projData.stats);
fid=fopen(statsFile,'w+');
for i=1:numel(statsName)
    fprintf(fid,'%s\t%g\n',statsName{i},statsData{i});
end
fclose(fid);


cd(homeDir)
