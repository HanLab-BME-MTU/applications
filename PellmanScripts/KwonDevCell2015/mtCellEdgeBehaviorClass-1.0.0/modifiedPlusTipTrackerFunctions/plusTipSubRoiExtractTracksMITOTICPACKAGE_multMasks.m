function [projData,M]=plusTipSubRoiExtractTracksMITOTICPACKAGE_multMasks(MD,subRoiDir,varargin)
% this fn is called by plusTipSubRoiTool
%
% UPDATES 03/2011: This Function was MODIFIED from the 1st release of
% plusTipTracker.  
% Now partitions pause and shrinkage information from the
% given subRegions.  Also fixes small bug in very code where the
% Fraction option was ALWAYS bi-passed, even if it was
% selected in the GUI, due to a small capitalization error.
% 
% Furthermore, it has internal option to partition subtracks based on if
% their start site (onlyInitiate option) or end site (onlyTarget)
% is located within the region of interest, or those tracks specifically 
% crossing over the borders of the region of interest 
%  Maria Bagonis (MB) 
%
% UPDATES 04/2012 Bugs fixes/cleaning and more extensive options for 
% subtrack extraction. 
% Maria Bagonis (MB) 
% UPDATES 06/2016 (MB)
% Includes extraction from multiple masks using the multMask1 extractType 
% Includes an option to calculate dwell time from small masks around the 
% terminal point. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% 
% REQUIRED:                  
% subRoiDir: the pathname pointing to the subRoi directory  
%                The subRoiDirectory needs to include 
%                the roiMask.tif (typcially created by the subRoiTool)
%                
% OPTIONAL: (made optional to not change previous input struct)
%
%  You can filter the growth subtracks for analysis based on time in region
%  There are two methods by which one can threshold
% 'fraction': the fraction of the total lifetime of the subtrack or 
% 'seconds': the absolute time in region 
% 
% Default = 50 % of the total lifetime of the subtrack needs to be in the subregion 
% to be included in statistics. 
%   
%
% timeUnitsGrowth: The metric by which you would you like to partition your growth subtrack data: 
%                  'fraction' (of total growth time) or 'second' (absolute value)
%                  Default = 'fraction'
%              
% timeValGrowth: The minimum time requirement in region which a growth
%                subtrack must meet to be included in the analysis
%                For 'fraction' needs to be between 0-1,
%                For 'second' give absolute value in seconds
%                Default = 0.5; 
%
%                %% BELOW NOT YET OPTIONS IN GUI %% 
% You can filter the inferred pause subtracks via a time requirement
% (though this is not really necessary, especially for pauses, as these events do not typically span a region)
%
% Default = 100 % of the total lifetime of the subtrack needs to be in the
% subregion to be included in the statistics
% 
% timeUnitsPause: The metric by which you would like to partition your
%                 pause subtrack data.
%                 'fraction' (of total pause time) or 'second' (absolute value)
%                 Default = 'fraction'
% 
% timeValPause:  The minimum time requirement in region 
%                which a pause subtrack must meet to be included in the
%                analysis
%                For timeUnitsPause = 'fraction': needs to be between 0-1,
%                For timeUnitsPause = 'seconds' give absolute value in seconds
%                Default = 1; 
%
% You can filter the inferred shrink subtracks via a time requirement
% (though this is not really necessary, pauses do not typically span a region)
% Default = 100 % of the total lifetime of the subtrack needs to be in the subregion 
% to be included in statistics. 
%
% timeUnitsShrink: The metric by which you would like to partition your
%                 shrink subtrack data.
%                 'fraction' (of total shrink time) or 'second' (absolute value)
%                 Default = 'fraction'
% 
% 
% timeValShrink: The minimum time requirement in region 
%                which an inferred shrinkage event must meet to be included in the
%                analysis
%                For timeUnitsPause = 'fraction': needs to be between 0-1,
%                For timeUnitsPause = 'seconds' give absolute value in seconds
%                Default = 1; 
% PARAMS: 
% extractTypeGrowth: 
%        'nonDirectional': either end or beginning of subtrack in region 
%        'subTrackEnd': end point of subtrack in region                       
%        'subTrackStart': start point of subtrack in region
%        'nucleation': only nucleation sites
%        'boundCrossIn': only those subtracks crossing into the region  
%        'boundCrossOut': only those subtracks crossing out of the region 
% 
%     
%%Input check
ip = inputParser;
ip.addRequired('subRoiDir',@(x)ischar(x) || isempty(x));

% make time requirements optional for now so don't have to change the
% input structure: 
ip.addOptional('timeUnitsGrowth','fraction', @(x)strcmpi(x,'fraction') ||strcmpi(x,'seconds')); 
ip.addOptional('timeValGrowth', 0.5 , @isscalar)

ip.addParamValue('timeUnitsPause','fraction',@ischar); 
ip.addParamValue('timeValPause',1,@isscalar); 
ip.addParamValue('timeUnitsShrink','fraction',@ischar); 
ip.addParamValue('timeValShrink',1',@isscalar); 



% Extract types (DEFAULT is nonDirectional any point within region)
ip.addParamValue('extractType','nonDirectional',@ischar) %  you can change the default extractType here
ip.addParamValue('remBegEnd', 1, @isscalar); 
ip.addParamValue('turnFiguresOn',1,@isscalar); 
ip.addParamValue('collectPlots',0,@isscalar); 
ip.addParamValue('turnHistOn',1,@isscalar); 
ip.addParamValue('mkPolPlot',0,@isscalar); 

% Added for Mitotic 
ip.addParamValue('bipolarMask',false); % flag to read in the opposite side 
                                       % so can filter tracks incoming from
                                       % this side 
ip.addParamValue('performDwellCalcs',true);                                       
ip.addParamValue('dwellRadius', 2); %  specifies the radius of the disk shaped structuring element used 
% for the dilation around the terminal point of the MT subTrack. 

ip.parse(subRoiDir,varargin{:});

timeUnitsGrowth = ip.Results.timeUnitsGrowth;
timeValGrowth = ip.Results.timeValGrowth;
timeUnitsPause = ip.Results.timeUnitsPause; 
timeValPause = ip.Results.timeValPause; 
timeUnitsShrink = ip.Results.timeUnitsShrink; 
timeValShrink = ip.Results.timeValShrink; 
remBegEnd = ip.Results.remBegEnd; 
extractType = ip.Results.extractType;
turnFiguresOn = ip.Results.turnFiguresOn;
collectPlots = ip.Results.collectPlots; 
turnHistOn = ip.Results.turnHistOn; 
mkPolPlot = ip.Results.mkPolPlot; 


if isempty(subRoiDir)
    error('please enter a valid subroi file name');
end

if isempty(strfind(subRoiDir,'sub')) % maybe change
    warning('the file you are analyzing may not be a subRoi');
end

if ~ismember(lower(timeUnitsGrowth),{'fraction','seconds'})
    error('plusTipSubRoiTool: timeUnits must be fraction or seconds')
end
if isempty(timeValGrowth)
    error('plusTipSubRoiTool: time input missing')
end
if isempty(strcmpi(timeUnitsGrowth,'fraction')) && ~(timeValGrowth >0 && timeValGrowth<=1)
    error('plusTipSubRoiTool: timeUnitsGrowth is fraction, timeVal must be in 0-1')
end

%%
%% Extract Type
onlyTarget = 0 ;
onlyInitiate = 0;
onlyNuc = 0;
boundCrossIn = 0 ;
boundCrossOut = 0;
nonDirectional =0; 
switch extractType
    case 'subTrackEnd'
        onlyTarget = 1;
    case 'subTrackStart'
        onlyInitiate = 1;
    case 'nucleation'
        onlyNuc = 1;
        onlyInitiate = 1;
    case 'boundCrossIn'
        boundCrossIn = 1;
    case 'boundCrossOut' 
        boundCrossOut =1; 
    case 'nonDirectional'
        nonDirectional =1; 
        multMasks =1;
    case 'multMasks1' 
        multMasks = 1 ;
        inOnly =0   ;
    case 'multMask2' 
        multMasks =1 ; 
        inOnly =1; 
end

if turnFiguresOn == 1
    visible = 'On';
else
    visible= 'Off';
end


%% Load Necessary Files AND Initiate
load(MD.processes_{1,3}.outFilePaths_{1});
sourceProjData = projData; 
%% Read in Multiple Masks if applicable : 
% NOTE: this is the main difference compared to the 
% plusTipSubRoiExtractTracksFORPACKAGE in the current plusTipTracker/UTrack release

if multMasks == 1
   % maskDir = [upThree filesep 'SegmentationPackage' filesep 'masks' filesep 'masks_for_channel_1']; 
maskDir = [subRoiDir filesep 'masks']; 
masks = mitoticSearchFiles('.tif',[],maskDir,0);
[upOne, ~, sub] = mitoticUpDirectory(subRoiDir,1); 

% For bipolar mitotic masks based on the pole axis
% we delete MT tracks originating from the opposite side of the cell
if ip.Results.bipolarMask
    if strcmp(sub,'1')
        otherSide ='2';
    else
        otherSide ='1';
    end
    masksOtherSideDir = [upOne filesep 'sub_' otherSide filesep 'masks' ];
    masksOtherSide = mitoticSearchFiles('.tif',[], masksOtherSideDir ,0);
end

end % if multMasks

roiMask=logical(imread([masks{1,2} filesep masks{1,1}])); % sometimes mask not saved correctly make sure logical 
wholeCellRoiMask=logical(imread([masks{1,2} filesep masks{1,1}]));
[imL,imW]=size(roiMask);

% if there is an exclude mask, use it; otherwise, use the inverse of the
% whole cell mask and write that into the sub_x folder
if exist([subRoiDir filesep 'excludeMask.tif'],'file')~=0
    excludeMask=imread([subRoiDir filesep 'excludeMask.tif']);
else
    excludeMask=~wholeCellRoiMask;
    imwrite(excludeMask,[subRoiDir filesep 'excludeMask.tif']);
end

%combine the exclude mask and the roiMask
roiMask = roiMask.*~excludeMask; 

%%%%%%%%% Initiate Filenames For Figure Making %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [up3 roi numSub] = getFilenameBody(upTwo);
% [collectedDataPath folder numProj] = getFilenameBody(up3);
% collectedDataPath = [collectedDataPath filesep 'collectedSubRoiPlots'];
% roi = [roi numSub];
% name = [folder numProj];
% projName = [name ' ' roi];
% projNameTitle = regexprep(projName,'_',' ');
projNameTitle = []; 
% delete old figures
delete([subRoiDir filesep 'extract_plots' filesep '*.fig']);
delete([subRoiDir filesep 'extract_plots' filesep '*.eps']); 

% make directories if they don't exist 
if ~isdir([subRoiDir filesep 'meta'])
    mkdir([subRoiDir filesep 'meta'])
end 

if ~isdir([subRoiDir filesep 'feat'])
    mkdir([subRoiDir filesep 'feat'])

    % copy over movieInfo file
anDir = upTwo; 
 [a,b,c]=copyfile([anDir filesep 'feat' filesep 'movieInfo.mat'],[subRoiDir filesep 'feat' filesep 'movieInfo1.mat']);
 [a,b,c]=movefile([subRoiDir filesep 'feat' filesep 'movieInfo1.mat'],[subRoiDir filesep 'feat' filesep 'movieInfo.mat']);

end 

if ~isdir([subRoiDir filesep 'extract_plots'])
    mkdir([subRoiDir filesep 'extract_plots']); 
end 
%% Collect Data to Be Subdivided

% get all the original merged tracks and convert frames to seconds and
% pixels to microns

% Reads in merged output

dataMatMergeAll = sourceProjData.mergedDataMatAllSubTracksConverted;
dataMatNoBegEnd = sourceProjData.dataMat_FOR_STATS; 
nTracks = length(unique(dataMatNoBegEnd(:,1))); 

% get growth track indices only and coordinates for those tracks only
idx=find(dataMatMergeAll(:,5)==1);
dataMatMerge=dataMatMergeAll(idx,:); % just growth subtracks

%Note: xMat: x-coordinate of detected particle: row # = subtrack ID, col.# = frame number
% 1 indicates that it will use merged data
[xMat,yMat]=plusTipGetSubtrackCoords(sourceProjData,idx,1);

%%
%%%% DIVIDE GROWTH SUBTRACK COORDINATES TO SUBROI AND NON-SUBROI %%%%%%%%%%

%Convert from xy indices to pixel index
%pixIdx maintains structure of xMat/yMat ie row = subTrack ID
% and column = frame number.  Instead of x/y coordinates now point
% denoted by pixel number in a imL by imW image.
pixIdxNaN=sub2ind([imL,imW],ceil(yMat-.5),ceil(xMat-.5));

pixIdx = pixIdxNaN;
% assume that the first pixel in the image will be a zero
% but make sure this is so by making the first roiMask pixel = 0
pixIdx(isnan(pixIdx))=1; % set nan (no coord in frame) to 1. 
roiMask(1,1)=0; % anything not in track will automatically equal zero
% as roiMask(1) = 0 

% find which of the TRACK indices are in the roiMask
% IN is a matrix of row = subtrack ID and col = frame number
% places a 1 in those frames where the point is inside the
% roiMask and zero where it is not
% will later truncate this based on partitioning of data

IN=roiMask(pixIdx); % here is where you would need to apply a different 
                    % mask if you were going to update frame-to-frame
                    % could do pixIdx each frame/roiMask each frame then 
                    % segregate
                    INOrig = IN; 
lifeSec=dataMatMerge(:,6); % track lifetimes, seconds
% %the lifetime(ie number of frames)
% Find life inside the subRoi: these values will be used for temporal
% criteria below
%inside the ROI is just the sum of each row (sum in the 2nd D)-1
insideSec=sourceProjData.secPerFrame.*(sum(IN,2)-1); %

IN=swapMaskValues(IN,0,NaN); % 1 for the sections of the tracks inside the sub-roi

if ~isempty(IN(~isnan(IN)))  == 1 
OUT=swapMaskValues(IN); % exchange 1 and NaN 
else 
    
    OUT = ones(size(IN));
end 



%% START PARTITIONING TRACKS BASED ON USER SPECIFICATION   %%%%%%

    %%% GET DIRECTIONAL INFORMATION %%%% 
    
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
  
    % get outside mask
            wholeMasks = mitoticSearchFiles('.tif',[],[MD.outputDirectory_ filesep 'masks'], 0); 
    
    if multMasks == 1 ;
        % don't count the last frame as not good for transition or lifetime
        % measurements
        for iFrame = 1:length(masks)-1 
            display(['Extracting MT Tracks for Frame' num2str(iFrame)])
            roiMask = logical(imread([masks{iFrame,2} filesep masks{iFrame,1}]));% load mask correspoinding to just the current frame
            
            % 
            if ip.Results.bipolarMask
                roiMaskOtherSide = logical(imread([ masksOtherSide{iFrame,2} filesep masksOtherSide{iFrame,1}]));
            else
                roiMaskOtherSide = zeros(size(roiMask));
            end
            
            fileRoiMaskOutside = [wholeMasks{iFrame,2} filesep wholeMasks{iFrame,1}]; 
            
            roiMaskOutside = ~double(imread(fileRoiMaskOutside)); 
            
            figure('visible','off')
            
            imshow(roiMask,[]); 
            hold on 
            
            endIn_iFrame = lastPtFrIdx==iFrame;
            % get the logical indexing of those tracks the start and end 
            % in the region
            endInRegion_iFrame = roiMask(pixIdxLast);
            beginInRegion_iFrame = roiMask(pixIdxFirst);
            % get the logical indexing of those tracks that start on other 
            % side
            startInOtherSide_iFrame = roiMaskOtherSide(pixIdxFirst); 
            
             % find tracks greater than 5 frames in length 
            framesTrackTotal = sum(~isnan(pixIdxNaN),2);
          
            greaterThan5 = framesTrackTotal>0; % previously only greater than 5 
           
            outsideStart = roiMaskOutside(pixIdxFirst); 
            
         
            % 1) find tracks that end iFrame
            % 2) find tracks that end in mask current frame
            
            % 3) take out the tracks that likewise start in the region
            % (cross-over only) 
            % 4) take out tracks less than 5 frames (optional can be set)
            % 5) take out tracks that start in other subRegion (optional
            % can be set) 
            % 6) take out tracks that begin outside (means your mask was bad) 
            % 7) take out those those tracks  that cross-over multiple
            % times ( have only 1 cross-over) -  (likely a lateral moving
            % track)
            if inOnly ==1
                idx_subTrackToKeep= find( endIn_iFrame & endInRegion_iFrame & beginInRegion_iFrame & greaterThan5);
            else
                idx_subTrackToKeep= find( endIn_iFrame & endInRegion_iFrame & ~beginInRegion_iFrame &  ...
                    ~startInOtherSide_iFrame & greaterThan5 & ~ outsideStart);
                
            end
            
            % create the "IN" logical matrix for each track
            % place a 1 if the xtrack is IN the region and a 0 if the
            % track is out of the region
            xMat_iFrame = xMat(idx_subTrackToKeep,:);
            yMat_iFrame = yMat(idx_subTrackToKeep,:);
           
            % get the coordinates in the region
            pixIdxNaN_iFrame=sub2ind([imL,imW],ceil(yMat_iFrame-.5),ceil(xMat_iFrame-.5));
            
            pixIdx_iFrame = pixIdxNaN_iFrame;
            
            % assume that the first pixel in the image will be a zero
            % but make sure this is so by making the first roiMask pixel = 0
            pixIdx_iFrame(isnan(pixIdx_iFrame))=1; % set nan (no coord in frame) to 1.
            roiMask(1,1)=0; % anything not in track will automatically equal zero
            % as roiMask(1) = 0
            
            % find which of the TRACK indices are in the roiMask
            % IN is a matrix of row = subtrack ID and col = frame number
            % places a 1 in those frames where the point is inside the
            % roiMask and zero where it is not
            % will later truncate this based on partitioning of data
            
            IN_iFrame=roiMask(pixIdx_iFrame);
            % filter out tracks going in and out 
            numTrans = sum(abs(diff(IN_iFrame,1,2)),2);
            idx_subTrackToKeepFinal = idx_subTrackToKeep(numTrans<=2); % filter out tracks that go in and out of region (want clean transitions)
            IN_iFrame = IN_iFrame(numTrans<=2,:); 
            
            idx_subTrackToKeep_ByFrame{iFrame}  = idx_subTrackToKeepFinal;
            
            xMat_iFrame = xMat(idx_subTrackToKeepFinal,:);
            yMat_iFrame = yMat(idx_subTrackToKeepFinal,:);
            
            % plot sanity checks 
             if ~isempty(xMat_iFrame) 
                arrayfun(@(i) plot(xMat_iFrame(i,:),yMat_iFrame(i,:),'r'),1:size(xMat_iFrame,1)); 
                if ~isdir([subRoiDir filesep 'extract_plots' filesep 'multMaskExtract']); 
                    mkdir([subRoiDir filesep 'extract_plots' filesep 'multMaskExtract']);
                end 
              
                saveas(gcf,[subRoiDir filesep 'extract_plots' filesep 'multMaskExtract' filesep num2str(iFrame,' %03d') '.tif']);   
                close gcf 
            end 
            
            insideSec_ByFrame{iFrame}=sourceProjData.secPerFrame.*(sum(IN_iFrame,2)-1);
            
            IN_eachFrame{iFrame}=swapMaskValues(IN_iFrame,0,NaN);
            OUT_eachFrame{iFrame} = swapMaskValues(IN_eachFrame{iFrame});
           
            xMatIn_ByFrame{iFrame}  = xMat_iFrame.*IN_eachFrame{iFrame} ; % save these in a cell as well 
            yMatIn_ByFrame{iFrame} = yMat_iFrame.*IN_eachFrame{iFrame};
          
            xMatOut_ByFrame{iFrame} = xMat_iFrame.*OUT_eachFrame{iFrame}; 
            yMatOut_ByFrame{iFrame} = yMat_iFrame.*OUT_eachFrame{iFrame}; 
            
            speedIn_iFrame =nanmean(sqrt(diff(xMatIn_ByFrame{iFrame}.*IN_eachFrame{iFrame},[],2).^2+diff(yMatIn_ByFrame{iFrame}.*IN_eachFrame{iFrame},[],2).^2),2);
            speedInByFrame{iFrame} = pixPerFrame2umPerMin(speedIn_iFrame,projData.secPerFrame,projData.pixSizeNm);
            
            dispByFrame{iFrame} = speedInByFrame{iFrame}.*insideSec_ByFrame{iFrame}/60; % in um
            velsByFrame{iFrame} = speedInByFrame{iFrame}; % velocity (um/min)
            
            
        end % iFrame
        
    
        
        IN_MultMasks = vertcat(IN_eachFrame{:}); 
        OUT_MultMasks = vertcat(OUT_eachFrame{:}); 
        % need displacement each frame and coordinates each frame 
        projData.dispByFrame = dispByFrame; % will used this to decide whether to partition as lateral or not 
        projData.velsByFrame = speedInByFrame;  
        
        
        projData.IN_eachFrame = IN_eachFrame; 
        projData.OUT_eachFrame = OUT_eachFrame; 
        projData.xMatIn_ByFrame = xMatIn_ByFrame;  
        projData.yMatIn_ByFrame = yMatIn_ByFrame;
        projData.trckIdxInMultMasks_ByFrame = idx_subTrackToKeep_ByFrame; 
        projData.trckIdxInMultMasks = vertcat(idx_subTrackToKeep_ByFrame{:}); 
        projData.xMatOut_ByFrame = xMatOut_ByFrame; 
        projData.yMatOut_ByFrame = yMatOut_ByFrame; 
        
             
        else 
    %%
    trckIdxInRegionPreTempFilt = intersect(inIncludeRegionFirst,inIncludeRegionLast);
    trckIdxBoundCrossInPreTempFilt = setdiff(inIncludeRegionLast,trckIdxInRegionPreTempFilt);
    trckIdxBoundCrossOutPreTempFilt = setdiff(inIncludeRegionFirst,trckIdxInRegionPreTempFilt);
    
   
    xMatBoundCrossInPreTempFilt = xMat(trckIdxBoundCrossInPreTempFilt,:); 
    yMatBoundCrossInPreTempFilt = yMat(trckIdxBoundCrossInPreTempFilt,:); 
    
    xMatBoundCrossOutPreTempFilt = xMat(trckIdxBoundCrossOutPreTempFilt,:); 
    yMatBoundCrossOutPreTempFilt = yMat(trckIdxBoundCrossOutPreTempFilt,:); 
    
    xMatInRegionPreTempFilt = xMat(trckIdxInRegionPreTempFilt,:); 
    yMatInRegionPreTempFilt = yMat(trckIdxInRegionPreTempFilt,:); 
    
     % if nucleation only: just make the first subtrack specification
    % more specific (ie not a start linked to an fgap or bgap)
    % a mark of 1 in the 8th column of dataMatMerge indicates that the growth subtrack 
    % was a nucleation event
    if onlyNuc == 1
        nucIdx = find(dataMatMerge(:,8) == 1);
        inIncludeRegionFirst = intersect(inIncludeRegionFirst,nucIdx);
       
    end
   
    
    
    %%% APPLY TEMPORAL FILTERS TO GROWTH SUBTRACKS %%%
    
    if  strcmpi(timeUnitsGrowth,'fraction') ==1
        trckIdxInFirst=intersect(find(insideSec./lifeSec>=timeValGrowth),inIncludeRegionFirst);
    else
        trckIdxInFirst=intersect(find(insideSec>=timeValGrowth),inIncludeRegionFirst);
    end

    xMatFirst = xMat(trckIdxInFirst,:);
    yMatFirst = yMat(trckIdxInFirst,:);
    
    if strcmpi(timeUnitsGrowth,'fraction') == 1
        trckIdxInLast= intersect(find(insideSec./lifeSec>= timeValGrowth),inIncludeRegionLast);
    else 
        trckIdxInLast = intersect(find(insideSec>=timeValGrowth),inIncludeRegionLast); 
    end 
    
    xMatLast = xMat(trckIdxInLast,:);
    yMatLast = yMat(trckIdxInLast,:);
    
    if strcmpi(timeUnitsGrowth,'fraction') == 1 
        trckIdxNonDir = find(insideSec./lifeSec>=timeValGrowth); 
    else 
        trckIdxNonDir = find(insideSec>=timeValGrowth); 
    end 
    
    %%% DIVIDE BORDER CROSSING SUBTRACKS THAT STAY IN REGION
    trckIdxInRegion = intersect(trckIdxInFirst,trckIdxInLast);
    trckIdxBoundCrossIn = setdiff(trckIdxInLast,trckIdxInRegion);
    trckIdxBoundCrossOut = setdiff(trckIdxInFirst,trckIdxInRegion);
    
   
    
    xMatBoundCrossIn = xMat(trckIdxBoundCrossIn,:);
    yMatBoundCrossIn = yMat(trckIdxBoundCrossIn,:);
    

    xMatBoundCrossOut = xMat(trckIdxBoundCrossOut,:);
    yMatBoundCrossOut = yMat(trckIdxBoundCrossOut,:);
    
    xMatInRegion = xMat(trckIdxInRegion,:);
    yMatInRegion = yMat(trckIdxInRegion,:);

    end %  if multMasks
     
%% PLOTS %%%%
     if multMasks == 0                    
  if strcmpi(timeUnitsGrowth,'fraction'); 
            x = '';
        else 
            x = ' sec'; 
  end 
    forTitleTime = ['LIFETIME CRITERIA: ' upper(timeUnitsGrowth) ' GREATER OR EQUAL TO ' num2str(timeValGrowth) x]; 
  

% First Make an Overview Plot of All Tracks in Region So User Can Determine
% the Tracks of Interest

figure('Visible',visible);
colormap('gray'); 
imagesc(roiMask); 
axis 'off'; 
hold on; 

plot(xMatBoundCrossInPreTempFilt',yMatBoundCrossInPreTempFilt','b');
plot(xMatBoundCrossOutPreTempFilt',yMatBoundCrossOutPreTempFilt','g'); 
plot(xMatInRegionPreTempFilt',yMatInRegionPreTempFilt','m'); 
 
 forTitle1 = 'Blue: Growth SubTracks Crossing Into SubRoi ';
 forTitle2 = 'Green: Growth SubTracks Crossing Out of SubRoi';
 forTitle3 = 'Magenta: Growth SubTracks Initiated AND Terminated In SubRoi';
 forTitle4 = 'Yellow: Growth SubTracks Passing Through Region'; 

title({'All GrowthSubTracks In Region (Before Temporal/Directional Filtering) ColorCoded by Direction of Growth';...
    forTitle1; forTitle2; forTitle3;forTitle4}); 
filename = 'growthSubTrack_RegionOverview_BeforeTemporalAndDirectionalFiltering';
saveas(gcf,[subRoiDir filesep 'extract_plots' filesep filename '.eps'],'psc2');
    
%%  COLLECT GROWTH SUBTRACK PLOTS IF USER DESIRES
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
        %mkdir([collectedDataPathTracks  filesep 'eps']);
    end
    
    if ~isdir([collectedDataPathTracks filesep 'eps']); 
        mkdir([collectedDataPathTracks filesep 'eps']); 
    end 
    
    projName = [name roi];
    saveas(gcf,[collectedDataPathTracks filesep 'tifs' filesep projName '.tif']);
    saveas(gcf,[collectedDataPathTracks filesep 'eps' filesep projName '.eps'],'psc2');
end % if collectPlots


close(gcf)

     end % if multMask == 0 
%%%% Make plots of tracks after user specified filtering %%%% 


if( onlyTarget==1 || onlyInitiate ==1)  % historical reasons these are clustered together
    % their output plots are designed similarly
    
    % TARGETING AND INITIATION
    if  onlyTarget == 1  % include those tracks with subTrack end sites in analysis
        trckIdxIn = trckIdxInLast;
        
        xMatIn = xMatLast;
        yMatIn = yMatLast;
        
       
        forTitle1 = 'Growth SubTracks Terminated In SubRoi (Red)';
      
      
        filename = 'PartitionOption_GrowthSubTracksTerminatingInSubRoi';
        
    elseif onlyInitiate == 1 % include those tracks with subtrack initiation sites in analysis
        trckIdxIn = trckIdxInFirst;
       
        xMatIn = xMatFirst;
        yMatIn = yMatFirst;
        
        % Again nucleation will have more stringent criteria for
        % for initiation ( performed above)
        % needs different labels
        if onlyNuc ==1
            x = ' Nucleated ';
           
        else
            x = ' Initiated ';
         
        end
        
        forTitle1 = ['Growth SubTracks', x,  'In SubRoi (Red)'];
        filename = 'PartitionOption_GrowthSubTracksInitiatedInSubRoi';
        
    end % if only target
    
    % MAKE PLOTS
    % plot all member tracks in red on top of the mask
    figure('Visible',visible);
    imagesc(roiMask);
    colormap('gray'); 
    axis off; 
    hold on;
    plot(xMatIn',yMatIn','r') %  plot all tracks targeted/initiate in region- those
    % that are both will be overwritten in red.
   
    title({projNameTitle; forTitle1; forTitleTime});
    saveas(gcf,[subRoiDir filesep 'extract_plots' filesep filename '.eps'],'psc2');
    close(gcf)
    
    
    % BOUNDARY CROSSING ONLY
elseif boundCrossIn == 1 % boundary crossing in
    trckIdxIn = trckIdxBoundCrossIn;
    xMatIn = xMatBoundCrossIn;
    yMatIn = yMatBoundCrossIn;
    
   
    
    figure('Visible',visible);
    colormap('gray'); 
    imagesc(roiMask);
    axis off; 
    hold on;
    plot(xMatIn',yMatIn','r') 
    title({projNameTitle; 'Growth SubTracks Crossing Into Region (red): Included in Analysis'; forTitleTime});
    saveas(gcf,[subRoiDir filesep 'extract_plots' filesep 'PartitionOption_GrowthSubTracksCrossingIntoRegion.eps'],'psc2');
    close(gcf)
    
elseif boundCrossOut == 1
    trckIdxIn = trckIdxBoundCrossOut; 
    xMatIn  = xMatBoundCrossOut; 
    yMatIn = yMatBoundCrossOut; 
    
    figure('Visible',visible); 
    colormap('gray') 
    imagesc(roiMask);
    axis off;
    hold on; 
    plot(xMatIn',yMatIn','r'); 
   
    title({projNameTitle; 'Growth SubTracks Crossing Out of Region (red): Included In Analysis';...
    forTitleTime}); 
    saveas(gcf,[subRoiDir filesep 'extract_plots' filesep 'PartitionOption_GrowthSubTracksCrossingOutOfRegion.eps'],'psc2'); 
    close(gcf)
    % NONDIRECTIONAL
elseif nonDirectional == 1   % criteria non-directional but has a time constraint
    
    % get those tracks located within the subregion based with a given
    % time criteria (but NO requirement for whether the track starts or
    % ends in the ROI of interest)
   
   trckIdxIn = trckIdxNonDir; 
   xMatIn = xMat(trckIdxIn,:); 
   yMatIn= yMat(trckIdxIn,:); 
   
   
  
    
    % plot all member tracks in red on top of the mask
    figure('Visible',visible);
    imagesc(roiMask);
    colormap('gray');
    axis off; 
    hold on;
    plot(xMatIn',yMatIn','r')
    title({projNameTitle; 'Growth SubTracks Included in SubRoi Regional Analysis (Red): No Direction Criteria';forTitleTime})
    %saveas(gcf,[subRoiDir filesep 'tracksInSubRoi' fileExt])
    saveas(gcf,[subRoiDir filesep 'extract_plots' filesep 'PartitionOption_GrowthSubTracksInSubRoi.eps'],'psc2');
    close(gcf)
    
    %plot all tracks excluded from analysis
%     figure('Visible',visible)
%     imshow(roiMask);
%     hold on;
%     plot(xMatOut',yMatOut','g');
%     title({projNameTitle; 'Tracks In SubRoi Region But Excluded From SubRoi Regional Analysis (Green)'});
    
%     %saveas(gcf,[subRoiDir filesep 'tracksExcludedFromSubRoi' fileExt]);
%     saveas(gcf,[subRoiDir filesep 'tracksExcludedFromSubRoi.eps'],'psc2' );
%     close(gcf)
    
elseif multMasks == 1

 xMatIn =  xMat(projData.trckIdxInMultMasks,:); 
 yMatIn =  yMat(projData.trckIdxInMultMasks,:);
 trckIdxIn = projData.trckIdxInMultMasks; 
 x = '?'; 
% no plots 

end %(if onlyTarget || onlyInitiate)
%% Added PerformDwellCalculationsOption

if ip.Results.performDwellCalcs;
    if ~ isempty(trckIdxIn)
    dwellMaskCum = zeros(size(roiMask)); 
    pixIdxBoundCrossInNaN =  pixIdxNaN(trckIdxIn,:) ; % all set up so that can calc which ones inside mask
    pixIdxBoundCrossIn = pixIdx(trckIdxIn,:);
    for iSubTrack = 1: length(xMatIn(:,1)) % just set entire thing up as a loop for now...
        % get each endpoint
        
        lastPtFrIdxForDwell= find(~isnan(pixIdxBoundCrossInNaN(iSubTrack,:)),1,'last');
        pixIdxForMask =  pixIdxBoundCrossIn(iSubTrack,lastPtFrIdxForDwell);
        endPtMask = zeros(size(roiMask)); % reinitiate mask
        endPtMask(pixIdxForMask) =1;
        endPtDilate =  imdilate(endPtMask,strel('disk',ip.Results.dwellRadius)); % create a mask dil
        endPtDilate(1,1) = 0;
        dwellMaskCum = dwellMaskCum|endPtDilate; 
        dwellMasks(:,:,iSubTrack) = endPtDilate;% save the masks for plotting
        INDwell = endPtDilate(pixIdxBoundCrossIn(iSubTrack,:));% of all the pixels which ones are in region
        % find first and last 
     
        dwell(iSubTrack)=sourceProjData.secPerFrame.*(sum(INDwell,2)-1);
        INDwell=swapMaskValues(INDwell,0,NaN);
        xCoordsDwell(iSubTrack,:) =  xMatIn(iSubTrack,:).*INDwell;
        yCoordsDwell(iSubTrack,:) = yMatIn(iSubTrack,:).*INDwell;
        
        if ~isempty(INDwell(~isnan(INDwell)))  == 1
            OUTDwell=swapMaskValues(INDwell); % exchange 1 and NaN
        else
            
            OUTDwell = ones(size(INDwell));
        end
        
        xCoordsDwellOut(iSubTrack,:) = xMatIn(iSubTrack,:).*OUTDwell;
        yCoordsDwellOut(iSubTrack,:) = yMatIn(iSubTrack,:).*OUTDwell;
        
    end 
        else 
            xCoordsDwellOut = NaN; 
            yCoordsDwellOut = NaN; 
            xCoordsDwell = NaN; 
            yCoordsDwell = NaN; 
            dwell = NaN; 
            dwellMasks(:,:,1) = zeros(size(roiMask)); 
            dwellMaskCum = zeros(size(roiMask)); 
    
    end
    save([subRoiDir filesep 'dwellMasks.mat'],'dwellMasks','dwellMaskCum');
  
   
end % if refinedwell
%% Extract Data Pause
 
%Get Pause Data
idxPause=find(dataMatMergeAll(:,5)==2);

if ~isempty(idxPause)
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
%the lifetime(ie number of frames)
%inside the ROI is just the sum of each row (sum in the 2nd D)-1
insideSecPause =sourceProjData.secPerFrame.*(sum(INPause,2)-1);


% Get Track IDs for those Growth Sub-Tracks in the Sub-ROI of Interest
if ~isempty(strcmpi(timeUnitsPause,'fraction'))
    trckIdxInPause=intersect(find(insideSecPause./lifeSecPause>=timeValPause),inIncludeRegionPause);
    %trckIdxOutPause=intersect(find(insideSecPause./lifeSecPause<timeVal & insideSecPause > 0),inIncludeRegionPause);
else
    trckIdxInPause=intersect(find(insideSecPause>=timeValPause),inIncludeRegionPause);
    %trckIdxOutPause=intersect(find(insideSecPause<timeVal & insideSecPause > 0),inIncludeRegion);
end


% limit data to these tracks
xMatPause=xMatPause(trckIdxInPause,:);
yMatPause=yMatPause(trckIdxInPause,:);


% plot all member tracks in red on top of the mask
figure('Visible',visible);
imagesc(roiMask);
colormap 'gray'; 
axis off
hold on;
plot(xMatPause',yMatPause','b')
forTitleTime = ['LIFETIME CRITERIA: ' upper(timeUnitsPause) ' GREATER OR EQUAL TO ' num2str(timeValPause) x];   
 title({forTitleTime; projNameTitle; 'Pauses In SubRoi'});
saveas(gcf,[subRoiDir filesep 'extract_plots' filesep 'PausesInSubRoi.eps'],'psc2')
close(gcf)

end % ~isempty
%% Extract Data Shrinkage

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
if ~isempty(strcmpi(timeUnitsShrink,'fraction'))
    trckIdxInShrink=intersect(find(insideSecShrink./lifeSecShrink>=timeValShrink),inIncludeRegionShrink);
else
    trckIdxInShrink=intersect(find(insideSecShrink>=timeValShrink),inIncludeRegionShrink);
    
end

% limit data to these tracks
xMatShrink=xMatShrink(trckIdxInShrink,:);
yMatShrink=yMatShrink(trckIdxInShrink,:);


% plot all member tracks in black on top of the mask
figure('Visible',visible);
imagesc(roiMask);
colormap 'gray'
axis off
hold on;
plot(xMatShrink',yMatShrink','k')
 if strcmpi(timeUnitsShrink,'fraction'); 
            x = '';
        else 
            x = ' sec'; 
 end 
forTitleTime = ['LIFETIME CRITERIA: ' upper(timeUnitsShrink) ' GREATER OR EQUAL TO ' num2str(timeValShrink) x]; 
  %title({projNameTitle;'Shrinkages In SubRoi';forTitleTime})
  title({'Shrinkages In SubRoi' ; forTitleTime}); 
%saveas(gcf,[subRoiDir filesep 'ShrinkagesInSubRoi' fileExt])
saveas(gcf,[subRoiDir filesep 'extract_plots' filesep 'ShrinkagesInSubRoi.eps'],'psc2');
close(gcf)


%% Limit IN to only those subtracks in region based on partitioning
IN = IN(trckIdxIn,:); 
OUT = OUT(trckIdxIn,:); % matrix with ones where the subtrack associated with the 
% region lies OUTSIDE the region (useful for looking at speed/lifetime 
% as the microtubule crosses some user-defined border. 
%% Write Stats to Output


projData.anDir=subRoiDir; % path to sub-roi
projData.imDir=MD.getChannelPaths{1};

    if multMasks ~= 1 % quick fix! 
    xCoordIn = xMatIn.*IN; % calc coordinates before removing beg and end for plotting
    yCoordIn = yMatIn.*IN; % calc coordinates before removing beg and end for plotting
  
    else 
        
        
        xCoordIn = xMatIn.*IN_MultMasks; 
        yCoordIn = yMatIn.*IN_MultMasks; 

        xCoordOut = xMatIn.*OUT_MultMasks;
        yCoordOut = yMatIn.*OUT_MultMasks;
    end
    
    
projData.xCoordInAllTracks = xCoordIn;  % save these coordinates before removing beg and end for plotting
projData.yCoordInAllTracks = yCoordIn;

projData.xCoordOutAllTracks= xCoordOut ;
projData.yCoordOutAllTracks = yCoordOut;

% Save the Dwell Coordinates if applicable : 
if ip.Results.performDwellCalcs
    projData.xCoordDwellOutAllTracks = xCoordsDwellOut;
    projData.yCoordDwellOutAllTracks = yCoordsDwellOut;
    
    projData.xCoordDwell = xCoordsDwell;
    projData.yCoordDwell = yCoordsDwell;
end

if ~isempty(trckIdxIn)
    projData.nTracksSubRoi = length(unique(dataMatMerge(trckIdxIn,1))); % number of
    % compound tracks represented (NOTE: depending on the user specifications the
    %complete compound track may not be in the region)
    if isfield(projData.stats,'percentTracksStartAndEndInRegion')
    projData.stats.percentTracksStartAndEndInRegion = length(trckIdxInRegion)/length(trckIdxIn)*100;
    end 
else
    projData.nTracksSubRoi = 0;
    projData.stats.percentTracksStartAndEndInRegion = NaN;
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
if multMasks == 0 
% specifically write those compound tracks with an Fgap or Bgap included in
% region
if ~isempty(idxPause)
    if ~isempty(trckIdxInPause)
        projData.tracksWithFgapSubRoi = dataMatMergePause(trckIdxInPause,1);
        projData.tracksWithBgapSubRoi = dataMatMergeShrink(trckIdxInShrink,1);
    else
        projData.tracksWithFgapSubRoi = 0;
        projData.tracksWithBgapSubRoi = 0;
    end
end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Density Calcs: Calculate Median Nearest Neighbor For SubRoi

[numSubTracks, numFrames] = size(xMatIn);

% find the number of detected particles per frame
cometNum = nan(numFrames,1); 
if ~isempty(xMatIn)
    for iFrame = 1: numFrames
    cometNum(iFrame) = sum(~isnan(xMatIn(:,iFrame)));  
    end 
    
end 


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

if projData.nTracksSubRoi~=0
    % get frame-to-frame displacement for growth only (not forward/backward gaps)
    frame2frameDispPix=sqrt(diff(projData.xCoord,1,2).^2+diff(projData.yCoord,1,2).^2);
    % get rid of NaNs and linearize the vector
    projData.frame2frameDispPix=frame2frameDispPix(~isnan(frame2frameDispPix(:)));
else
    projData.frame2frameDispPix=NaN;
end

if projData.nTracksSubRoi~=0
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


%% Create a DataMat of just SubRoi Tracks for Stats

dataGrowthROI = dataMatMerge(trckIdxIn,:);
if ~isempty(idxPause)
dataPauseROI = dataMatMergePause(trckIdxInPause,:);
else 
    dataPauseROI = [];
end 
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

% Include only those tracks in region that have been specified by partitioning
lifeSec = lifeSec(trckIdxIn);  
insideSec = insideSec(trckIdxIn); 
projData.insideSecAllTracks = insideSec; % save these values before removing subtracks
% at the beginning or end of movie (for plotting)
if ip.Results.performDwellCalcs
    projData.dwellAllTracks = dwell'; % for now this seems to be allrit 
end 

% calculate speed in before removing the beginning and the end tracks (for
% plots)
speedIn=nanmean(sqrt(diff(xMatIn.*IN,[],2).^2+diff(yMatIn.*IN,[],2).^2),2);
speedOut=nanmean(sqrt(diff(xMatIn.*OUT,[],2).^2+diff(yMatIn.*OUT,[],2).^2),2);

projData.speedInMicPerMinAllTracks=pixPerFrame2umPerMin(speedIn,projData.secPerFrame,projData.pixSizeNm);
projData.speedOutMicPerMinAllTracks=pixPerFrame2umPerMin(speedOut,projData.secPerFrame,projData.pixSizeNm);

% these values will be used for plotting

%%
% Remove Growth SubTracks that begin in frame
% first frame or end in last frame if user specifies
% Adjust lifeSec and insideSec column to include only those sub-tracks included
% in stats

if remBegEnd== 1
    [dataTotROI] = plusTipRemBegEnd(dataTotROI, projData);
    projData.remBegEnd = 'yes';
    
    % note for the user the percentage of growths at the start or end that
    % they are removing
    projData.startOrEnd=~isnan(xMatIn(:,1)) | ~isnan(xMatIn(:,end));
    %projData.percentGrowthAtStartOrEnd=sum(projData.startOrEnd)./projData.stats.nGrowths*100;
    
    xMatInAll = xMatIn;
    yMatInAll = yMatIn; 
    % modify lifeSec and insideLifeSec to only include those included in stats
    
    lifeSec = lifeSec(~projData.startOrEnd); % remove subtracks at the beg or end of  movie for stats
    insideSec = insideSec(~projData.startOrEnd);
    
    
    
    % modify the coordinate matrix and the lifeSec matrix so can calculate
    % some other stats from cropped data structure
    
    xMatIn = xMatIn(~projData.startOrEnd,:);
    yMatIn = yMatIn(~projData.startOrEnd,:);
    IN = IN(~projData.startOrEnd,:);
    OUT = OUT(~projData.startOrEnd,:);
    
    % plot so user knows all tracks in the beginning or end of movie that
    % have been removed
    figure('Visible',visible)
    imagesc(roiMask)
    colormap('gray')
    axis off
    hold on 
    plot(xMatInAll',yMatInAll','y');
    plot(xMatIn',yMatIn','r');
    title({'Note: You have chosen to remove subtracks that start in first frame and end in last frame';...
        'Yellow: Removed Growth Subtracks (Subtracks That Start In First Frame or End in Last Frame )';'Red: Maintained Growth Subtracks- Final For Analysis'})
    saveas(gcf,[subRoiDir filesep 'extract_plots' filesep 'removeBegEndPlot.eps'],'psc2'); 
    close(gcf); 
    
    
    
    
    %
%     speedIn=nanmean(sqrt(diff(xMatIn.*IN,[],2).^2+diff(yMatIn.*IN,[],2).^2),2);
%     speedOut=nanmean(sqrt(diff(xMatIn.*OUT,[],2).^2+diff(yMatIn.*OUT,[],2).^2),2);
    
%     speedInMicPerMin=pixPerFrame2umPerMin(speedIn,projData.secPerFrame,projData.pixSizeNm);
%     speedOutMicPerMin=pixPerFrame2umPerMin(speedOut,projData.secPerFrame,projData.pixSizeNm);
    speedInMicPerMin = projData.speedInMicPerMinAllTracks(~projData.startOrEnd);
    speedOutMicPerMin = projData.speedOutMicPerMinAllTracks(~projData.startOrEnd); 
    
else
    if isfield(projData,'dataMatCropped_MicPerMin_Sec_Mic')== 1;
        projData =  rmfield(projData,'dataMatCropped_MicPerMin_Sec_Mic'); % a field in an older version
    end
    projData.remBegEnd = 'no';
    
    speedInMicPerMin = projData.speedInMicPerMinAllTracks;
    speedOutMicPerMin = projData.speedOutMicPerMinAllTracks;
end

%% FROM NOW ON ALL x/yMatIn,IN,OUT ARE FILTERED TO REMOVE BEGINNING END TRACKS IF THAT OPTION ON!!!%%%%%
%%
% WRITE NEW FIELDS NOT IN SOURCE (WHOLE CELL) PROJDATA.MAT
if projData.nTracksSubRoi~=0
    pixSizMic=projData.pixSizeNm/1000; % size of a pixel in microns
    pixAreaSqMic=pixSizMic^2; % area of a pixel in square microns
    cellAreaSqMic=sum(wholeCellRoiMask(:)).*pixAreaSqMic;
    % calculate area in square microns of this roi and its
    % percentage of the full cell's area
    subRoiAreaSqMic=sum(roiMask(:)).*pixAreaSqMic;
    percentRoiArea=100*(subRoiAreaSqMic/cellAreaSqMic);
    projData.subRoiAreaSqMic=subRoiAreaSqMic;
    projData.percentRoiArea=percentRoiArea; 
    projData.stats.CometDensityPerFrame_mean = mean(cometNum./subRoiAreaSqMic); % cometNum is a vector of with the number of comets in region per frame
    projData.stats.CometDensityPerFrame_std = std(cometNum./subRoiAreaSqMic);
   
    
    % record in the projData a list of the lifetime of each growth subtrack
    % ONLY within the subregion of interest. (the subroi partitioning
    % allows you to place a long track to a specific subregion but this value
    % will calculate two local parameters(speed and growth lifetime)
    % within this region)
    
%     projData.lifeSec=lifeSec; % total lifetime (seconds) % only those cacluated in stats
%     projData.insideSec=insideSec; % lifetime within sub-roi (seconds)
    percentLifeInside=100*(insideSec./lifeSec); % percent time within sub-roi
%     projData.speedInMicPerMin = speedInMicPerMin;
%     projData.speedOutMicPerMin = speedOutMicPerMin; 
%     
    numInsideWholeLife = length(find(percentLifeInside == 100));
    
    projData.stats.percentTracksStartAndEndInRegion = (numInsideWholeLife/...
        length(percentLifeInside))*100;
    
    % format a table for regional output 
    titles = cell(1,6); 
    titles{1} = 'ID of CompoundTrack to Which SubTrack Belongs';
    titles{2} = 'Percent Time of Growth SubTrack Inside Region';
    titles{3} = 'Lifetime of Growth SubTrack Inside Region (s)'; 
    titles{4} = 'Total Lifetime of SubTrack (s)';
    titles{5} = 'Averge Frame-to-Frame Velocity of Growth SubTrack In Region (um/min)'; 
    titles{6} = 'Average Frame-to-Frame Velocity of Growth SubTrack Out of Region (um/min)';
    titles{7} = 'Ratio: Speed IN / Speed OUT';
    
    ratio = speedInMicPerMin./speedOutMicPerMin; 
    idxGrowth = dataTotROI(dataTotROI(:,5) == 1,1); 
    inRegionValues= [idxGrowth percentLifeInside insideSec lifeSec speedInMicPerMin speedOutMicPerMin ratio];
    
    
    
    inRegionValues = num2cell(inRegionValues); 
    projData.dataMatSubRoi_CompareGrowthInToOut = [titles ; inRegionValues]; 
    
    
    % record the median and mean values of the lifetime of the
    % the subdivided growth events
    
    projData.stats.growth_lifetime_median_INSIDE_REGION = median(insideSec);
    projData.stats.growth_lifetime_mean_INSIDE_REGION = mean(insideSec);
    
    projData.stats.growth_speed_median_INSIDE_REGION = nanmedian(speedInMicPerMin);
    projData.stats.growth_speed_mean_INSIDE_REGION = nanmean(speedInMicPerMin);
    
    %% if user desires filter based on displacement
   
    dispInside = speedInMicPerMin.*(insideSec/60); 
    
    insideSecLowDisp = insideSec(dispInside<2.5);  % take with only less then x displacement
  
    projData.stats.growth_lifetime_mean_INSIDE_REGION_lowDisp = mean(insideSecLowDisp); 
    projData.stats.growth_lifetime_median_INSIDE_REGION_lowDisp = median(insideSecLowDisp); 
    
else % output if no tracks
    
%     projData.lifeSec=NaN;
%     projData.insideSec=NaN;
%     projData.percentLifeInside=NaN;
%     
%     projData.speedInMicPerMin=NaN;
%     projData.speedOutMicPerMin=NaN;
%     
    
    projData.startOrEnd=NaN;
    projData.percentGrowthAtStartOrEnd=NaN;
    projData.percentTracksStartAndEndInRegion = NaN;
    projData.stats.growth_speed_mean_INSIDE_REGION = NaN;
    projData.stats.growth_speed_median_INSIDE_REGION = NaN;
    projData.stats.growth_lifetime_median_INSIDE_REGION = NaN;
    projData.stats.growth_lifetime_mean_INSIDE_REGION = NaN;
    projData.stats.CometDensityPerFrame_mean = NaN;
    projData.stats.CometDensityPerFrame_std = NaN; 
end



%% Calculate the stats
% Note the index of one at the end tells the program it is calling this function
% from subRoiAnalysis and to not calculate certain params that require
% intact compound tracks

[projData,M]=plusTipDynamParam(dataTotROI,projData,0,1);
%[projData,M,idxPer,idxPar]=plusTipDynamParam(dataTotROI,projData,0,1);

% update saved M
projData.M = M;

%% Two extra subRoi parameters : (Only for Us)
if projData.nTracksSubRoi ~=0; 
    projData.stats.percentWholeCellTracksNucInSubRegion = projData.stats.numNucleationEvents/nTracks; 
    projData.stats.percentWholeCellTracksTermInSubRegion = projData.stats.nGrowthTermEvents/nTracks; 
else 
    projData.stats.percentWholeCellTracksNucInSubRegion = NaN; 
    projData.stats.percentWholeCellTracksTermInSubRegion = NaN; 
end 

%% Make Distribtuion Plots
histDir = [subRoiDir filesep 'meta' filesep 'histograms'];
if turnHistOn == 1
    if ~isempty(xMatIn)
        if projData.nTracksSubRoi ~= 0
            plusTipMakeHistograms(M,[subRoiDir filesep 'meta' filesep 'histograms']);
        end
    end
end




%% Make Polarity Plots
if mkPolPlot == 1
    if isempty(xMatIn); % skip
        projData.subTrackAnglesIndMean = NaN;
        projData.subTrackAnglesIndStd = NaN;

        projData.stats.polarCoordMeanOfAllSubtracks = NaN; 
        projData.stats.polarCoordMedianOfAllSubtracks = NaN; 
        projData.stats.polarCoordStdOfAllSubtracks = NaN;
        
        
    else
        if ~isdir(histDir);
            mkdir(histDir);
        end
       
        [anglesFinal,subTrackAnglesIndMean,subTrackAnglesIndStd] = plusTipPlotTrackAngles2(xMatIn,yMatIn);
        
        projNameTit = regexprep(projName,'_',' ');
        
        saveFig1 =figure('Visible',visible);
        %normalize
        [tout rout]  =  rose(anglesFinal); % based on all frame to frame angles
        rho = rout/sum(rout);
        polar(tout,rho);
        projData.anglesFrame2Frame = anglesFinal;
        
        title({projNameTit; strrep(sub,'_','');  'Rose Plot of All Frame-to-Frame Growth Displacement (Converted to Theta)'});
        saveas(saveFig1,[histDir filesep 'angles_histogram.eps'], 'psc2');
       
        if collectPlots == 1 
            
        collectPathOrientTif = [collectedDataPath filesep projName filesep 'tifs'];
        collectPathOrientEPS = [collectedDataPath filesep projName filesep 'eps']; 
        
        if ~isdir(collectPathOrientTif)
            mkdir(collectPathOrientTif);
            mkdir(collectPathOrientEPS);
        end 
        
        saveas(saveFig1,[collectPathOrientEPS filesep sub 'Rose Plot Frame to Frame.eps'],'psc2');
        saveas(saveFig1,[collectPathOrientTif filesep sub 'Rose Plot Frame to Frame.tif']); 
        
        end
        
        
        close(saveFig1);
        
      
        projData.subTrackAnglesIndMean = [dataTotROI(dataTotROI(:,5)==1,1), subTrackAnglesIndMean]; 
        projData.subTrackAnglesIndStd = [dataTotROI(dataTotROI(:,5)==1,1),subTrackAnglesIndStd]; 


        projData.stats.polarCoordMeanOfAllSubtracks = nanmean(projData.subTrackAnglesIndMean(:,2)); 
        projData.stats.polarCoordMedianOfAllSubtracks = nanmedian(projData.subTrackAnglesIndMean(:,2)); 
        projData.stats.polarCoordStdOfAllSubtracks = nanstd(projData.subTrackAnglesIndMean(:,2)); 

        xCoordIn = xMatIn.*IN;
        yCoordIn = yMatIn.*IN; 
        [anglesFinalInside,subTrackAnglesIndMeanInside] = plusTipPlotTrackAngles2(xCoordIn,yCoordIn); 



        projData.subTrackAnglesIndMean_INSIDE_REGION = [dataTotROI(dataTotROI(:,5)==1,1),subTrackAnglesIndMeanInside]; 

        
        
        
        
    end
end 



%% angle data



 





% Now write the interpolated coordinates corresponding to pause (necessary
% for plotting functions)
if ~isempty(idxPause)
    for iSub=1:length(trckIdxInPause)
        k=trckIdxInPause(iSub);
        
        projData.xCoord(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3))=sourceProjData.xCoord(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3));
        projData.yCoord(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3))=sourceProjData.yCoord(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3));
        
        %     projData.featArea(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3))=sourceProjData.featArea(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3));
        %     projData.featInt(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3))=sourceProjData.featInt(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3));
        
        projData.frame2frameVel_micPerMin(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3)-1)=sourceProjData.frame2frameVel_micPerMin(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3)-1);
        projData.segGapAvgVel_micPerMin(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3)-1)=sourceProjData.segGapAvgVel_micPerMin(dataMatMergePause(k,1),dataMatMergePause(k,2):dataMatMergePause(k,3)-1);
    end
end
%%
%Now write the interpolated coordinates corresponding to shrink (necessary
%for plotting functions)

for iSub=1:length(trckIdxInShrink)
    k=trckIdxInShrink(iSub);
    
    projData.xCoord(dataMatMergeShrink(k,1),dataMatMergeShrink(k,2):dataMatMergeShrink(k,3))=sourceProjData.xCoord(dataMatMergeShrink(k,1),dataMatMergeShrink(k,2):dataMatMergeShrink(k,3));
    projData.yCoord(dataMatMergeShrink(k,1),dataMatMergeShrink(k,2):dataMatMergeShrink(k,3))=sourceProjData.yCoord(dataMatMergeShrink(k,1),dataMatMergeShrink(k,2):dataMatMergeShrink(k,3));
    
   
    projData.frame2frameVel_micPerMin(dataMatMergeShrink(k,1),dataMatMergeShrink(k,2):dataMatMergeShrink(k,3)-1)=sourceProjData.frame2frameVel_micPerMin(dataMatMergeShrink(k,1),dataMatMergeShrink(k,2):dataMatMergeShrink(k,3)-1);
    projData.segGapAvgVel_micPerMin(dataMatMergeShrink(k,1),dataMatMergeShrink(k,2):dataMatMergeShrink(k,3)-1)=sourceProjData.segGapAvgVel_micPerMin(dataMatMergeShrink(k,1),dataMatMergeShrink(k,2):dataMatMergeShrink(k,3)-1);
end

%% Record SubRoi Extraction Parameters 
projData.subRoiParams.extractTypeGrowth = extractType; 
projData.subRoiParams.timeUnitsGrowth = timeUnitsGrowth; 
projData.subRoiParams.timeValGrowth = timeValGrowth; 
projData.subRoiParams.timeUnitsPause = timeUnitsPause; 
projData.subRoiParams.timeValPause = timeValPause; 
projData.subRoiParams.timeUnitsShrink = timeUnitsShrink; 
projData.subRoiParams.timeValShrink = timeValShrink; 

%% SAVE
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

end