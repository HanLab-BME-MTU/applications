function [dataMatMerge,dataMatReclass,dataMatCrpSecMic,percentFgapsReclass,percentBgapsSlow,percentBgapsShortTime, percentBgapsVeryFast,avgLatSec, tracksBGapShortTime,tracksBGapSlow,tracksBGapVeryFast,cutoffValueFGap,fgapReclassScheme,bgapReclassScheme]=plusTipMergeSubtracks(projData,dataMat)
% plusTipMergeSubtracks merges growth fgaps with the flanking growth phases
%
% NOTE: THIS IS A WORKING COPY: where I am attempting to fix problems with 
% the reclassification scheme: therefore there might and likely is bugs
% talk to me if you intend to use MB 03/27/11
% 
% SYNOPSIS  : [dataMatMerge,dataMatReclass,percentFgapsReclass]=...
%                   plusTipMergeSubtracks(projData,dataMat)
%
% INPUT
% projData  : structure containing frame rate, pix size info, etc. 
%             if dataMatrix isn't given, then projData should have the
%             field nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix
%             which contains reclassified fgaps and bgaps.
% dataMat   : mainly to be used in plusTipPostTracking, where it is the
%             data matrix prior to reclassification.
%             
% OUTPUT
% dataMatMerge        : matrix where fgaps that should be reclassified (as
%                       growth, because their speeds are
%                       >=70% the speed at the end of the growth phase
%                       prior) are consolidated with flanking growth
%                       subtracks. also bgaps that should be reclassified
%                       as pause, because their speeds are lower than the
%                       95th percentile of fgap speeds (after reclass),
%                       have their indices changed to 2.
% dataMatReclass      : matrix where fgaps that should be reclassified are
%                       not consolidated but the track type is changed to
%                       5, and bgaps that should be reclassified have their
%                       track type changed to 6.
% dataMatCrpSecMic    : like dataMatMerge except the growth phases and
%                       linked fgaps/bgaps from the first and last frames
%                       of the movie have been removed. also, column 6
%                       represents lifetime in seconds and 7 represents
%                       displacement in microns. all speeds and displacements
%                       are positive (which makes more sense anyway for
%                       these measurements)
% percentFgapsReclass : percentage of fgaps that get reclassified as
%                       continuation of growth
% percentBgapsReclass : percentage of bgaps that get reclassified as
%                       pause
%
% 
% this function is called by:
% plusTipPostTracking
% plusTipPoolGroupData
% plusTipParamPlot
% plusTipGetSubtrackCoords
% plusTipSpeedMovie
%% Specify Reclassification Schemes (eventually put into input of function)

% CHOICES FOR FGAP RECLASSIFICATION

localFGapReclassNew = 1; % quick fix for now make option of reclass scheme later
useFullGrowthSubTrack = 1; % easiest without artifacts to just use full growth 
% subtrack for now, in the end would like to eliminate the last 2-3 frames 
% before the pause event (the latency time associated with dissociation of
% the comet)



% Choices for bgap reclassification

unimodalReclass = 0;  % if 1 perform unimodalReclassification of fgaps (right now
% implementation is an iterative scheme -- works but not sure logicaly justified)

% note choices of fgap reclassification scheme for output
if (localFGapReclassNew == 1 && useFullGrowthSubTrack == 1)
    fgapReclassScheme = 'Local Scheme: Full Growth Subtrack';
elseif unimodalReclass == 1
    fgapReclassScheme = 'New Iterative UniModal';
end 

% adjust some of the output based on which reclassification was chosen                  
if unimodalReclass ~= 1
    cutoffValueFGap = NaN; % there is no global cut-off value for fgap (we use local cut-offs instead)
else 
end 
    

% CHOICES FOR BGAP RECLASSIFICATION
nonsymmetric = 0; % perform new reclassification scheme for bgaps taking into 
                  % consideration the assymetrical effects of the latency                 
                  % of comet formation. particularly necessary at fast 
                  % sampling time scales, otherwise too many bgaps 
                  % will be reclassified as pause events: default uses the 
                  % the 95 percentile of "true fgaps" speeds as a cut off value (all bgaps with a 
                  % speed slower than this value are considered pauses. 

%note choices of bgap reclassification scheme for output                 
if nonsymmetric == 1
bgapReclassScheme = 'New Non-Symmetric';
else 
    bgapReclassScheme = '95th percentile fgap speeds';
end; 

              
if nonsymmetric == 0
    tracksBGapShortTime = NaN; % this value is only applicable if you reclassify based on the nonsymmetric method
    tracksBGapSlow = NaN;
    tracksBGapVeryFast= NaN;
    percentBgapsShortTime = NaN;
    percentBgapsSlow = NaN;
    percentBgapsVeryFast = NaN;
else 
end 
   

%% Check Inputs and Perform Local Fgap Reclassification If Desired

if nargin<2
    % dataMat is the output matrix in projData, where the gaps to
    % conolidate are labeled as trackType=5
    dataMat=projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix;
    dataMatInput = dataMat; % save an unadulteraterated version of the data
    % the set of all fgaps in this case have trackType of 2 or 5
    fgapIdx=find(dataMat(:,5)==2 | dataMat(:,5)==5);
    
    % these are the fgaps to consolidate
    growthFgapIdx=find(dataMat(:,5)==5);
    
    % Calculate percentage of fgaps reclassified 
    if isempty(fgapIdx)
    percentFgapsReclass=NaN;
    else
    percentFgapsReclass=100*length(growthFgapIdx)/length(fgapIdx);
    end
else
   
    dataMatInput = dataMat; % save an unadulterated version of the data
    % dataMat is an input, we need to find where the growthFgaps are

    % look at fgaps - if their speed << growth phase just prior, then it's a
    % true pause.  it it is greater than 70% growth speed, reclassify it as
    % continuation of growth.
    fgapIdx=find(dataMat(:,5)==2);
    
    beforeFgapIdx=fgapIdx-1;
    % get speed of growth phase at last 2-3 time pts (depending on how long
    % growth phase lasted) just prior to pause
    eF=dataMat(beforeFgapIdx,3);% get the end frame of subtrack right before pause
    beforeFgapSpeed=zeros(length(fgapIdx),1); 
    
   % perform new local reclassification scheme for the entire growth
   % subtrack
   
   if localFGapReclassNew == 1 && useFullGrowthSubTrack == 1
    % need first start frame of subtrack
        for iGap = 1:length(fgapIdx)
           beforeFgapSpeed(iGap) = dataMat(beforeFgapIdx(iGap),4); % just use the average 
           % frame-to-frame velocity of the entire growth subtrack for now
        end % end for iGap
           
    
   else  % default is to proceed with old reclassification scheme
    % uses avg frame to frame velocity from 2-3 frames of the growth 
    % subtrack just before the pause.  Problem with this scheme is that 
    % if the velocity of the comet significantly slows right before pausing
    % (will especially happen if there is a delay in the dissociation of
    % the comet upon a pause) potentially viable pause information may
    % become reclassified as an undetected growth event.  
    
    %NOTE it appears the old scheme seems to have a bit of a bug in that if a 
    % bgap proceeds a short growth subtrack it will result in those frame
    % to frame velocities being included in the beforeFgapSpeed estimation
    % making this estimate negative and thus the fgap will automatically be
    % reclassified. (MB: 03/11)
    
    for iGap=1:length(fgapIdx)
        idxRange=[max(1,eF(iGap)-3):eF(iGap)-1]; % not sure why she uses max here
        beforeFgapSpeed(iGap)=nanmean(projData.frame2frameVel_micPerMin(dataMat(beforeFgapIdx(iGap),1),idxRange)); 
      
    end
    
   end % end iflocalFGapReclassNew
    % these are the fgaps to consolidate
    growthFgapIdx=fgapIdx(dataMat(fgapIdx,4)>0.7.*beforeFgapSpeed);
    
    
    % Calculate percentage of fgaps reclassified 
    if isempty(fgapIdx)
    percentFgapsReclass=NaN;
    else
    percentFgapsReclass=100*length(growthFgapIdx)/length(fgapIdx);
    end
   
end % end if nargin


% Save Local Reclassification fgaps (performed above):

% Have already performed the localreclassification above, if unimodalReclass 
% not selected-  then use this scheme

if unimodalReclass ~= 1
dataMatReclass=dataMat;
dataMatReclass(growthFgapIdx,5)=5; % these are the fgaps to be reclassified as growths             
else % we'll perform the unimodal reclassification of fgpas (a global reclassification scheme) below
end 


%% Unimodal Thesholding Reclassification Scheme for FGAPS 
% Below is the new iterative unimodal reclassification scheme which works 
% but not convinced is justified.

if unimodalReclass == 1  
    
    
    %Initiate a data variable with  unique name for manipulation 
    dataMatReclassUniMode = dataMatInput;
  
    
    % the set of all fgaps in this case have trackType of 2 or 5
    fgapIdx=find(dataMatInput(:,5)==2 | dataMatInput(:,5)==5);
    
     
    dataMatReclassUniMode(:,6)=dataMatReclassUniMode(:,6).* projData.secPerFrame; % convert lifetimes to seconds
    dataMatReclassUniMode(:,7)=dataMatReclassUniMode(:,7).*(projData.pixSizeNm/1000); % convert displacements to microns

    delta = 1; %  put in a delta value to initiate 
    cutoffValueFGap = 0; % first cutoffValueFGap = 0  
    
    convgThresh = 0.0001;
    
    while delta > convgThresh
    
    % Reassign all fgaps so that they are again all classified as "true
    % pauses" = ID of 2 (just in case some form of reclassification has 
    % already been performed undo this)
    dataMatReclassUniMode(fgapIdx,5) = 2;   
        
        
    cutoffValueFGapBC = cutoffValueFGap;   % for the first iteration this will be zero
    
    % Collect all fgap speeds
    fgapSpeeds = dataMatReclassUniMode(fgapIdx,4);
    %fgapSpeeds = fgapSpeeds(fgapSpeeds > 0); % do not include negative values in histogram 
        
    % perform unimodal thresholding using all fgap speeds to obtain 
    % maximum fgap speed corresponding to a pause event (above this thresh
    % are likely undetected growth events)
    [cutoffIdx, cutoffValueFGap, sp, axesH] = cutFirstHistMode(fgapSpeeds,1);
     
    % calculate the difference between the value obtained before and after 
    % correction of fgap velocity for comet latency
    delta = abs(cutoffValueFGapBC - cutoffValueFGap); % when this value = the convergence thresh the loop will stop 
    % use previous values for growthFgapIdx etc as should be basically the same as 
    % last iteration

    % Mark the fgaps for reclassification to growth under unimodal scheme
    % Criteria for fgap reclass: If velocity  of the fgap is greater than cuttoff value from
    % unimodal threshholding of the fgap velocity histogram
    % the fgap is reclassified as an undetected growth event (subtrack ID = 5) 
    growthFgapIdx = find(dataMatReclassUniMode(:,4) > cutoffValueFGap & dataMatReclassUniMode(:,5) == 2);
    
    % Mark all pauses that should be reclassified as growth 
    % column 5 subTrack ID from 2 --> 5) 
    dataMatReclassUniMode(growthFgapIdx,5) = 5; 
    
   
    
    avgVelGrowth = mean(dataMatReclassUniMode(dataMatReclassUniMode(:,5) == 1 | dataMatReclassUniMode(:,5) == 5, 4)); 
    
    
    % use "true" fgaps (subtrack ID = 2) obtained from unimodal thresh   
    % to calculate the avg latency of the comet formation 
    avgDispPause = mean(dataMatReclassUniMode(dataMatReclassUniMode(:,5) == 2,7)); % in microns

    % calc avg latency of comet formation  (in min)
    avgLat = avgDispPause/avgVelGrowth; % in minutes
    avgLatSec = avgLat*60; % in seconds
    

    
    
    % Correct the fgap velocity based on estimated avg latency of comet
    % formationv
    fgapTimeCorr = (dataMatReclassUniMode(fgapIdx,6) - avgLatSec)/60; % correct lifetime for all fgaps (pre-any reclassification) ans in min
   
    fgapDispCorr = (dataMatReclassUniMode(fgapIdx,7) - avgDispPause);
    
    for i = 1:length(fgapIdx)
    dataMatReclassUniMode(fgapIdx(i),4) = fgapDispCorr(i)/fgapTimeCorr(i); % recalculate fgap velocity um/min
    end % end for i
    end % end while  
    
    %[cutoffIdx, cutoffValueFGap, sp, axesH] =
    %cutFirstHistMode(fgapSpeeds); % his of the final fgap speeds corrected
    
   truePauseIdx = find((dataMatReclassUniMode(:,5) == 2 | dataMatReclassUniMode(:,5) == 5)...
       & dataMatReclassUniMode(:,4) < cutoffValueFGap);
   
   growthFgapIdx = find((dataMatReclassUniMode(:,5) == 2 | dataMatReclassUniMode(:,5) == 5) ... 
       & dataMatReclassUniMode(:,4) > cutoffValueFGap);
    
    for i = 1:length(fgapIdx)
        dataMatReclassUniMode(fgapIdx(i),4) = dataMatReclassUniMode(fgapIdx(i),7)/(dataMatReclassUniMode(fgapIdx(i),6)/60);
    end % end for i 
    
    figure; hist(dataMatReclassUniMode(truePauseIdx,4));
     
    % Rename for consistency with output
    dataMatReclass = dataMatReclassUniMode;
    
    % Officially reclassify all pauses that are marked for growth
    % (column 5 subTrack ID  from 2 --> 1)
    dataMatReclassUniMode(growthFgapIdx,5) = 1;
    % rename for input into merge reclassification scheme below 
    dataMat = dataMatReclassUniMode;
    
    
    % Calculate percentage of fgaps reclassified 
    if isempty(fgapIdx)
    percentFgapsReclass=NaN;
    else
    percentFgapsReclass=100*length(growthFgapIdx)/length(fgapIdx);
    end
    
     
else 
end % if unimodal reclass == 1
    
    
    
    



%% Merge Reclassifications (ie subTrack ID 5 --> 1 need to incorporate reclassified pauses in growth stats)
% merging will be performed regardless of reclassification scheme


% these are the affected track numbers for growth fgaps
tracks2check=unique(dataMat(growthFgapIdx,1));
rows2remove=[];

for i=1:length(tracks2check)
    % rows of dataMat corresponding to i track
    subIdx=find(dataMat(:,1)==tracks2check(i));
    % rows of dataMat corresponding to fgaps in track to consolidate
    fgap2remIdx=intersect(growthFgapIdx,subIdx);
    % rows of dataMat corresonding to bgaps or real pauses in track
    sepIdx=union(subIdx(dataMat(subIdx,5)==3),setdiff(intersect(subIdx,fgapIdx),fgap2remIdx))';
    % split the track based on bgaps, so that all fgaps that
    % should be consolidated together can be done dataMat the same time
    sIdx=[subIdx(1); sepIdx(:)];
    eIdx=[sepIdx(:); subIdx(end)];
    % loop through groups of subtracks (split by bgaps)
    for j=1:length(sIdx)
        % pTemp contains the ones to consolidate in this section
        pTemp=intersect(fgap2remIdx,sIdx(j):eIdx(j));
        if ~isempty(pTemp)
            fIdx=min(pTemp)-1; % first row - prior to first fgap
            lIdx=max(pTemp)+1; % last row - after final fgap

            dataMat(fIdx,3)=dataMat(lIdx,3); % change end frame
            dataMat(fIdx,6)=sum(dataMat(fIdx:lIdx,6)); % sum lifetimes
            dataMat(fIdx,7)=sum(dataMat(fIdx:lIdx,7)); % sum total displacements
            dataMat(fIdx,4)= pixPerFrame2umPerMin(dataMat(fIdx,7)/dataMat(fIdx,6),projData.secPerFrame,projData.pixSizeNm); % find new average velocity

            % keep track of which are the extra rows
            rows2remove=[rows2remove fIdx+1:lIdx];
        end
    end
end


% remove the extra rows
dataMat(rows2remove,:)=[];

% Data with pauses reclassified as growth merged in the growth velocity 
% stats (for output)
dataMatMerge=dataMat;


% perform lifetime and displacement unit conversions
dataMat(:,6)=dataMat(:,6).* projData.secPerFrame; % convert lifetimes to seconds
dataMat(:,7)=dataMat(:,7).*(projData.pixSizeNm/1000); % convert displacements to microns

% NOTE: dataMat now has all fgaps that are reclassified merged
% with their preceding growth subtrack.  Therefore, the indexing of 
% dataMat and dataMatReclass (before merging) will be different!
%% Old Reclassification Scheme for BGaps--- fix this later not using now!!!! 

if nonsymmetric ~= 1  % if unimodal is not checked proceed with old reclassification scheme for  bgaps
% get 95th-percentile of true pause speeds
fgapMaxSpeed=prctile(dataMatReclass(dataMatReclass(:,5)==2,4),95);

% get indices of bgaps where the speeds is less than the 95% cutoff
bgapAllIdx=find(dataMatReclass(:,5)==3 | dataMatReclass(:,5)==6);
bgap2pauseIdx=find((dataMatReclass(:,5)==3 & abs(dataMatReclass(:,4))<fgapMaxSpeed) | dataMatReclass(:,5)==6);

% bspeeds=abs(dataMatReclass(bgapAllIdx,4));
% [cutoffIndex, cutoffValue] = cutFirstHistMode(bspeeds,1);

dataMat(bgap2pauseIdx,5)=2; % reassign dataMat to have type 2 
dataMat(bgap2pauseIdx,:)=abs(dataMat(bgap2pauseIdx,:)); % get abs value since now pause
dataMatReclass(bgap2pauseIdx,5)=6; % reassign reclassified matrix to have new type of 6

% Calculate Fraction Reclassified 
% fraction of bgaps that are changed to pause, because their speeds are
% slower than the cutoff for fgap pauses
if isempty(bgapAllIdx)
    percentBgapsReclass=NaN;
else
    percentBgapsReclass=100*length(bgap2pauseIdx)/length(bgapAllIdx);
end

% fraction of fgaps that were consolidated into growth, since their speeds
% were more than 50% of the growth speed just prior to the gap
if isempty(fgapIdx)
    percentFgapsReclass=NaN;
else
    percentFgapsReclass=100*length(growthFgapIdx)/length(fgapIdx);
end

else % 
end % if nonsymmetric ~= 1

%% NEW bgap Reclassification Scheme (added by MB 03/16/11)

    
% calculate necessary avg values from above merged data (where all 
% reclassified pauses have been merged with preceding growth 
% subtrack

avgVelGrowth = mean(dataMat(dataMat(:,5) == 1,4));
avgDispPause = mean(dataMat(dataMat(:,5) == 2,7)); % in microns


% calc avg latency of comet formation  (in min)
avgLat = avgDispPause/avgVelGrowth; % in minutes
avgLatSec = avgLat*60; % in seconds

 
if nonsymmetric == 1 % proceed with the bgap reclassification based on new method which accounts for assymetric effect of 
                      % latency of comet formation on fgaps and bgaps

% initiate a new data variable from above merged data

dataMatCorrBgapStats = dataMat; % merged data
dataMatReclassCorrBgapStats = dataMatReclass; % non-merged data (retains
% info regarding fgap reclassification)

% correct bgap lifetime based on latency of comet formation (min)
% contributions to the time of bgap come from the time in shrinkage + 
% the time corresonding to the latency of comet formation
% attempt to correct for this by estimating a time for the latency of comet
% formation. 
bGapsIdx = find(dataMat(:,5) == 3);
useCos = 0; 
if useCos == 1
    
    xCoord = projData.xCoord; 
    yCoord = projData.yCoord;
   
    
    beforeBgapIdx = bGapsIdx -1;
    
    afterBgapIdx = bGapsIdx +1;
    
    beforeFrame1 = dataMat(beforeBgapIdx,3);
    beforeFrame2 = beforeFrame1-2;
    
    afterFrame1 = dataMat(afterBgapIdx,2); 
    
    for i = 1:length(bGapsIdx);
        trackID = dataMat(beforeBgapIdx(i),1);
        beforeX1 = xCoord(trackID,beforeFrame1(i));
        beforeY1 = yCoord(trackID,beforeFrame1(i));
        
        beforeX2 = xCoord(trackID,beforeFrame2(i)); 
        beforeY2 = yCoord(trackID,beforeFrame2(i));
        
        afterX1 = xCoord(trackID,afterFrame1(i));
        afterY1 = yCoord(trackID,afterFrame1(i));
        
        % first disp vector
        v1x = beforeX2 - beforeX1;
        v1y = beforeY2- beforeY1;
        v1mag = sqrt(v1x.^2+v1y.^2);
        
        % second disp vector 
        v2x = afterX1 - beforeX1;
        v2y = afterY1 - beforeY1;
        v2mag = sqrt(v2x.^2+v2y.^2);

        % cos of angle between consecutive vectors (displacements)
        cosV12=(v1x.*v2x+v1y.*v2y)./(v1mag.*v2mag);
        
        

        dataMatCorrBgapStats(bGapsIdx(i),7) = abs(dataMatCorrBgapStats(bGapsIdx(i),7)*cosV12)   + avgDispPause;
    end % end for i
        else 
            
        dataMatCorrBgapStats(dataMatCorrBgapStats(:,5) == 3, 7) = abs(dataMat(dataMat(:,5) == 3,7)) + avgDispPause;
            
end 
    
    

%if avgLat > dataMat(dataMat(:,5) == 3, 6)/60
dataMatCorrBgapStats(dataMatCorrBgapStats(:,5) == 3, 6) = dataMat(dataMat(:,5) == 3,6)/60 - avgLat; % convert seconds to minutes
%else 
    %dataMatCorrBgapStats(dataMatCorrBgapStats(:,5) == 3, 6) = dataMat(dataMat(:,5) == 3 ,6)/60;
%end 




% correct bgap displacment (um) based on the distance the microtubule
% likely traveled during comet formation. assuming two things here
% that comet formation after rescue is roughly on the order of comet
% formation after a pause.  correction based on average comet latency as 
% calculated from the "true" fgap distribution of displacements observed
 

%if avgLat > dataMat(dataMat(:,5) == 3, 6)/60
% contributions of bgap displacement are displacement during shrinkage - 
% the forward displacement of the microtubule that is undetected due to the
% latency of comet formation, if the latency of comet formation NOT
% neglible than attempt to correct the displacment value based on an
% estimate of the displacment the microutubule grew undetected 
% during comet formation 

%dataMatCorrBgapStats(dataMatCorrBgapStats(:,5) == 3, 7) = abs(dataMat(dataMat(:,5) == 3,7)) + avgDispPause;
%else 
%    dataMatCorrBgapStats(dataMatCorrBgapStats(:,5) == 3,6) = abs(dataMat(dataMat(:,5) == 3, 7));
%end 



% calculate the corrected bgap velocity in um/sec


% nonconsol data
bGapsIdxReclass = find(dataMatReclass(:,5) == 3);

for i = 1: length(bGapsIdx)
    dataMatCorrBgapStats(bGapsIdx(i),4) = dataMatCorrBgapStats(bGapsIdx(i),7)/(dataMatCorrBgapStats(bGapsIdx(i),6));
end 


% good to make histogram to trouble-shoot. 
bgapSpeedsCorr = dataMatCorrBgapStats(bGapsIdx,4);

% nonconsol data 
for i = 1:length(bGapsIdxReclass)
    dataMatReclassCorrBgapStats(bGapsIdxReclass(i),4) = bgapSpeedsCorr(i);
end 


% perform this to quickly get a hist of Corrected bgapSpeeds (will remove
% in the future)
%[cutoffIdx1, cutoffValue1, sp1, axesH1] = cutFirstHistMode(bgapSpeedsCorr);

% Assume shrinkage speeds are typically faster than the 0.7*average growth
% velocity (can modify this thresh to make more permissive/less
% permissive)

bGap2Pause = find(dataMatCorrBgapStats(:,5) == 3 & dataMatCorrBgapStats(:,4) < 0.7*avgVelGrowth & dataMatCorrBgapStats(:,4) >0);

%find the idx of those bgaps that when velocity is corrected 
% would yield reasonable shrinkage velocities (default between 0.7 of 
% avgGrowthVelocity and 4* the avg growth velocity)

bGapRealShrink = find(dataMatCorrBgapStats(:,5) == 3 & dataMatCorrBgapStats(:,4) > 0.7*avgVelGrowth & dataMatCorrBgapStats(:,4) < 4*avgVelGrowth);

% nonconsol data
bGap2PauseReclass = find(dataMatReclassCorrBgapStats(:,5) == 3 & dataMatReclassCorrBgapStats(:,4) < 0.7*avgVelGrowth & dataMatReclassCorrBgapStats(:,4) >0);
bGapRealShrinkReclass = find(dataMatReclassCorrBgapStats(:,5) == 3 & dataMatReclassCorrBgapStats(:,4) > 0.7*avgVelGrowth & dataMatReclassCorrBgapStats(:,4)<4*avgVelGrowth);


% Some shrinkage events are shorter than the avg latency of comet formation
% therefore their corrected bgap velocity is negative...for now leave these 
% as shrinkage events but document: we will likely make these shrinkage
% events as a separate category. 

bGapShortTime = find(dataMatCorrBgapStats(:,5) == 3 & dataMatCorrBgapStats(:,4) < 0);


bGapShortTimeReclass = find(dataMatReclassCorrBgapStats(:,5) == 3 & dataMatReclassCorrBgapStats(:,4) < 0);

% some shrinkage events when corrected yield velocities that are much
% higher than that anticipated based on the avg growth velocity 
% of the microtubule: segregate these out from other shrinkage events

bGapVeryFast = find(dataMatCorrBgapStats(:,5) == 3 & dataMatCorrBgapStats(:,4) > 4*avgVelGrowth);

bGapVeryFastReclass = find(dataMatReclassCorrBgapStats(:,5) == 3 & dataMatReclassCorrBgapStats(:,4) > 4*avgVelGrowth);


% Get the compound track ID of those tracks with the various types of bgaps 
% so can make movies of these later
% Would like to write this info to projData but need to find in what
% function she writes this 
%projData.tracksWithBgaps = dataMatCorrBgapStats(bGapRealShrink,1);
%projData.tracks2Pause = dataMatCorrBgapStats(bGap2Pause,1);
tracksBGapShortTime = dataMatCorrBgapStats(bGapShortTime,1);
tracksBGapSlow = dataMatCorrBgapStats(bGap2Pause,1);
tracksBGapVeryFast = dataMatCorrBgapStats(bGapVeryFast,1); 


% Mark bgaps for reclassification (col 5 ID 3(bgap) --> 6(bgap to be
% reclassified as a 'pause' (ie a very slow shrinkage event)

dataMatReclass(bGap2PauseReclass,5) = 6; % note these values now go into 

% Mark bgaps that are very fast shrinkage events as 7 
dataMatReclass(bGapVeryFastReclass,5) = 7 ; 


dataMatReclass(bGapShortTimeReclass,5) = 8; % keep these values separate for now 
% these are values where the bgap time is on the order of the average
% time of comet latency, making time in shrinkage on the order of zero,
% these event could represent a detection error OR a significant change 
% in the latency of comet formation to significantly faster values as compared to 
% recovery after a pause



% For now do NOT merge any of the bgap reclassifications into pauses, 
% but keep as separate categories for analysis

dataMat(bGap2Pause,5) = 6; % slow moving shrinks Do NOT put into pause stats- just segregate can later get stats on these!!! 
%dataMat(bGap2Pause,:) = abs(dataMat(bGap2Pause,:));
dataMat(bGapVeryFast,5) = 7;
dataMat(bGapShortTime,5) = 8; % shrinks on order of latency time Do NOT reclassify just segregate
dataMatMerge = dataMat;


% Calculate Fraction Segregated
% fraction of bgaps that are changed to 'pause', because their speeds are
% slower than the cutoff for fgap pauses
if isempty(bGapsIdx)
    percentBgapsSlow=NaN;
else
    percentBgapsSlow=100*length(bGap2Pause)/length(bGapsIdx);
end

% Calculate Fraction of bgaps on the order of comet latency
if isempty(bGapsIdx) 
    percentBgapsShortTime = NaN;
    
else 
    percentBgapsShortTime = 100*length(bGapShortTime)/length(bGapsIdx); 
end 

% Calulate Fration of bgaps that when corrected are VERY FAST
if isempty(bGapsIdx)
    percentBgapsVeryFast = NaN;
else 
    percentBgapsVeryFast = 100*length(bGapVeryFast)/length(bGapsIdx);
end 




elseif (unimodalReclass == 1 && nonsymmetric ~=1 ) 
     
      tracksBGapShortTime = NaN;
      tracksBGapVeryFast = NaN; 
      
       % If new bgap Reclassification Scheme NOT chosen: Reclassify bGaps to pauses based on symmetry assumption: ie that
       % the maximum pause velocity determined from above unimodal thresholding
       % defines the maximum bgap velocity that can be defined as a pause
       % (I personally this this assumption is not valid-as the effect of the latency of 
       % comet formation on the velocity of a fgap and bgap is not the same...  in the paper 
       % it supposedly corresponded well with the gfp tubulin data, perhaps it depends on time between sampling and if 
       % you can consider this effect negligible.  However if this is the
       % case the pause velocities should be very near zero (unless there
       % is significant undetected growths due to the comet moving in and
       % out of the plane). With a sampling rate of 0.75 sec we see a major
       % number of shrinkage events being reclassified at pauses. 
       % (Comment by MB: 03/16/11)
    
        % reclassify based on symmetry with fgaps: a true shrinkage is one that has a greater velocity 
        % observed velocity than the max pause speed (as determined by unimodal thresholding of all fgap 
        % velocities - all other shrinkage events are considered pauses
        bGapsIdx = find(dataMat(:,5) == 3);
        bGap2PauseReclass = find(dataMatReclass(:,5) == 3 & abs(dataMatReclass(:,4)) < cutoffValueFGap);% find the index of those bgaps that are likely pauses
        
        tracksBGapSlow = dataMatReclass(bGap2PauseReclass,1);
        
        % reclassification simply marked (no consolidation)
        dataMatReclass(bGap2PauseReclass,5) = 6;  
        dataMatReclass(bGap2PauseReclass,:) = abs(dataMatReclass(bGap2PauseReclass,:)); % now a pause so change neg to pos value 
        
        % reclassification consolidated bgaps now officially a pause
        % (subTrack ID  in col5 = 2) input for below 
       %bGap2Pause = find(dataMat(:,5) ==3 & abs(dataMat(:,4)) < cutoffValueFGap);
       %dataMat(bGap2Pause,5) = 2; % NOTE: might NOT want to not include in fgap speed calculation when comparing these
       %dataMat(bGap2Pause,:) = abs(dataMat(bGap2Pause,:)); % now a pause so change neg to pos value (valid for interpreting stats of fgap speeds and displacements?)
       
       %Rename for Output 
       dataMatMerge = dataMat;
       
       % Calculate Fraction Reclassified 
       % fraction of bgaps that are changed to pause, because their speeds are
       % slower than the cutoff for fgap pauses
        if isempty(bGapsIdx)
            percentBgapsReclass=NaN;
        else
            percentBgapsReclass=100*length(bGap2Pause)/length(bGapsIdx);
        end
    
    else 
        
end 

%% Remove Growths Initiated in First Frame or Ending in Last Frame From Stats
% do this so one does not bias growth lifetime/displacement data

subIdx2rem=[];
% get index of growth and following fgap or bgap (if it exists) that
% begin in the first frame

sF = projData.detectionFrameRange(1,1);
eF = projData.detectionFrameRange(1,2);

%sF=min(dataMat(:,2));

% compound track IDs of all subtracks with nminimum starting frame number

fullIdx2rem=unique(dataMat(dataMat(:,2)==sF,1)); 
for iTr=1:length(fullIdx2rem)
    subIdx=find(dataMat(:,1)==fullIdx2rem(iTr));
    if (length(subIdx)>1) && (dataMat(subIdx(2),5) > 1) % if there is a forward backward gap linked to growth in first frame
        subIdx2rem=[subIdx2rem; subIdx(1:2)]; % Don't remove entire compound track only the fgap or bgap it's linked to
    else
        subIdx2rem=[subIdx2rem; subIdx(1)];
    end
end

% get index of growth and preceeding fgap or bgap
% (if it exists) that end in the last frame

%eF=max(dataMat(:,3));
fullIdx2rem=unique(dataMat(dataMat(:,3)==eF,1));
for iTr=1:length(fullIdx2rem)
    subIdx=find(dataMat(:,1)==fullIdx2rem(iTr));
    if (length(subIdx)>1) && (dataMat(subIdx(end-1),5) > 1)
        subIdx2rem=[subIdx2rem; subIdx(end-1:end)]; % take out the last two of list  
    else
        subIdx2rem=[subIdx2rem; subIdx(end)];
    end
end
% remove both classes for statistics
dataMat(subIdx2rem,:)=[];



dataMatCrpSecMic=abs(dataMat);


