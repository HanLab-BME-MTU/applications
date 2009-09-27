function [dataMatMerge,dataMatReclass,dataMatCrpMinMic,percentFgapsReclass,percentBgapsReclass]=plusTipMergeSubtracks(projData,dataMat)
% plusTipMergeSubtracks merges growth fgaps with the flanking growth phases
%
% SYNOPSIS  : [dataMatMerge,dataMatReclass,percentFgapsReclass]=...
%                   plusTipMergeSubtracks(projData,dataMat)
%
% INPUT
% projData  : structure containing frame rate, pix size info, etc. 
%             if dataMatrix isn't given, then projData should have the
%             field nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix
%             which contains reclassified fgaps.
% dataMat   : mainly to be used in plusTipPostTracking, where it is the
%             data matrix prior to reclassification.
%             
% OUTPUT
% dataMatMerge        : matrix where fgaps that should be reclassified (as
%                       growth, because their speeds are
%                       >=50% the speed at the end of the growth phase
%                       prior) are consolidated with flanking growth
%                       subtracks
% dataMatReclass      : matrix where fgaps that should be reclassified are
%                       not consolidated but the track type is changed to 5
% dataMatCrpMinMic    : like dataMatMerge except the growth phases and
%                       linked fgaps/bgaps from the first and last frames
%                       of the movie have been removed. also, column 6
%                       represents lifetime in minutes and 7 represents
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


if nargin<2
    % dataMat is the output matrix in projData, where the gaps to
    % conolidate are labeled as trackType=5
    dataMat=projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix;
    
    % the set of all fgaps in this case have trackType of 2 or 5
    fgapIdx=find(dataMat(:,5)==2 | dataMat(:,5)==5);
    
    % these are the fgaps to consolidate
    growthFgapIdx=find(dataMat(:,5)==5);
else
    % dataMat is an input, we need to find where the growthFgaps are

    % look at fgaps - if their speed << growth phase just prior, then it's a
    % true pause.  it it is greater than 90% growth speed, reclassify it as
    % continuation of growth.
    fgapIdx=find(dataMat(:,5)==2);
    
    beforeFgapIdx=fgapIdx-1;
    % get speed of growth phase at last 2-3 time pts (depending on how long
    % growth phase lasted) just prior to pause
    eF=dataMat(beforeFgapIdx,3);
    beforeFgapSpeed=zeros(length(fgapIdx),1);
    for iGap=1:length(fgapIdx)
        idxRange=[max(1,eF(iGap)-3):eF(iGap)-1];
        beforeFgapSpeed(iGap)=nanmean(projData.frame2frameVel_micPerMin(dataMat(beforeFgapIdx(iGap),1),idxRange));
    end
    % these are the fgaps to consolidate
    growthFgapIdx=fgapIdx(dataMat(fgapIdx,4)>0.7.*beforeFgapSpeed);
end


% same as input, except growth fgaps now have 5 for trackType
dataMatReclass=dataMat;
dataMatReclass(growthFgapIdx,5)=5;

% get 95th-percentile of true pause speeds
fgapMaxSpeed=prctile(dataMatReclass(dataMatReclass(:,5)==2,4),95);

% get indices of bgaps where the speeds is less than the 95% cutoff
bgapAllIdx=find(dataMatReclass(:,5)==3);
bgap2pauseIdx=find(dataMatReclass(:,5)==3 & abs(dataMatReclass(:,4))<fgapMaxSpeed);

% bspeeds=abs(dataMatReclass(bgapAllIdx,4));
% [cutoffIndex, cutoffValue] = cutFirstHistMode(bspeeds,1);

dataMat(bgap2pauseIdx,5)=2; % reassign dataMat to have type 2 
dataMat(bgap2pauseIdx,:)=abs(dataMat(bgap2pauseIdx,:)); % get abs value since now pause
dataMatReclass(bgap2pauseIdx,5)=6; % reassign reclassified matrix to have new type of 6

% fraction of bgaps that are changed to pause, because their speeds are
% slower than the cutoff for fgap pauses
percentBgapsReclass=100*length(bgap2pauseIdx)/length(bgapAllIdx);

% fraction of fgaps that were consolidated into growth, since their speeds
% were more than 50% of the growth speed just prior to the gap
percentFgapsReclass=100*length(growthFgapIdx)/length(fgapIdx);


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

dataMatMerge=dataMat;



dataMat(:,6)=dataMat(:,6).* projData.secPerFrame; % convert lifetimes to seconds
dataMat(:,7)=dataMat(:,7).*(projData.pixSizeNm/1000); % convert displacements to microns

subIdx2rem=[];
% get index of growth and following fgap or bgap (if it exists) that
% begin in the first frame
sF=min(dataMat(:,2));
fullIdx2rem=unique(dataMat(dataMat(:,2)==sF,1));
for iTr=1:length(fullIdx2rem)
    subIdx=find(dataMat(:,1)==fullIdx2rem(iTr));
    if length(subIdx)>1
        subIdx2rem=[subIdx2rem; subIdx(1:2)];
    else
        subIdx2rem=[subIdx2rem; subIdx(1)];
    end
end
% get index of growth and preceeding fgap or bgap
% (if it exists) that end in the last frame
eF=max(dataMat(:,3));
fullIdx2rem=unique(dataMat(dataMat(:,3)==eF,1));
for iTr=1:length(fullIdx2rem)
    subIdx=find(dataMat(:,1)==fullIdx2rem(iTr));
    if length(subIdx)>1
        subIdx2rem=[subIdx2rem; subIdx(end-1:end)];
    else
        subIdx2rem=[subIdx2rem; subIdx(end)];
    end
end
% remove both classes for statistics
dataMat(subIdx2rem,:)=[];

dataMatCrpMinMic=abs(dataMat);


