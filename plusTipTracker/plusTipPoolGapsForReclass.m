function [projData ] = plusTipPoolGapsForReclass( groupList, meta2Use, saveCopy, metaOldNewName, makeHistogram )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%% 




%% Choices for Input
bgapUniModeThreshNoCorrect = 0;

bgapUniModeThreshCorrect = 0;

bgapReclassFluctRadius = 0;

onlyFgap2GrowthReclass = 1;

useFirstInList = 1;
%%

    
projGroupName=groupList(:,1);
projGroupDir= groupList(:,2);

% fix the names if there are spaces or hyphens and append prefix 'grp'
projGroupName=cellfun(@(x) strrep(x,'-','_'),projGroupName,'uniformoutput',0);
projGroupName=cellfun(@(x) strrep(x,' ','_'),projGroupName,'uniformoutput',0);
projGroupName=cellfun(@(x) ['grp_' x],projGroupName,'uniformoutput',0);


% count unique groups and keep them in order of the original list
[btwGrpNames,m,projGroupIdx] = unique(projGroupName);
[b,idx]=sort(m);
btwGrpNames=btwGrpNames(idx);

groupNames = cell(length(btwGrpNames),1);
for iGroup = 1:length(btwGrpNames)
    groupNames(iGroup,1) = btwGrpNames(iGroup,1);
end 




  

for iGroup = 1:length(btwGrpNames)
    
 

tempIdx= strmatch(btwGrpNames(iGroup),projGroupName,'exact');

if (iGroup > 1 && useFirstInList == 1); % do nothing already calculated cut-offs 

else % recalculate cut-off based on each group
    
for iProj = 1:length(tempIdx)
    
    temp = load([projGroupDir{tempIdx(iProj)} filesep meta2Use filesep 'projData']); 


    projData = temp.projData;

    dataMat = projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix;

    if iProj == 1  
    dataMatPooled = dataMat;
        else 
    dataMatPooled = [dataMatPooled ; dataMat] ;
    end 

end % for iProj


 % perform unimodal thresholding using all fgap speeds to obtain 
    % maximum fgap speed corresponding to a pause event (above this thresh
    % are likely undetected growth events)

fgapIdx = find(dataMatPooled(:,5) == 2 | dataMatPooled(:,5) == 5);



fgapSpeeds = dataMatPooled(fgapIdx,4);

%gapSpeed95th = prctile(fgapSpeeds,95);


%fgapSpeeds = fgapSpeeds(fgapSpeeds<fgapSpeed95th);

[cutoffIdx, cutoffValueFGap, sp, axesH,maxBinValue] = cutFirstHistMode(fgapSpeeds,makeHistogram);
     
   
    cutOffValueFGap_VelMicPerMin = cutoffValueFGap;
    maxBinValue_VelMicPerMin = maxBinValue;
    numFgapsTotal = length(fgapIdx);
    fgapReclassScheme = 'Unimodal Thresholding: Pooled Data';
    
    if makeHistogram == 1
       [pathup1 dummy1 dummy2 dummy3] =  getFilenameBody(projGroupDir{tempIdx(iProj)});
       [pathup2 name dummy2 dummy3] = getFilenameBody(pathup1);
       
       uniModalFigureDir = [pathup2 filesep 'PoolUnimodalThreshFigures'];  
       
       if ~isdir(uniModalFigureDir)     
       mkdir(uniModalFigureDir);
       end
       
        saveas(gcf,[uniModalFigureDir filesep name 'UnimodalThresh.fig']); 
        close(gcf);
    end
%% BGap2Pause No Reclass
if onlyFgap2GrowthReclass == 1 % 
    bgapThresh = 0;
    bgapReclassScheme = 'No Bgap2Pause Reclass';
    cutOffValueBGap_VelMicPerMin = 0;
end 
%%  Bgap Reclassification Using Unimodal Fgap Thresh

if bgapUniModeThreshNoCorrect == 1
    
    bgapThresh = cutoffValueFGap;
    bgapReclassScheme = 'Unimodal Fgap Thresh- No Correct: Pooled Data'; 
    cutOffValueBGap_VelMicPerMin = cutoffValueFGap;
    
end

%%  Bgap Reclassification Using "Corrected" Unimodal Fgap Thesh

if bgapUniModeThreshCorrect == 1 
    
    widthDist = cutoffValueFGap - maxBinValue;
    
    
   
    if widthDist < maxBinValue % distribution does not go into the negative values
        bgapThresh = 0;
    else 
        bgapThresh = widthDist - maxBinValue;
    end 
    
    bgapReclassScheme = 'Unimodal Fgap Thresh- Correct For Comet Latency: Pooled Data';
    cutOffValueBGap_VelMicPerMin = bgapThresh;
    maxBinValue_VelMicPerMin = maxBinValue;
    fgapDisbWidth_VelMicPerMin = widthDist;

                                                        
end 

end % if 

%% Using Above Thresholds Correct Individual Projects


for iProj = 1:length(tempIdx)
    
    temp = load([projGroupDir{tempIdx(iProj)} filesep meta2Use filesep 'projData']); 
    projData = temp.projData;
    dataMat =  projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix;
    
   
    fgapIdx = find(dataMat(:,5) == 2 | dataMat(:,5) == 5);
    bgapAllIdx = find(dataMat(:,5) == 3 | dataMat(:,5) == 6);
    
    
    dataMat(fgapIdx,5) = 2;  % make sure to undo any previous reclass scheme 
    dataMat(bgapAllIdx,5) = 3; 
    
    
    growthFgapIdx = find((dataMat(:,5) == 2 | dataMat(:,5) == 5) & abs(dataMat(:,4)) > cutoffValueFGap);
    dataMatReclass = dataMat;
    dataMatReclass(growthFgapIdx,5) =  5; 
    
if bgapReclassFluctRadius == 1
         bgap2pauseIdx = find((dataMat(:,5) == 3 & abs(dataMat(:,7)) < projData.trackingParameters.fluctRadius)); 
     else 
    
   
    bgap2pauseIdx = find((dataMat(:,5) == 3 | dataMat(:,5) == 6) & abs(dataMat(:,4)) < bgapThresh);
    dataMatReclass(bgap2pauseIdx,5) = 6;
    dataMat(bgap2pauseIdx,5) = 2;
    projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix = dataMatReclass;    
end 

%% SAVE INFO REGARDING THE RECLASS SCHEME EMPLOYED %%
 if useFirstInList == 1
     projData.fgapReclassScheme = [fgapReclassScheme, ': Use ', btwGrpNames(1), ' To Set Thresh '];
     projData.bgapReclassScheme = [bgapReclassScheme, ': Use ' btwGrpNames(1), ' To Set Thresh '];
 
 else 
    projData.fgapReclassScheme = fgapReclassScheme;
    projData.bgapReclassScheme = bgapReclassScheme;
 end
 
 projData.cutOffValueFGap_VelMicPerMin = cutOffValueFGap_VelMicPerMin;
 projData.numFGapsUniModalHist = numFgapsTotal;
 
 if (bgapUniModeThreshCorrect == 1 || bgapUniModeThreshNoCorrect ==1 ||onlyFgap2GrowthReclass == 1 ) 
 projData.cutOffValueBGap_VelMicPerMin= cutOffValueBGap_VelMicPerMin;
 end
 
 if bgapUniModeThreshCorrect == 1
 projData.maxBinValue_VelMicPerMin = maxBinValue_VelMicPerMin;
 projData.fgapDisbWidth_VelMicPerMin = fgapDisbWidth_VelMicPerMin;
 end 
     
  
 projData.tracksWithFGapPause = dataMatReclass(dataMatReclass(:,5) == 2,1); 
 projData.tracksWithFGapGrowth = dataMatReclass(growthFgapIdx,1);
 projData.tracksWithBGapPause = dataMatReclass(bgap2pauseIdx,1);
 projData.tracksWithBGapShrink = dataMatReclass(dataMatReclass(:,5) == 3,1);
 
 
 
%%  


% MERGE RECLASSIFICATIONS (ie subTrack ID 5 --> 1 need to incorporate reclassified pauses in growth stats)
% merging will be performed regardless of reclassification scheme unless we
% are pooling data
    
    
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
    
% Conversions and Calculation of Comet Latency
% calculate necessary avg values from above merged data (where all 
% reclassified pauses have been merged with preceding growth 
% subtrack

% perform lifetime and displacement unit conversions
dataMat(:,6)=dataMat(:,6).* projData.secPerFrame; % convert lifetimes to seconds
dataMat(:,7)=dataMat(:,7).*(projData.pixSizeNm/1000); % convert displacements to microns

%% Remove Growths Initiated in First Frame or Ending in Last Frame From Stats
% do this so one does not bias growth lifetime/displacement data (might not
% be what we want for 

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
dataMatCrpSecMic=dataMat; % NOTE: Kathyrn makes these all absolute values 
% I think for stats it is better to keep sign (MB) 

%% Calculate Stats from the matrix where the tracks at the beginning and end have been removed
[projData,projData.M] = plusTipDynamParam(dataMatCrpSecMic,projData,0,0);

MCell = num2cell(projData.M);
titles = cell(1,9);
titles{1,1} = 'growth speed (microns/min)';
titles{1,2} = 'fgap speed (microns/min)';
titles{1,3} = 'bgap speed (microns/min)';
titles{1,4} = 'growth lifetimes (sec)';
titles{1,5} = 'fgap lifetimes (sec)';
titles{1,6} = 'bgap lifetimes (sec)';
titles{1,7} = 'growth displacements (microns)';
titles{1,8} = 'fgap displacements (microns)';
titles{1,9} = 'bgap displacementes (microns)';


projData.MCell = [titles;MCell];

%avgVelGrowth = mean(dataMat(dataMat(:,5) == 1,4));

%avgDispPause = mean(dataMat(dataMat(:,5) == 2,7)); % in microns
%absAvgDispPause = mean(abs(dataMat(dataMat(:,5) == 2,7))); % in microns consider bgaps positive
%avgDispPauseNoBgapReclass= mean(dataMatReclass(dataMatReclass(:,5) == 2, 7)).*projData.pixSizeNm/1000;


% calc avg latency of comet formation  (in min)
%avgLat = avgDispPause/avgVelGrowth; % in minutes
%projData.avgCometLatSec = avgLat*60; % in seconds

%absAvgLat = absAvgDispPause/avgVelGrowth;
%projData.avgCometLatSecAbs = absAvgLat*60;

%avgLatNoBgapReclass = avgDispPauseNoBgapReclass/avgVelGrowth; % in minutes
%projData.avgCometLatSecNoBgapReclass = avgLatNoBgapReclass*60;


%Calculate Fraction of SubTracks Reclassified 
% fraction of bgaps that are changed to pause, because their speeds are
% slower than the cutoff for fgap pauses
    
if isempty(bgapAllIdx)
    projData.percentBgapsReclass=NaN;
else
    projData.percentBgapsReclass=100*length(bgap2pauseIdx)/length(bgapAllIdx);
end

% fraction of fgaps that were consolidated into growth, since their speeds
% were more than 50% of the growth speed just prior to the gap
if isempty(fgapIdx)
    projData.percentFgapsReclass=NaN;
else
    projData.percentFgapsReclass=100*length(growthFgapIdx)/length(fgapIdx);
end
dirNameNew = [projGroupDir{tempIdx(iProj)} filesep 'meta'];
if saveCopy  == 1
% save a copy of the meta Dir Before Pooling
dirBeforeRename = [ projGroupDir{tempIdx(iProj)} filesep meta2Use];
dirAfterRename = [projGroupDir{tempIdx(iProj)} filesep metaOldNewName];

if isdir(dirAfterRename)
    rmdir(dirAfterRename,'s');
end 

movefile(dirBeforeRename, dirAfterRename);
else % don't save a copy 
   % rmdir(dirBeforeRename)
end 

mkdir(dirNameNew);

% rewrite the projData
save([projGroupDir{tempIdx(iProj)} filesep 'meta' filesep 'projData'],'projData')    

end % for iProj
end % for iGroup
end 
