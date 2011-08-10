function [projData,M]=plusTipDynamParam(dataMatCrpSecMic,projData,fromPoolGroupData,subRoiAnalysis)
% plusTipDynamParam: generic function for calculating dynamics parameters
%
% SYNOPSIS: [projData.stats,M]=plusTipDynamParam(dataMatCrpSecMic)
%
% INPUT:
% dataMatCrpSecMic : matrix produced by plusTipMergeSubtracks, or a
%                    concatenated matrix from multiple movies
%
% OUTPUT:
% projData.stats : structure containing parameters (see plusTipPostTracking for
%         list) based on the subset of tracks starting after the first
%         frame and ending before the last frame
% M     : n x 9 matrix, where the columns are
%           1. growth speed (microns/min)
%           2. fgap speed
%           3. bgap speed
%           4. growth lifetimes (sec)
%           5. fgap lifetimes
%           6. bgap lifetimes
%           7. growth displacements (microns)
%           8. fgap displacements
%           9. bgap displacements
%
% functions that call this one:
% plusTipPostTracking
% plusTipPoolGroupData


%subRoiAnalysis = 0; 
%subRoiAnalysis = 0; 
% Put a Copy of the Data Mat Where Those Tracks Beginning in the First
% Frame or Ending in the Last Frame Have Been Removed
projData.dataMat_FOR_STATS = dataMatCrpSecMic;  

% growth speed, lifetime, displacement
gIdx=find(dataMatCrpSecMic(:,5)==1);
gs=dataMatCrpSecMic(gIdx,4);
gl=dataMatCrpSecMic(gIdx,6);
gd=dataMatCrpSecMic(gIdx,7);

% fgap speed, lifetime, displacement
fIdx=find(dataMatCrpSecMic(:,5)==2);
fs=dataMatCrpSecMic(fIdx,4);
fl=dataMatCrpSecMic(fIdx,6);
fd=dataMatCrpSecMic(fIdx,7);

% bgap speed, lifetime, displacement
bIdx=find(dataMatCrpSecMic(:,5)==3);
bs=dataMatCrpSecMic(bIdx,4);
bl=dataMatCrpSecMic(bIdx,6);
bd=dataMatCrpSecMic(bIdx,7);

% put populations into a matrix backfilled with NaNs
M=nan(max([length(gs) length(fs) length(bs)]),9);
M(1:length(gs),1)=gs;
M(1:length(fs),2)=fs;
M(1:length(bs),3)=bs;
M(1:length(gl),4)=gl;
M(1:length(fl),5)=fl;
M(1:length(bl),6)=bl;
M(1:length(gd),7)=gd;
M(1:length(fd),8)=fd;
M(1:length(bd),9)=bd;

%% PARAMETERS RELATED TO COMET DENSITY 

projData.stats.medNNdistWithinFrameMic = projData.medNNdistWithinFramePix*projData.pixSizeNm/1000;


%% PARAMETERS RELATED TO GROWTH

if isempty(gIdx)
    projData.stats.nGrowths=0;
    % median/mean for growth speeds (microns/minute)
    projData.stats.growth_speed_median = NaN;
    projData.stats.growth_speed_mean = NaN;
    projData.stats.growth_speed_std = NaN;
    projData.stats.growth_speed_mean_robust = NaN;
    
    % median/mean for growth lifetime (sec)
    projData.stats.growth_lifetime_median = NaN;
    projData.stats.growth_lifetime_mean = NaN;
    projData.stats.growth_lifetime_std = NaN;
    projData.stats.growth_lifetime_mean_robust = NaN;
    
    % median/mean for growth displacement (microns)
    projData.stats.growth_length_median = NaN;
    projData.stats.growth_length_mean = NaN;
    projData.stats.growth_length_std = NaN;
    projData.stats.growth_length_mean_robust = NaN;
    
else
    projData.stats.nGrowths=length(gIdx);
    
    % median/mean for growth speeds (microns/minute)
    projData.stats.growth_speed_median = median(gs);
    projData.stats.growth_speed_mean = mean(gs) ;
    projData.stats.growth_speed_std  = std(gs);
    projData.stats.growth_speed_mean_robust  = trimmean(gs,10);
    
    % median/mean for growth lifetime (sec)  
    projData.stats.growth_lifetime_median = median(gl);
    projData.stats.growth_lifetime_mean = mean(gl);
    projData.stats.growth_lifetime_std = std(gl);
    projData.stats.growth_lifetime_mean_robust  = trimmean(gl,10);
    
    % median/mean for growth displacement (microns)
    projData.stats.growth_length_median = median(gd);
    projData.stats.growth_length_mean = mean(gd);
    projData.stats.growth_length_std = std(gd);
    projData.stats.growth_length_mean_robust = trimmean(gd,10);
    
end

projData.stats.percentFgapsReclass = projData.percentFgapsReclass;

%% PARAMETERS RELATED TO FGAPS 

if isempty(fIdx)
    projData.stats.nFgaps=0;
    projData.stats.nFgaps2nGrowths = 0;
    % median/mean for fgap speeds (microns/minute)
    projData.stats.fgap_speed_median = NaN;
    projData.stats.fgap_speed_mean = NaN;
    projData.stats.fgap_speed_std = NaN;
    projData.stats.fgap_speed_mean_robust = NaN;
    
    % median/mean for fgap lifetime (sec)
    projData.stats.fgap_lifetime_median = NaN;
    projData.stats.fgap_lifetime_mean = NaN;
    projData.stats.fgap_lifetime_std = NaN;
    projData.stats.fgap_lifetime_mean_robust = NaN;
    
    % median/mean for fgap displacement (microns)
    projData.stats.fgap_length_median = NaN;
    projData.stats.fgap_length_mean = NaN ;
    projData.stats.fgap_length_std = NaN;
    projData.stats.fgap_length_mean_robust = NaN;
    
    % fgap frequency is the inverse of the average growth time (sec) or
    % displacement (microns) prior to fgap
    %projData.stats.fgap_freq_time=NaN;
    %projData.stats.fgap_freq_length=NaN;
    
    % 
    projData.stats.fgap_VelGrowthBefore_MicPerMin_mean = NaN;
    projData.stats.fgap_VelGrowthBefore_MicPerMin_SE =  NaN;
    projData.stats.fgap_LifetimeGrowthBefore_Sec_mean= NaN;
    projData.stats.fgap_LifetimeGrowthBefore_Sec_SE = NaN;
    projData.stats.fgap_LengthGrowthBefore_Mic_mean = NaN;
    projData.stats.fgap_LengthGrowthBefore_Mic_SE = NaN;

else
    projData.stats.nFgaps=length(fIdx);
    projData.stats.nFgaps2nGrowths = projData.stats.nFgaps/projData.stats.nGrowths;
    % median/mean for fgap speeds (microns/minute)
    projData.stats.fgap_speed_median = median(fs);
    projData.stats.fgap_speed_mean = mean(fs);
    projData.stats.fgap_speed_std = std(fs);
    projData.stats.fgap_speed_mean_robust = trimmean(fs,10);
    
    % median/mean for fgap lifetime (sec)
    projData.stats.fgap_lifetime_median = median(fl);
    projData.stats.fgap_lifetime_mean = mean(fl);
    projData.stats.fgap_lifetime_std = std(fl);
    projData.stats.fgap_lifetime_mean_robust = trimmean(fl,10);
    
    % median/mean for fgap displacement (microns)
    projData.stats.fgap_length_median = median(fd);
    projData.stats.fgap_length_mean = mean(fd);
    projData.stats.fgap_length_std = std(fd);
    projData.stats.fgap_length_mean_robust = trimmean(fd,10);
    
    % fgap frequency is the inverse of the average growth time (sec) or
    % displacement (microns) prior to fgap
  
    beforeFgapIdx=fIdx-1; % get idx of the growth subtrack right before event
    
    %A fix for if an error which may occur if there dataMats from sub-rois
    % are input
    
    if beforeFgapIdx(1) == 0; % This means the first track was a gap event  this should of course not happen with                       
        beforeFgapIdx(1) = []; % simply remove the 0 idx as it is nonsensical and will give error. 
    end 
    
    
    
    %Calculate Average Time and Displacement Parameters for Growth
    %Preceding Fgap 
     
    velocityBeforeFgap = dataMatCrpSecMic(beforeFgapIdx,4);
    lifetimeBeforeFgap = dataMatCrpSecMic(beforeFgapIdx,6);
    lengthBeforeFgap =  dataMatCrpSecMic(beforeFgapIdx,7);
    
    projData.stats.fgap_VelGrowthBefore_MicPerMin_mean = mean(velocityBeforeFgap);
    projData.stats.fgap_VelGrowthBefore_MicPerMin_SE =  std(velocityBeforeFgap)/sqrt(length(velocityBeforeFgap));
    projData.stats.fgap_LifetimeGrowthBefore_Sec_mean= mean(lifetimeBeforeFgap);
    projData.stats.fgap_LifetimeGrowthBefore_Sec_SE = std(lifetimeBeforeFgap)/sqrt(length(lifetimeBeforeFgap));
    projData.stats.fgap_LengthGrowthBefore_Mic_mean = mean(lengthBeforeFgap);
    projData.stats.fgap_LengthGrowthBefore_Mic_SE = std(lengthBeforeFgap)/sqrt(length(lengthBeforeFgap));
    
    %Convert these values to frequency
    %projData.stats.fgap_freq_time=1/mean(dataMatCrpSecMic(beforeFgapIdx,6));
    %projData.stats.fgap_freq_length=1/mean(dataMatCrpSecMic(beforeFgapIdx,7));
    
    %Get the Individual Frequencies Corresponding to the Growth Before EACH
    % pause, the average and the std of these values
    %freq=1./dataMatCrpSecMic(beforeFgapIdx,6);
    %projData.stats.fgap_freq_time_mean_SE=[mean(freq) std(freq)/sqrt(length(freq))];
    %freq=1./dataMatCrpSecMic(beforeFgapIdx,7);
    %projData.stats.fgap_freq_length_mean_SE=[mean(freq) std(freq)/sqrt(length(freq))];
         % all track indices where there is either a forward or backward gap

    
end %isempty

%% PARAMETERS RELATED TO BGAPS

if isempty(bIdx)
    projData.stats.nBgaps=0;
    projData.stats.nBgaps2nGrowths = 0;
    % median/mean for bgap speeds (microns/minute)
    projData.stats.bgap_speed_median = NaN;
    projData.stats.bgap_speed_mean = NaN;
    projData.stats.bgap_speed_std = NaN;
    projData.stats.bgap_speed_mean_robust = NaN;
    
    % median/mean for bgap lifetime (sec)
    projData.stats.bgap_lifetime_median = NaN;
    projData.stats.bgap_lifetime_mean= NaN;
    projData.stats.bgap_lifetime_std = NaN;
    projData.stats.bgap_lifetime_mean_robust = NaN;
    
    % median/mean for bgap displacement (microns)
    projData.stats.bgap_length_median = NaN;
    projData.stats.bgap_length_mean = NaN;
    projData.stats.bgap_length_std = NaN;
    projData.stats.bgap_length_mean_robust = NaN;
    
    % bgap frequency is the inverse of the average growth time (sec) or
    % displacement (microns) prior to bgap
    %projData.stats.bgap_freq_time=NaN;
    %projData.stats.bgap_freq_length=NaN;
    
    projData.stats.bgap_VelGrowthBefore_MicPerMin_mean = NaN;
    projData.stats.bgap_VelGrowthBefore_MicPerMin_SE = NaN;
    projData.stats.bgap_LifetimeGrowthBefore_Sec_mean = NaN ;
    projData.stats.bgap_LifetimeGrowthBefore_Sec_SE = NaN;
    projData.stats.bgap_LengthGrowthBefore_Mic_mean = NaN ;
    projData.stats.bgap_LengthGrowthBefore_Mic_SE = NaN; 
    

else
    projData.stats.nBgaps=length(bIdx);
    projData.stats.nBgaps2nGrowths = projData.stats.nBgaps/projData.stats.nGrowths;
    % median/mean for bgap speeds (microns/minute)
    projData.stats.bgap_speed_median = median(bs);
    projData.stats.bgap_speed_mean = mean(bs);
    projData.stats.bgap_speed_std = std(bs);
    projData.stats.bgap_speed_mean_robust = trimmean(bs,10);
    % median/mean for bgap lifetime (sec)
    projData.stats.bgap_lifetime_median = median(bl);
    projData.stats.bgap_lifetime_mean = mean(bl);
    projData.stats.bgap_lifetime_std = std(bl);
    projData.stats.bgap_lifetime_mean_robust = trimmean(bl,10);
    
    % median/mean for bgap displacement (microns)
    projData.stats.bgap_length_median = median(bd);
    projData.stats.bgap_length_mean = mean(bd);
    projData.stats.bgap_length_std = std(bd);
    projData.stats.bgap_length_mean_robust = trimmean(bd,10);
    
    
    
    %%%%% STAT VALUES PRIOR TO BGAP %%%%
    
    beforeBgapIdx=bIdx-1;
    
    
    if beforeBgapIdx(1) == 0; % This means the first track was a gap event 
                              % this should of course not happen with
                              % normal tracking, but tracks can be sub-divided
                              % by spatial regions and therefore the pause
                              % and not the growth track just before 
                              % might be first in dataMatCrpSecMat
        beforeBgapIdx(1) = []; % simply remove the 0 idx as it is nonsensical and will give error. 
    end 
    
    %Calculate Average Time and Displacement Parameters for Growth
    %Preceding Bgap
    velocityBeforeBgap = dataMatCrpSecMic(beforeBgapIdx,4);
    lifetimeBeforeBgap = dataMatCrpSecMic(beforeBgapIdx,6);
    lengthBeforeBgap =  dataMatCrpSecMic(beforeBgapIdx,7);
    
    projData.stats.bgap_VelGrowthBefore_MicPerMin_mean = mean(velocityBeforeBgap);
    projData.stats.bgap_VelGrowthBefore_MicPerMin_SE = std(velocityBeforeBgap)/sqrt(length(velocityBeforeBgap));
    projData.stats.bgap_LifetimeGrowthBefore_Sec_mean = mean(lifetimeBeforeBgap) ;
    projData.stats.bgap_LifetimeGrowthBefore_Sec_SE = std(lifetimeBeforeBgap)/sqrt(length(lifetimeBeforeBgap));
    projData.stats.bgap_LengthGrowthBefore_Mic_mean = mean(lengthBeforeBgap) ;
    projData.stats.bgap_LengthGrowthBefore_Mic_SE = std(lengthBeforeBgap)/sqrt(length(lengthBeforeBgap));
   
    
    %Convert these Average Values to frequencies
    %projData.stats.bgap_freq_time=1/mean(lifetimeBeforeBgap);
    %projData.stats.bgap_freq_length=1/mean(lengthBeforeBgap);
    
    % Find Freqency Values for each subtrack and find avg and std
   % freq=1./dataMatCrpSecMic(beforeBgapIdx,6);
   % projData.stats.bgap_freq_time_mean_SE=[mean(freq) std(freq)/sqrt(length(freq))];
    
    %freq =1./dataMatCrpSecMic(beforeBgapIdx,7);
    %projData.stats.bgap_freq_length_mean_SE=[mean(freq) std(freq)/sqrt(length(freq))];
  
end % isempty

%% PARAMETERS OF GROWTH SUBTRACKS PRECEDING TERMINAL EVENT
    
  if subRoiAnalysis == ~1  % if calling from original analysis  
       % plusTipPoolGroupData and plusTipSubRoiExtractTracks will use value
       % of 1 


    % Create data mat that has only growth before term events
    gapIdx = find(dataMatCrpSecMic(:,5)~=1); % find all gap subtracks 
    
   


    beforeGapIdx = gapIdx(gapIdx~=1)-1; % find growth subtracks before these gaps
    % SB: removed test below as it crashes when gapIdx is empty
%     if beforeGapIdx(1) == 0; % This means the first track was a gap event 
%                               % this should of course not happen with
%                               % normal tracking, but tracks can be sub-divided
%                               % by spatial regions and therefore the pause
%                               % and not the growth track just before 
%                               % might be first in dataMatCrpSecMat
%         beforeGapIdx(1) = []; % simply remove the 0 idx as it is nonsensical and will give error. 
%     end 
    
    % Note these values might not be correct for  sub roi manipulations
    % Have to go back and check 
    termGrowthOnly = dataMatCrpSecMic; % initiate for manipulation 
    termGrowthOnly(beforeGapIdx,:) = []; %  remove all growth subtracks before fgaps
    gapIdxNew = termGrowthOnly(:,5) ~= 1; % get the new Idx for gap subtracks 
    termGrowthOnly(gapIdxNew,:) = []; % remove all gap substracks: now should have only growth subtracks before a terminal event
    
    % collect the parameter under question for all subtracks
    velocityBeforeTermEvent = termGrowthOnly(:,4);
    lifetimeBeforeTermEvent = termGrowthOnly(:,6);
    lengthBeforeTermEvent = termGrowthOnly(:,7); 
    
    %%% AVERAGE VELOCITY, TIME, AND DISPLACEMENT PRECEDING TERMINAL EVENT%%% 
    
    projData.stats.term_VelGrowthBefore_MicPerMin_mean = mean(velocityBeforeTermEvent);
    projData.stats.term_VelGrowthBefore_MicPerMin_SE = std(velocityBeforeTermEvent)/sqrt(length(velocityBeforeTermEvent)) ;
    projData.stats.term_LifetimeGrowthBefore_Sec_mean = mean(lifetimeBeforeTermEvent);
    projData.stats.term_LifetimeGrowthBefore_Sec_SE = std(lifetimeBeforeTermEvent)/sqrt(length(lifetimeBeforeTermEvent));
    projData.stats.term_LengthGrowthBefore_Mic_mean = mean(lengthBeforeTermEvent); 
    projData.stats.term_LengthGrowthBefore_Mic_SE = std(lengthBeforeTermEvent)/sqrt(length(lengthBeforeTermEvent));
    
    %%% FREQUENCIES (INVERSES OF AVG TIME AND DISPLACEMENT %%%%
  
    projData.stats.term_freq_time=1/mean(termGrowthOnly(:,6));
    projData.stats.term_freq_length=1/mean(termGrowthOnly(:,7));
    
    %%% INDIVIDUAL FREQUENCIES %%%
    
    %freq=1./termGrowthOnly(:,6);
   % projData.stats.fgap_freq_time_mean_SE=[mean(freq) std(freq)/sqrt(length(freq))];
    %freq=1./termGrowthOnly(:,7);
    %projData.stats.fgap_freq_length_mean_SE=[mean(freq) std(freq)/sqrt(length(freq))];
   
    %%% RATIOS OF GROWTH VELOCITY, GROWTH LIFETIME, OR GROWTH LENGTH JUST 
    % BEFRORE AN FGAP/BGAP/TERM EVENT

    % Velocity Ratios 
    projData.stats.ratio_preFgapVel2preTermVel = projData.stats.fgap_VelGrowthBefore_MicPerMin_mean/projData.stats.term_VelGrowthBefore_MicPerMin_mean;
    projData.stats.ratio_preBgapVel2preTermVel = projData.stats.bgap_VelGrowthBefore_MicPerMin_mean/projData.stats.term_VelGrowthBefore_MicPerMin_mean;
    projData.stats.ratio_preFgapVel2preBgapVel= projData.stats.fgap_VelGrowthBefore_MicPerMin_mean/projData.stats.bgap_VelGrowthBefore_MicPerMin_mean;
    
    % Lifetime Ratios
    projData.stats.ratio_preFgapLife2preTermLife  = projData.stats.fgap_LifetimeGrowthBefore_Sec_mean/projData.stats.term_LifetimeGrowthBefore_Sec_mean;
    projData.stats.ratio_preBgapLife2preTermLife= projData.stats.bgap_LifetimeGrowthBefore_Sec_mean/projData.stats.term_LifetimeGrowthBefore_Sec_mean;
   projData.stats.ratio_preFgapLife2preBgapLife  = projData.stats.fgap_LifetimeGrowthBefore_Sec_mean/projData.stats.bgap_LifetimeGrowthBefore_Sec_mean;
    
    % Displacment (length) ratios
    projData.stats.ratio_preFgapDisp2preTermDisp = projData.stats.fgap_LengthGrowthBefore_Mic_mean/projData.stats.term_LengthGrowthBefore_Mic_mean;
    projData.stats.ratio_preFgapDisp2preTermDisp = projData.stats.bgap_LengthGrowthBefore_Mic_mean/projData.stats.term_LengthGrowthBefore_Mic_mean;
    projData.stats.ratio_preFgapDisp2preBgapDisp = projData.stats.fgap_LengthGrowthBefore_Mic_mean/projData.stats.bgap_LengthGrowthBefore_Mic_mean;
    
  else % don't do these calcs because it will just error 
  end 
     
%% MISC PARAMETERS

% percent of time spent in growth, fgap, and bgap
totalTime=sum(gl)+sum(fl)+sum(bl);
if totalTime==0
    projData.stats.percentTimeGrowth=NaN;
    projData.stats.percentTimeFgap  =NaN;
    projData.stats.percentTimeBgap  =NaN;
else
    projData.stats.percentTimeGrowth=100*(sum(gl)/totalTime);
    projData.stats.percentTimeFgap  =100*(sum(fl)/totalTime);
    projData.stats.percentTimeBgap  =100*(sum(bl)/totalTime);
end

if projData.stats.nFgaps+projData.stats.nBgaps==0
    % percent nFgaps/nGaps
    projData.stats.percentGapsForward = NaN;
    % percent nBgaps/nGaps
    projData.stats.percentGapsBackward= NaN;
else
    % percent nFgaps/nGaps
    projData.stats.percentGapsForward = 100*(projData.stats.nFgaps/(projData.stats.nFgaps+projData.stats.nBgaps));
    % percent nBgaps/nGaps
    projData.stats.percentGapsBackward= 100*(projData.stats.nBgaps/(projData.stats.nFgaps+projData.stats.nBgaps));
end


if isempty(gIdx)
    projData.stats.percentGrowthLinkedForward  = NaN;
    projData.stats.percentGrowthLinkedBackward = NaN;
    projData.stats.percentGrowthTerminal       = NaN;

else
    % get rid of the last one if it's a growth event, since there is no "next
    % index" for that one
    nextIdx=gIdx+1;
    if nextIdx(end)>size(dataMatCrpSecMic,1)
        gIdx(end)=[];
        nextIdx(end)=[];
    end
    % get rid of any indices where the next one is an fgap or bgap but it
    % doesn't come from the same track
    rmIdx=find(dataMatCrpSecMic(gIdx,1)~=dataMatCrpSecMic(nextIdx,1) & (dataMatCrpSecMic(nextIdx,5)==2 | dataMatCrpSecMic(nextIdx,5)==3));
    gIdx(rmIdx)=[];
    nextIdx(rmIdx)=[];

    % percent of growths ending in fgap, bgap, or nothing
    f=sum(dataMatCrpSecMic(nextIdx,5)==2);
    b=sum(dataMatCrpSecMic(nextIdx,5)==3);
    u=length(gIdx)-(f+b);

    projData.stats.percentGrowthLinkedForward  = 100*(f/length(gIdx));
    projData.stats.percentGrowthLinkedBackward = 100*(b/length(gIdx));
    projData.stats.percentGrowthTerminal       = 100*(u/length(gIdx));
end





%% Compound Versus Single Track projData.stats
    % all track indices where there is either a forward or backward gap
    tracksWithFgap=unique(dataMatCrpSecMic(fIdx,1));
    tracksWithBgap=unique(dataMatCrpSecMic(bIdx,1));
    tracksWithGap = unique(dataMatCrpSecMic([fIdx ;bIdx],1));

     numGrowthSubTracksAll = length(find(dataMatCrpSecMic(:,5) == 1));
     numTracksTotal = length(unique(dataMatCrpSecMic(:,1)));
     projData.stats.ratio_TotalTracks2NumGrowthSubtracks =  numTracksTotal/numGrowthSubTracksAll;
  
     projData.stats.NumTracksWithfGap2TotalTracks_Per = length(tracksWithFgap)/numTracksTotal*100;
     projData.stats.NumTracksWithbGap2TotalTracks_Per = length(tracksWithBgap)/numTracksTotal*100;
     projData.stats.NumTracksWithGap2TotalTracks_Per =  length(tracksWithGap)/numTracksTotal*100;
     
   % Segregate Tracks That Are Exclusively From a Compound Track
   
   % Get indices of gap subtracks and those tracks before and after a gap
   % Currently Includes Undefined Gaps
   
   % Fix for if have divided the tracks by subRoi and have the first member
   % of dataMatCrpSecMic a pause or a shrink:
   
   %if dataMatCrpSecMic(gIdx(1),5) ~=1 && dataMatCrpSecMic(gIdx(end),5) ~= 1
    %   compIdx = [gapIdx(1); (gapIdx(1)-1)]; % just get the growth subtrack after this first gap
     %  compIdx = [compIdx; gapIdx(2:end) ; (gapIdx(2:end) + 1); (gapIdx(2,end) -1)]; % all the rest are fine
   
   %elseif dataMatCrpSecMic(gIdx(1),5) == 1 && dataMatCrpSecMic(gIdx(end),5) ~= 1
    %   compIdx = [gapIdx(1:end-1); (gapIdx(1:end-1) - 1)]; % just get the growth subtrack before the d
   %    compIdx = [compIdx; gapIdx(end) ; gapIdx(end) +1]; % all the rest are fine 
    
  % elseif dataMatCrpSecMic(gIdx(1),5) == 1 && dataMatCrpSecMic(gIdx(end,5) == 1
    %    compIdx = [gapIdx(1); (gapIdx(1)-1)];
    %   compIdx = [gapIdx(1:end-1) ; (gapIdx(1:end-1)-1)];
    %   compIdx = [gapIdx(1:
       
   %For now during subROI analysis the full compound track may not be 
   % partitioned in the same ROI making the below stats buggy.  
   % until fix simply do not do this analysis for subROI regions. 
   
   if (subRoiAnalysis ~=1 && fromPoolGroupData ~= 1) % if not calling this from original analysis 
       % skip below 
       % plusTipPoolGroupData and plusTipSubRoiExtractTracks 
       % so will skip the below analysis otherwise buggy MB : quick fix Check it later. 
       
   compIdx= [gapIdx ; (gapIdx+1) ; (gapIdx -1)];
   compIdx = unique(sort(compIdx));
   compDataMat = dataMatCrpSecMic(compIdx,:);
  
   
   % Compound Track No Unclassified Gaps : use this one for stats
   uGapIdx = find(compDataMat(:,5) == 4); 
   compIdxUGap = [uGapIdx ; (uGapIdx +1); (uGapIdx-1)];
   compDataMat(compIdxUGap,:) = [];
   
   projData.compDataMat = compDataMat;
   
   
   
   % Segregate Tracks That Are Exclusively From Single Tracks
   singleDataMat = dataMatCrpSecMic;
   singleDataMat(compIdx,:) = [];
   projData.singleDataMat = singleDataMat;
   
   % Number of Compound to Number of Single Tracks 
   numCompTracks = length(unique(compDataMat(:,1)));
   numSingleTracks = length(singleDataMat(:,1));
   
   projData.stats.ratio_Compound2SingleTracks = numCompTracks/numSingleTracks;
   
   % Get Growth Parameters Exclusively From Either the Single or Compound
   % Track
   vc = mean(compDataMat(compDataMat(:,5) ==1,4));
   vs = mean(singleDataMat(singleDataMat(:,5) == 1,4));
 
   
   lc = mean(compDataMat(compDataMat(:,5) == 1,6));
   ls = mean(singleDataMat(singleDataMat(:,5) ==1,6));
  
   
   dc = mean(compDataMat(compDataMat(:,5) == 1,7));
   ds = mean(singleDataMat(singleDataMat(:,5) == 1,7));
  
   
   
   projData.stats.VelGrowthInCompTrack_mean_MicPerMin = vc;  
   projData.stats.VelGrowthInSingleTrack_mean_MicPerMin = vs;
  
   
   projData.stats.LifeGrowthInCompTrack_mean_Sec = lc;
   projData.stats.LifeGrowthInSingleTrack_mean_Sec = ls;
   
   projData.stats.DispGrowthInCompTrack_mean_Mic = dc;
   projData.stats.DispGrowthInSingleTrack_mean_Mic = ds;
   
     
   projData.stats.ratio_VelGrowthComp2Single = vc/vs;
   projData.stats.ratio_LifeGrowthComp2Single = lc/ls;
   projData.stats.ratio_DispGrowthComp2Single = dc/ds;
   
   
       % do not perform this analysis and remove these fields as they
       % represent the original analysis and this can be confusing. 
      projData =  rmfield(projData,'singleDataMat');
      projData =  rmfield(projData,'compDataMat');
      projData.stats =  rmfield(projData.stats,'ratio_Compound2SingleTracks');
      projData.stats =  rmfield(projData.stats,'VelGrowthInCompTrack_mean_MicPerMin');
      projData.stats =  rmfield(projData.stats, 'VelGrowthInSingleTrack_mean_MicPerMin');
      projData.stats =  rmfield(projData.stats, 'LifeGrowthInCompTrack_mean_Sec');
      projData.stats =  rmfield(projData.stats, 'LifeGrowthInSingleTrack_mean_Sec');
      projData.stats =  rmfield(projData.stats, 'DispGrowthInCompTrack_mean_Mic');
      projData.stats =  rmfield(projData.stats, 'DispGrowthInSingleTrack_mean_Mic');
      projData.stats =  rmfield(projData.stats, 'ratio_VelGrowthComp2Single');
      projData.stats =  rmfield(projData.stats, 'ratio_LifeGrowthComp2Single');
      projData.stats =  rmfield(projData.stats, 'ratio_DispGrowthComp2Single');
       
      
      
   end 
       
   % 
   

   
%% Percent Time 
   
% calculate the average percentage of time a MT spends in fgap
idx=unique(tracksWithFgap);
if isempty(idx)
    projData.stats.avgIndivPercentTimeFgap=NaN;
else
    idxCell=mat2cell(idx,ones(length(idx),1),1);
    % sub track indices of the full tracks
    subIdxAll=cellfun(@(x) find(dataMatCrpSecMic(:,1)==x),idxCell,'uniformoutput',0);
    % full track lifetimes and displacements
    ltfAll=cell2mat(cellfun(@(x) sum(dataMatCrpSecMic(x,6)),subIdxAll,'uniformoutput',0));

    % sub track indices of the fgaps
    subIdxFgaps=cellfun(@(x) find(dataMatCrpSecMic(:,1)==x & dataMatCrpSecMic(:,5)==2),idxCell,'uniformoutput',0);
    % fgap lifetimes and displacements
    ltfFgaps=cell2mat(cellfun(@(x) sum(dataMatCrpSecMic(x,6)),subIdxFgaps,'uniformoutput',0));

    % average percent of time spent in fgap for individual MT
    projData.stats.avgIndivPercentTimeFgap=100*mean(ltfFgaps./ltfAll);

    clear idx idxCell subIdxAll
end


% calculate the average percentage of time a MT spends in bgap
idx=unique(tracksWithBgap);
if isempty(idx)
    projData.stats.avgIndivPercentTimeBgap=NaN;
else

    idxCell=mat2cell(idx,ones(length(idx),1),1);
    % sub track indices of the full tracks
    subIdxAll=cellfun(@(x) find(dataMatCrpSecMic(:,1)==x),idxCell,'uniformoutput',0);
    % full track lifetimes and displacements
    ltfAll=cell2mat(cellfun(@(x) sum(dataMatCrpSecMic(x,6)),subIdxAll,'uniformoutput',0));

    % sub track indices of the bgaps
    subIdxBgaps=cellfun(@(x) find(dataMatCrpSecMic(:,1)==x & dataMatCrpSecMic(:,5)==3),idxCell,'uniformoutput',0);
    % bgap lifetimes and displacements
    ltfBgaps=cell2mat(cellfun(@(x) sum(dataMatCrpSecMic(x,6)),subIdxBgaps,'uniformoutput',0));

    % average percent of time spent in bgap for individual MT
    projData.stats.avgIndivPercentTimeBgap=100*mean(ltfBgaps./ltfAll);

    clear idx idxCell subIdxAll
end

%% Dynamicity
% calculate dynamicity (mic/min)
idx=unique([tracksWithFgap; tracksWithBgap]);
if isempty(idx)
    projData.stats.dynamicity=NaN;
else

    idxCell=mat2cell(idx,ones(length(idx),1),1);
    % sub track indices of the full tracks
    subIdx=cellfun(@(x) find(dataMatCrpSecMic(:,1)==x),idxCell,'uniformoutput',0);
    % full track lifetimes and displacements
    ltf=cell2mat(cellfun(@(x) sum(dataMatCrpSecMic(x,6)),subIdx,'uniformoutput',0));
    disp=cell2mat(cellfun(@(x) sum(abs(dataMatCrpSecMic(x,7))),subIdx,'uniformoutput',0));
    % collective displacement of all gap-containing MTs over their collective lifetime
    projData.stats.dynamicity=sum(disp)/(sum(ltf)/60);
end
%% Comet Latency Parameters

projData.stats.avgComLatSec = projData.stats.fgap_length_mean/projData.stats.growth_speed_mean*60;

%% Growth, fgap, bgap events density
try %#ok<TRYNC>
    area=projData.roiArea*(projData.pixSizeNm/1e3)^2; % area in microns
    time= projData.nFrames*projData.secPerFrame; % time in seconds
    projData.stats.growth_density=projData.stats.nGrowths/area/time;
    projData.stats.fgap_density=projData.stats.nFgaps/area/time;
    projData.stats.bgap_density=projData.stats.nBgaps/area/time;
end
end

