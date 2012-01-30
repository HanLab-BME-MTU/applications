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

% Put a Copy of the Data Mat Used to Calcualte All Stats in the projData

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
% already calculated just move over into stats
% if from subRoi calculate in subRoiExtractTracks (more useful to do that 
% there simply because of data struct paritioning)
if (isfield(projData,'medNNdistWithinFramePix') == 1 && fromPoolGroupData ~=1 ... 
    && subRoiAnalysis ~= 1); 
    projData.stats.medNNdistWithinFrameMic = projData.medNNdistWithinFramePix*projData.pixSizeNm/1000;
else 
end

% field not relevent for pooled data
%if fromPoolGroupData == 1
 %   if isfield(projData.stats,'medNNdistWithinFrameMic') ==1 ; 
 % projData.stats =  rmfield(projData.stats,'medNNdistWithinFrameMic'); 
  %  end 
%end 


%% PARAMETERS RELATED TO GROWTH

if isempty(gIdx)
    projData.stats.nGrowths=0;
    % median/mean for growth speeds (microns/minute)
    projData.stats.growth_speed_median = NaN;
    projData.stats.growth_speed_mean = NaN;
    projData.stats.growth_speed_std = NaN;
   
    
    % median/mean for growth lifetime (sec)
    projData.stats.growth_lifetime_median = NaN;
    projData.stats.growth_lifetime_mean = NaN;
    projData.stats.growth_lifetime_std = NaN;
  
    
    % median/mean for growth displacement (microns)
    projData.stats.growth_length_median = NaN;
    projData.stats.growth_length_mean = NaN;
    projData.stats.growth_length_std = NaN;
   
else
    projData.stats.nGrowths=length(gIdx);
    
    % median/mean for growth speeds (microns/minute)
    projData.stats.growth_speed_median = median(gs);
    projData.stats.growth_speed_mean = mean(gs) ;
    projData.stats.growth_speed_std  = std(gs);
  
    
    % median/mean for growth lifetime (sec)  
    projData.stats.growth_lifetime_median = median(gl);
    projData.stats.growth_lifetime_mean = mean(gl);
    projData.stats.growth_lifetime_std = std(gl);
    
    % median/mean for growth displacement (microns)
    projData.stats.growth_length_median = median(gd);
    projData.stats.growth_length_mean = mean(gd);
    projData.stats.growth_length_std = std(gd);
    
    
   
end

% this is just an old field in an older version: simply remove it so it
% will not give an error
if isfield(projData.stats,'percentFgapsReclass')==1; 
projData.stats = rmfield(projData.stats,'percentFgapsReclass'); 
else 
end 
%% PARAMETERS RELATED TO FGAPS 

if isempty(fIdx)
    projData.stats.nFgaps=0;
   
    % median/mean for fgap speeds (microns/minute)
    projData.stats.fgap_speed_median = NaN;
    projData.stats.fgap_speed_mean = NaN;
    projData.stats.fgap_speed_std = NaN;
    
    
    % median/mean for fgap lifetime (sec)
    projData.stats.fgap_lifetime_median = NaN;
    projData.stats.fgap_lifetime_mean = NaN;
    projData.stats.fgap_lifetime_std = NaN;
   
    
    % median/mean for fgap displacement (microns)
    projData.stats.fgap_length_median = NaN;
    projData.stats.fgap_length_mean = NaN ;
    projData.stats.fgap_length_std = NaN;
    
    
    % fgap frequency is the inverse of the average growth time (sec) or
    % displacement (microns) prior to fgap
    projData.stats.fgap_freq_time_mean=NaN;
    projData.stats.fgap_freq_length_mean=NaN;
    projData.stats.fgap_freq_time_SE = NaN; 
    projData.stats.fgap_freq_length_SE = NaN;
    
    % 
    projData.stats.GrowthSpeedBeforeFgap_MicPerMin_mean = NaN;
    projData.stats.GrowthSpeedBeforeFgap_MicPerMin_SE =  NaN;
    projData.stats.GrowthLifetimeBeforeFgap_Sec_mean= NaN;
    projData.stats.GrowthLifetimeBeforeFgap_Sec_SE = NaN;
    projData.stats.GrowthLengthBeforeFgap_Mic_mean = NaN;
    projData.stats.GrowthLengthBeforeFgap_Mic_SE = NaN;

else
    projData.stats.nFgaps=length(fIdx);
    
    % median/mean for fgap speeds (microns/minute)
    projData.stats.fgap_speed_median = median(fs);
    projData.stats.fgap_speed_mean = mean(fs);
    projData.stats.fgap_speed_std = std(fs);
  
    
    % median/mean for fgap lifetime (sec)
    projData.stats.fgap_lifetime_median = median(fl);
    projData.stats.fgap_lifetime_mean = mean(fl);
    projData.stats.fgap_lifetime_std = std(fl);
    
    
    % median/mean for fgap displacement (microns)
    projData.stats.fgap_length_median = median(fd);
    projData.stats.fgap_length_mean = mean(fd);
    projData.stats.fgap_length_std = std(fd);  
    
    %Calculate Average Time and Displacement Parameters for Growth
    %Preceding Fgap 
    
    beforeFgapIdx = find(dataMatCrpSecMic(:,9) == 2); 
    
    velocityBeforeFgap = dataMatCrpSecMic(beforeFgapIdx,4);
    lifetimeBeforeFgap = dataMatCrpSecMic(beforeFgapIdx,6);
    lengthBeforeFgap =  dataMatCrpSecMic(beforeFgapIdx,7);
    
    projData.stats.GrowthSpeedBeforeFgap_MicPerMin_mean = mean(velocityBeforeFgap);
    projData.stats.GrowthSpeedBeforeFgap_MicPerMin_SE =  std(velocityBeforeFgap)/sqrt(length(velocityBeforeFgap));
    projData.stats.GrowthLifetimeBeforeFgap_Sec_mean= mean(lifetimeBeforeFgap);
    projData.stats.GrowthLifetimeBeforeFgap_Sec_SE = std(lifetimeBeforeFgap)/sqrt(length(lifetimeBeforeFgap));
    projData.stats.GrowthLengthBeforeFgap_Mic_mean = mean(lengthBeforeFgap);
    projData.stats.GrowthLengthBeforeFgap_Mic_SE = std(lengthBeforeFgap)/sqrt(length(lengthBeforeFgap));
  
    
    %Get the Individual Frequencies Corresponding to the Growth Before EACH
    % pause, the average and the std of these values
    freq=1./dataMatCrpSecMic(beforeFgapIdx,6);
    projData.stats.fgap_freq_time_mean = mean(freq); 
    projData.stats.fgap_freq_time_SE = std(freq)/sqrt(length(freq));
    
    freq=1./dataMatCrpSecMic(beforeFgapIdx,7);
    projData.stats.fgap_freq_length_mean = mean(freq);  
    projData.stats.fgap_freq_length_SE = std(freq)/sqrt(length(freq));
         % all track indices where there is either a forward or backward gap

    
end %isempty

%% PARAMETERS RELATED TO BGAPS

if isempty(bIdx)
    projData.stats.nBgaps=0;
  
    % median/mean for bgap speeds (microns/minute)
    projData.stats.bgap_speed_median = NaN;
    projData.stats.bgap_speed_mean = NaN;
    projData.stats.bgap_speed_std = NaN;
   
    
    % median/mean for bgap lifetime (sec)
    projData.stats.bgap_lifetime_median = NaN;
    projData.stats.bgap_lifetime_mean= NaN;
    projData.stats.bgap_lifetime_std = NaN;

    
    % median/mean for bgap displacement (microns)
    projData.stats.bgap_length_median = NaN;
    projData.stats.bgap_length_mean = NaN;
    projData.stats.bgap_length_std = NaN;
  
    
    % bgap frequency is the inverse of the average growth time (sec) or
    % displacement (microns) prior to bgap
    projData.stats.bgap_freq_time_mean=NaN;
    projData.stats.bgap_freq_length_mean=NaN;
    projData.stats.bgap_freq_time_SE = NaN;
    projData.stats.bgap_freq_length_SE = NaN;
    
    projData.stats.GrowthSpeedBeforeBgap_MicPerMin_mean = NaN;
    projData.stats.GrowthSpeedBeforeBgap_MicPerMin_SE = NaN;
    projData.stats.GrowthLifetimeBeforeBgap_Sec_mean = NaN ;
    projData.stats.GrowthLifetimeBeforeBgap_Sec_SE = NaN;
    projData.stats.GrowthLengthBeforeBgap_Mic_mean = NaN ;
    projData.stats.GrowthLengthBeforeBgap_Mic_SE = NaN; 
    

else
    projData.stats.nBgaps=length(bIdx);

    % median/mean for bgap speeds (microns/minute)
    projData.stats.bgap_speed_median = median(bs);
    projData.stats.bgap_speed_mean = mean(bs);
    projData.stats.bgap_speed_std = std(bs);
    
    % median/mean for bgap lifetime (sec)
    projData.stats.bgap_lifetime_median = median(bl);
    projData.stats.bgap_lifetime_mean = mean(bl);
    projData.stats.bgap_lifetime_std = std(bl);
 
    
    % median/mean for bgap displacement (microns)
    projData.stats.bgap_length_median = median(bd);
    projData.stats.bgap_length_mean = mean(bd);
    projData.stats.bgap_length_std = std(bd);
 
    
    
    %%%%% STAT VALUES PRIOR TO BGAP %%%%
    
    beforeBgapIdx= find(dataMatCrpSecMic(:,9) == 3); 
    
    %Calculate Average Time and Displacement Parameters for Growth
    %Preceding Bgap
    
    velocityBeforeBgap = dataMatCrpSecMic(beforeBgapIdx,4);
    lifetimeBeforeBgap = dataMatCrpSecMic(beforeBgapIdx,6);
    lengthBeforeBgap =  dataMatCrpSecMic(beforeBgapIdx,7);
    
    projData.stats.GrowthSpeedBeforeBgap_MicPerMin_mean = mean(velocityBeforeBgap);
    projData.stats.GrowthSpeedBeforeBgap_MicPerMin_SE = std(velocityBeforeBgap)/sqrt(length(velocityBeforeBgap));
    projData.stats.GrowthLifetimeBeforeBgap_Sec_mean = mean(lifetimeBeforeBgap) ;
    projData.stats.GrowthLifetimeBeforeBgap_Sec_SE = std(lifetimeBeforeBgap)/sqrt(length(lifetimeBeforeBgap));
    projData.stats.GrowthLengthBeforeBgap_Mic_mean = mean(lengthBeforeBgap) ;
    projData.stats.GrowthLengthBeforeBgap_Mic_SE = std(lengthBeforeBgap)/sqrt(length(lengthBeforeBgap));
   
    
    %Convert these Average Values to frequencies
    %projData.stats.bgap_freq_time=1/mean(lifetimeBeforeBgap);
    %projData.stats.bgap_freq_length=1/mean(lengthBeforeBgap);
    
    % Find Freqency Values for each subtrack and find avg and std
      freq=1./dataMatCrpSecMic(beforeBgapIdx,6);
      projData.stats.bgap_freq_time_mean = mean(freq); 
      projData.stats.bgap_freq_time_SE = std(freq)/sqrt(length(freq));
    
    freq =1./dataMatCrpSecMic(beforeBgapIdx,7);
    projData.stats.bgap_freq_length_mean = mean(freq);
    projData.stats.bgap_freq_length_SE= std(freq)/sqrt(length(freq));
  
end % isempty

%% PARAMETERS OF GROWTH SUBTRACKS PRECEDING TERMINAL EVENT 
% 


  beforeTermIdx  = find(dataMatCrpSecMic(:,9) == 1); 

  if isempty(beforeTermIdx)
      projData.stats.GrowthSpeedBeforeTermEvent_MicPerMin_mean = NaN; 
      
      projData.stats.GrowthSpeedBeforeTermEvent_MicPerMin_mean = NaN;
      projData.stats.GrowthSpeedBeforeTermEvent_MicPerMin_SE = NaN ;
      projData.stats.GrowthLifetimeBeforeTermEvent_Sec_mean = NaN;
      projData.stats.GrowthLifetimeBeforeTermEvent_Sec_SE = NaN;
      projData.stats.GrowthLengthBeforeTermEvent_Mic_mean = NaN; 
      projData.stats.GrowthLengthBeforeTermEvent_Mic_SE = NaN;
      projData.stats.term_freq_time_mean =NaN;
      projData.stats.term_freq_time_SE = NaN;
       projData.stats.term_freq_length_mean = NaN;
       projData.stats.term_freq_length_SE = NaN; 
      
      
  else 
      
      
    % collect the parameter under question for all subtracks
    velocityBeforeTermEvent = dataMatCrpSecMic(beforeTermIdx,4);
    lifetimeBeforeTermEvent = dataMatCrpSecMic(beforeTermIdx,6);
    lengthBeforeTermEvent = dataMatCrpSecMic(beforeTermIdx,7); 
    
    %%% AVERAGE VELOCITY, TIME, AND DISPLACEMENT PRECEDING TERMINAL EVENT%%% 
    
    projData.stats.GrowthSpeedBeforeTermEvent_MicPerMin_mean = mean(velocityBeforeTermEvent);
    projData.stats.GrowthSpeedBeforeTermEvent_MicPerMin_SE = std(velocityBeforeTermEvent)/sqrt(length(velocityBeforeTermEvent)) ;
    projData.stats.GrowthLifetimeBeforeTermEvent_Sec_mean = mean(lifetimeBeforeTermEvent);
    projData.stats.GrowthLifetimeBeforeTermEvent_Sec_SE = std(lifetimeBeforeTermEvent)/sqrt(length(lifetimeBeforeTermEvent));
    projData.stats.GrowthLengthBeforeTermEvent_Mic_mean = mean(lengthBeforeTermEvent); 
    projData.stats.GrowthLengthBeforeTermEvent_Mic_SE = std(lengthBeforeTermEvent)/sqrt(length(lengthBeforeTermEvent));
    
    %%% FREQUENCIES (INVERSES OF AVG TIME AND DISPLACEMENT %%%%
  
    %projData.stats.term_freq_time=1/mean(termGrowthOnly(:,6));
    %projData.stats.term_freq_length=1/mean(termGrowthOnly(:,7));
    
    %%% INDIVIDUAL FREQUENCIES %%%
    
    freq=1./dataMatCrpSecMic(beforeTermIdx,6);
    projData.stats.term_freq_time_mean =mean(freq);
    projData.stats.term_freq_time_SE = std(freq)/sqrt(length(freq));
    
    freq=1./dataMatCrpSecMic(beforeTermIdx,7);
    projData.stats.term_freq_length_mean = mean(freq);
    projData.stats.term_freq_length_SE = std(freq)/sqrt(length(freq));
  end 
   
    %%% RATIOS OF GROWTH VELOCITY, GROWTH LIFETIME, OR GROWTH LENGTH JUST 
    % BEFRORE AN FGAP/BGAP/TERM EVENT

    % Velocity Ratios 
    projData.stats.ratio_preFgapVel2preTermVel = projData.stats.GrowthSpeedBeforeFgap_MicPerMin_mean/projData.stats.GrowthSpeedBeforeTermEvent_MicPerMin_mean;
    projData.stats.ratio_preBgapVel2preTermVel = projData.stats.GrowthSpeedBeforeBgap_MicPerMin_mean/projData.stats.GrowthSpeedBeforeTermEvent_MicPerMin_mean;
    projData.stats.ratio_preFgapVel2preBgapVel= projData.stats.GrowthSpeedBeforeFgap_MicPerMin_mean/projData.stats.GrowthSpeedBeforeBgap_MicPerMin_mean;
    
    % Lifetime Ratios
    projData.stats.ratio_preFgapLife2preTermLife  = projData.stats.GrowthLifetimeBeforeFgap_Sec_mean/projData.stats.GrowthLifetimeBeforeTermEvent_Sec_mean;
    projData.stats.ratio_preBgapLife2preTermLife= projData.stats.GrowthLifetimeBeforeBgap_Sec_mean/projData.stats.GrowthLifetimeBeforeTermEvent_Sec_mean;
    projData.stats.ratio_preFgapLife2preBgapLife  = projData.stats.GrowthLifetimeBeforeFgap_Sec_mean/projData.stats.GrowthLifetimeBeforeBgap_Sec_mean;
    
    % Displacment (length) ratios
    projData.stats.ratio_preFgapDisp2preTermDisp = projData.stats.GrowthLengthBeforeFgap_Mic_mean/projData.stats.GrowthLengthBeforeTermEvent_Mic_mean;
    projData.stats.ratio_preFgapDisp2preTermDisp = projData.stats.GrowthLengthBeforeBgap_Mic_mean/projData.stats.GrowthLengthBeforeTermEvent_Mic_mean;
    projData.stats.ratio_preFgapDisp2preBgapDisp = projData.stats.GrowthLengthBeforeFgap_Mic_mean/projData.stats.GrowthLengthBeforeBgap_Mic_mean;
    
 
     
%% MISC PARAMETERS

% percent of time spent in growth, fgap, and bgap: 

totalTime=sum(gl)+sum(fl)+sum(bl);
if totalTime==0
    projData.stats.percentTimeGrowth=NaN;
    projData.stats.percentTimeFgap  =NaN;
    projData.stats.percentTimeBgap  =NaN;
else
    projData.stats.percentTimeGrowth= 100*(sum(gl)/totalTime);
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

beforeUndefinedIdx = find(dataMatCrpSecMic(:,9) == 4);

if isempty(gIdx)
    projData.stats.percentGrowthLinkedForward  = NaN;
    projData.stats.percentGrowthLinkedBackward = NaN;
    projData.stats.percentGrowthTerminal       = NaN;
else 
    if isempty(bIdx)
        projData.stats.percentGrowthLinkedBackward = NaN;
    else 
    projData.stats.percentGrowthLinkedBackward = 100*length(beforeBgapIdx)/length(gIdx);
    end
    
    if isempty(fIdx) 
        projData.stats.percentGrowthLinkedForward = NaN; 
    else 
    projData.stats.percentGrowthLinkedForward = 100*length(beforeFgapIdx)/length(gIdx); 
    end 
    
    if isempty(beforeTermIdx)
        projData.percentGrowthTerminal = NaN; 
    else 
    projData.stats.percentGrowthTerminal = 100*length(beforeTermIdx)/length(gIdx); 
    end 
    
    if isempty(beforeUndefinedIdx) 
    projData.stats.percentGrowthLinkedUndefinedGap = NaN; 
    else 
        projData.stats.percentGrowthLinkedUndefinedGap = 100*length(beforeUndefinedIdx)/length(gIdx); 
    end  
        
end 


%% More Microtubule Parameters:
% These parameters require in-tact compound microtubule tracks for correct 
% calculations and therefore may not always be correctly calculated 
% when the tracks are partioned based on spatial criteria. 
% The subRoi analysis is currently under construction and as this
% evolves for our own laboratory needs we will continue to update. 
% For now we remove these parameters from the subRoi statistics. 
   
   if (subRoiAnalysis == 1 || fromPoolGroupData == 1) % if not calling this from original analysis 
       % ie either pooling group data or subRoianalysis 
       
        % do not perform compound track analysis: Remove these fields as they
       % represent the original analysis and this can be confusing. 
      if isfield(projData,'singleDataMat') == 1 
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
      projData.stats = rmfield(projData.stats,'ratio_TotalTracks2NumGrowthSubtracks'); 
      projData.stats = rmfield(projData.stats,'NumTracksWithfGap2TotalTracks_Per'); 
      projData.stats = rmfield(projData.stats,'NumTracksWithbGap2TotalTracks_Per'); 
      projData.stats = rmfield(projData.stats,'NumTracksWithGap2TotalTracks_Per'); 
      projData.stats = rmfield(projData.stats, 'avgIndivPercentTimeFgap'); 
      projData.stats = rmfield(projData.stats,'avgIndivPercentTimeBgap'); 
      projData.stats  = rmfield(projData.stats, 'dynamicity'); 
      
      
      end 
       
   else % from whole cell analysis so perform this 
   %% Compound Vs. Single Track analysis
      numGrowthSubTracksAll = length(find(dataMatCrpSecMic(:,5) == 1));
   
      numTracksTotal = length(unique(dataMatCrpSecMic(gIdx,1)));
      projData.stats.ratio_TotalTracks2NumGrowthSubtracks =  numTracksTotal/numGrowthSubTracksAll;
     
      tracksWithFgap=unique(dataMatCrpSecMic(fIdx,1));
      tracksWithBgap=unique(dataMatCrpSecMic(bIdx,1));
      tracksWithGap = unique(dataMatCrpSecMic([fIdx ;bIdx],1));
      
      projData.stats.NumTracksWithfGap2TotalTracks_Per = length(tracksWithFgap)/numTracksTotal*100;
      projData.stats.NumTracksWithbGap2TotalTracks_Per = length(tracksWithBgap)/numTracksTotal*100;
      projData.stats.NumTracksWithGap2TotalTracks_Per =  length(tracksWithGap)/numTracksTotal*100;     
      
   %load compound track data before remove any data
   compDataMat = projData.compDataMat; 
   singleDataMat = projData.singleDataMat; 
   
   
   
   
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

  %% Percent Time Individual Tracks    
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

   end % if subRoiAnalysis
       
%% Comet Latency Parameters

projData.stats.avgComLatSec = projData.stats.fgap_length_mean/projData.stats.growth_speed_mean*60;

%% Growth, fgap, bgap events density
 if isfield(projData,'roiArea') == 1 || isfield(projData,'subRoiArea')
    if subRoiAnalysis == 1
    area = projData.subRoiAreaSqMic;
    else 
    area=projData.roiArea*(projData.pixSizeNm/1e3)^2; % area in microns
    end 
    time= projData.nFrames*projData.secPerFrame; % time in seconds
    projData.stats.growth_density=projData.stats.nGrowths/area/time;
    numNucEvents = sum(dataMatCrpSecMic(:,8)); 
    projData.stats.numNucleationEvents = numNucEvents; 
    projData.stats.nucleationDensity = numNucEvents/area/time; 
    projData.stats.fgap_density=projData.stats.nFgaps/area/time;
    projData.stats.bgap_density=projData.stats.nBgaps/area/time;
 end 

end

