function [projData]=plusTipPostTracking(runInfo,secPerFrame,pixSizeNm,timeRange)
% plusTipPostTracking extracts various statistics from EB3 tracks
%
%SYNOPSIS [projData]=plusTipPostTracking(runInfo,secPerFrame,pixSizeNm)
%
%INPUT  runInfo           : structure containing fields .anDir, which gives
%                           the full path to the roi_x directory
%                           and .imDir, which gives the full path to the
%                           folder containing the images for overlay.
%                           if given as [], program will query user for
%                           roi_x directory.
%       secPerFrame       : frame rate (seconds per frame)
%       pixSizeNm         : real-space pixel size in nanometers
%       timeRange         : frame range over which to run post-processing
%
%OUTPUT projData          : structure with the following fields, saved in
%                           a folder /roi_x/meta
%  .imDir
%       image directory path
%  .anDir
%       roi directory path
%  .trackingParameters
%       parameters used in the tracking step
%  .secPerFrame
%       frame rate (seconds per frame)
%  .pixSizeNm
%       real-space pixel size in nanometers
%  .numTracks
%       number of tracks
%  .numFrames
%       number of time points (frames)
%  .xCoord/yCoord
%       pixel coordinates of features in tracks
%  .featArea
%       area in pixels of feature from detection (for growth only)
%  .featInt
%       max intensity of feature from detection (for growth only)
%  .frame2frameDispPix
%       vector containing all frame-to-frame displacements from the initial
%       growth subtrack portions of all tracks
%  .pair2pairDiffPix
%       vector containing the difference in displacement between
%       frame pairs from the initial growth trajectories (before gap
%       closing) - provides a roughly normal distribution, the width of
%       which indicates how much intra-track speed variation there is
%  .pair2pairDiffMicPerMinStd
%       standard deviation of pair2pairDiffPix, converted to microns/min
%  .medNNdistWithinFramePix
%       median nearest-neighbor distance between all detected features in a
%       frame from movieInfo (loops through all frames to get composite
%       from full movie, not just first frame)
%  .meanDisp2medianNNDistRatio
%       mean displacement / median NN - 0.3 is upper bound tested in
%       Jaqaman Nature Methods 2008; this statistic is a measure of how
%       challenging the tracking is.  A higher number means it is more
%       challenging.
%  .frame2frameVel_micPerMin
%       intantaneous velocity based on interpolated positions during gaps,
%       converted to microns/minute
%  .segGapAvgVel_micPerMin
%       average velocity based on interpolated positions during gaps,
%       converted to microns/minute. takes into account reclassified fgaps
%  .percentFgapsReclass
%       percentage of fgaps that get reclassified as continuation of growth
%       because their speeds are >=90% the speed at the end of the growth
%       phase prior.  these reclassified fgaps get removed and the growth
%       phases before and after become joined into one
%  .tracksWithFgap
%       track numbers corresponding to tracks that contain forward gaps
%       (pause or out of focus events)
%  .tracksWithBgap
%       track numbers corresponding to tracks that contain backward gaps
%       (shrinkage events)
%  .stats
%       structure containing the following fields:
%
%           .nGrowths
%               total number of growth trajectories beginning after the
%               first frame and ending before the last frame
%           .growth_speed_median
%               median of all growth trajectory speeds (microns/min)
%           .growth_speed_mean_SE 
%               [mean SE] of all growth trajectory speeds (microns/min),
%               where SE is std/sqrt(N)
%           .growth_lifetime_median
%               median of all growth trajectory lifetimes (min)
%           .growth_lifetime_mean_SE
%               [mean SE] of all growth trajectory lifetimes (min),
%               where SE is std/sqrt(N)
%           .growth_length_median
%               median of all growth trajectory lengths (microns)
%           .growth_length_mean_SE
%               [mean SE] of all growth trajectory lengths (microns),
%               where SE is std/sqrt(N)
%
%           .nFgaps
%               total number of forward gap trajectories
%           .fgap_speed_median
%               median of all forward gap trajectory speeds (microns/min)
%           .fgap_speed_mean_SE 
%               [mean SE] of all forward gap trajectory speeds (microns/min),
%               where SE is std/sqrt(N)
%           .fgap_lifetime_median
%               median of all forward gap trajectory lifetimes (min)
%           .fgap_lifetime_mean_SE
%               [mean SE] of all forward gap trajectory lifetimes (min),
%               where SE is std/sqrt(N)
%           .fgap_length_median
%               median of all forward gap trajectory lengths (microns)
%           .fgap_length_mean_SE
%               [mean SE] of all forward gap trajectory lengths (microns),
%               where SE is std/sqrt(N)
%           .fgap_freq_time
%               [mean(1/D), SE] where D is the growth lifetime (min) just
%               prior to a forward gap, excluding growth trajectories
%               lasting only 3 frames and those beginning in the first
%               frame, and SE is std/sqrt(N)
%           .fgap_freq_length     
%               [mean(1/L), SE] where L is the growth length (microns) just
%               prior to a forward gap, excluding growth trajectories
%               lasting only 3 frames and those beginning in the first
%               frame, and SE is std/sqrt(N)
%
%           .nBgaps
%               total number of backward gap trajectories
%           .bgap_speed_median
%               median of all backward gap trajectory speeds (microns/min)
%           .bgap_speed_mean_SE 
%               [mean SE] of all backward gap trajectory speeds (microns/min),
%               where SE is std/sqrt(N)
%           .bgap_lifetime_median
%               median of all backward gap trajectory lifetimes (min)
%           .bgap_lifetime_mean_SE
%               [mean SE] of all backward gap trajectory lifetimes (min),
%               where SE is std/sqrt(N)
%           .bgap_length_median
%               median of all backward gap trajectory lengths (microns)
%           .bgap_length_mean_SE
%               [mean SE] of all backward gap trajectory lengths (microns),
%               where SE is std/sqrt(N)
%           .bgap_freq_time
%               [mean(1/D), SE] where D is the growth lifetime (min) just
%               prior to a backward gap, excluding growth trajectories
%               lasting only 3 frames and those beginning in the first
%               frame, and SE is std/sqrt(N)
%           .bgap_freq_length     
%               [mean(1/L), SE] where L is the growth length (microns) just
%               prior to a backward gap, excluding growth trajectories
%               lasting only 3 frames and those beginning in the first
%               frame, and SE is std/sqrt(N)
%
%           .projData.percentTimeGrowth
%               time all tracks spend in growth over time all tracks spend
%               in growth, fgap, or bgap
%           .projData.percentTimeFgap
%               time all tracks spend in fgap over time all tracks spend
%               in f, fgap, or bgap
%           .projData.percentTimeBgap
%               time all tracks spend in bgap over time all tracks spend
%               in growth, fgap, or bgap
%
%           .percentGapsForward 
%               percent of all gaps which are fgaps
%           .percentGapsBackward
%               percent of all gaps which are bgaps
%           .percentGrowthLinkedForward
%                percent of growth trajectories beginning after the first
%                frame and ending before the last frame that get linked to
%                fgaps
%           .percentGrowthLinkedBackward
%                percent of growth trajectories beginning after the first
%                frame and ending before the last frame that get linked to
%                bgaps
%           .percentGrowthTerminal
%                percent of growth trajectories beginning after the first
%                frame and ending before the last frame that do not get
%                linked to either fgaps or bgaps
%           .dynamicity
%               collective displacement of all gap-containing MTs over
%               their collective lifetime
%
%  .nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix
%       matrix containing the profile of all the tracks, where the columns
%       represent:
%           1. track number
%           2. start frame
%           3. end frame
%           4. speed (microns/min), negative for bgaps
%           5. subtrack type
%              (1=growth, 2=forward gap, 3=backward gap, 4=unclassifed gap)
%           6. lifetime (frames)
%           7. total displacement (pixels), negative for bgaps



% get runInfo in correct format
if nargin<1 || isempty(runInfo)
    homeDir=pwd;
    runInfo.anDir=uigetdir(pwd,'Please select analysis directory');
    cd([runInfo.anDir filesep '..'])
    runInfo.imDir=uigetdir(pwd,'Please select image directory');
    cd(homeDir)
else
    % adjust for OS
    if ~isfield(runInfo,'imDir') || ~isfield(runInfo,'anDir')
        error('--plusTipPostTracking: first argument should be a structure with fields imDir and anDir');
    else
        [runInfo.anDir] = formatPath(runInfo.anDir);
        homeDir=pwd;
        cd(runInfo.anDir)
        [runInfo.imDir] = formatPath(runInfo.imDir);
        cd(homeDir)
    end
end

if nargin<2 || isempty(secPerFrame)
    error('--plusTipPostTracking: frame rate missing')
end
if nargin<3 || isempty(pixSizeNm)
    error('--plusTipPostTracking: pixel size missing')
end
if nargin<4 || isempty(timeRange)
    timeRange=[];
end


% load movieInfo (detection result)
featDir  = [runInfo.anDir filesep 'feat'];
if ~isdir(featDir)
    error('--plusTipPostTracking: feat directory missing')
else
    if exist([featDir filesep 'movieInfo.mat'],'file')
        load([featDir filesep 'movieInfo.mat'])
    else
        error('--plusTipPostTracking: movieInfo missing...')
    end
end

% load tracksFinal (tracking result)
trackDir = [runInfo.anDir filesep 'track'];
if ~isdir(trackDir)
    error('--plusTipPostTracking: track directory missing')
else
    [listOfFiles]=searchFiles('.mat',[],trackDir,0);
    if ~isempty(listOfFiles)
        load([listOfFiles{1,2} filesep listOfFiles{1,1}])
        if ~exist('tracksFinal','var')
            error('--plusTipPostTracking: tracksFinal missing...');
        end
    else
        error('--plusTipPostTracking: tracksFinal missing...');
    end
end

% make a new "meta" directory if it does not exist
runInfo.metaDir = [runInfo.anDir filesep 'meta'];
if ~isdir(runInfo.metaDir)
    mkdir(runInfo.metaDir)
end

% get interpolated positions for gaps and calculate velocities
[trackedFeatureInfo,trackedFeatureInfoInterp,trackInfo,trackVelocities,timeRange]=...
    getVelocitiesFromMat(tracksFinal,movieInfo,3,timeRange);

%get number of tracks and number of time points
[numTracks,numTimePoints] = size(trackedFeatureInfo);
numTimePoints=numTimePoints/8;


% without interpolation yet
x = trackedFeatureInfo(:,1:8:end);
y = trackedFeatureInfo(:,2:8:end);

% initialize matrices for feature indices, area, and intensity
movieInfoIdx=nan(numTracks,numTimePoints);
featArea=nan(numTracks,numTimePoints);
featInt =nan(numTracks,numTimePoints);
for iFrame=1:numTimePoints
    % these are the track numbers which exist in iFrame
    existCoordIdx=find(~isnan(x(:,iFrame)));
    % these are the corresponding xy-coordinates
    xi=x(:,iFrame); xi(isnan(xi))=[];
    yi=y(:,iFrame); yi(isnan(yi))=[];

    if ~isempty(xi)
        % distance matrix reveals where features coincide with those recorded
        % in movieInfo
        D=createDistanceMatrix([xi,yi],[movieInfo(iFrame,1).xCoord(:,1),movieInfo(iFrame,1).yCoord(:,1)]);
        [r,c]=find(D==0); % r=track index, c=frame

        [newR,idx]=sort(r); % re-order based on track
        featIdx=c(idx); % movieInfo feature index, sorted to correspond to track indices

        % fill in movieInfoIdx with indices from features stored in movieInfo
        movieInfoIdx(existCoordIdx,iFrame)=featIdx;
        % fill in feature area (pixels) at corresponding features
        featArea(existCoordIdx,iFrame)=movieInfo(iFrame,1).amp(featIdx,1);
        % fill in feature max intensity at corresponding features
        featInt (existCoordIdx,iFrame)=movieInfo(iFrame,1).int(featIdx,1);
    end
end


% save misc info for output
projData.imDir = runInfo.imDir;
projData.anDir = runInfo.anDir;

projData.trackingParameters.timeWindow=gapCloseParam.timeWindow;
projData.trackingParameters.minTrackLen=gapCloseParam.minTrackLen;
projData.trackingParameters.minSearchRadius=costMatrices(1,1).parameters.minSearchRadius;
projData.trackingParameters.maxSearchRadius=costMatrices(1,1).parameters.maxSearchRadius;
projData.trackingParameters.maxFAngle=costMatrices(1,2).parameters.maxFAngle;
projData.trackingParameters.maxBAngle=costMatrices(1,2).parameters.maxBAngle;
projData.trackingParameters.backVelMultFactor=costMatrices(1,2).parameters.backVelMultFactor;
projData.trackingParameters.fluctRad=costMatrices(1,2).parameters.fluctRad;

projData.secPerFrame = secPerFrame;
projData.pixSizeNm = pixSizeNm;

projData.numTracks = numTracks;
projData.numFrames = numTimePoints;

% figure out which frames were used in detection
m=struct2cell(movieInfo); m=m(1,:); 
detExists=find(cellfun(@(x) ~isempty(x),m)); 
sF=min(detExists); eF=max(detExists);

% frame ranges for each step
projData.detectionFrameRange=[sF eF];
projData.trackingFrameRange=[costMatrices(1).parameters.startFrame costMatrices(1).parameters.endFrame];
projData.postTrackFrameRange = timeRange;

% coordinate/area/intensity info from detected features
projData.xCoord = trackedFeatureInfoInterp(:,1:8:end);
projData.yCoord = trackedFeatureInfoInterp(:,2:8:end);
projData.featArea = featArea;
projData.featInt = featInt;


% get frame-to-frame displacement for growth only (not forward/backward gaps)
frame2frameDispPix=sqrt(diff(x,1,2).^2+diff(y,1,2).^2);
% get rid of NaNs and linearize the vector
projData.frame2frameDispPix=frame2frameDispPix(~isnan(frame2frameDispPix(:)));

% get change in velocity between frame *pairs* for segments only
pair2pairDiffPix=diff(frame2frameDispPix,1,2);
% get rid of NaNs and linearize the vector
projData.pair2pairDiffPix=pair2pairDiffPix(~isnan(pair2pairDiffPix(:)));
% std (microns/min) of delta growthSpeed btw frames
projData.pair2pairDiffMicPerMinStd=std(pixPerFrame2umPerMin(projData.pair2pairDiffPix,projData.secPerFrame,projData.pixSizeNm));


% get all feature nearest neighbor distances from all frames in one vector
NNdist=nan(length(vertcat(movieInfo.xCoord)),1);
count=1;
for iFrame=5:length(movieInfo)

    xCoord = movieInfo(iFrame).xCoord;
    yCoord = movieInfo(iFrame).yCoord;
    
    if ~isempty(xCoord)
        xCoord=xCoord(:,1);
        yCoord=yCoord(:,1);

        D=createDistanceMatrix([xCoord yCoord],[xCoord yCoord]);
        [sD,idx]=sort(D,2);

        NNdist(count:count+length(xCoord)-1)=sD(:,2);

    end
    count=count+length(xCoord);
end

% median NN dist
projData.medNNdistWithinFramePix=nanmedian(NNdist);

% get mean displacement to median NN distance ratio
projData.meanDisp2medianNNDistRatio = mean(projData.frame2frameDispPix)/projData.medNNdistWithinFramePix;

% convert interpolated velocities to microns per minute
[projData.frame2frameVel_micPerMin]=pixPerFrame2umPerMin(trackVelocities.frame2frame,secPerFrame,pixSizeNm);

% concatenate all segments and gaps into n x 4 matrices, then add info:
% [trackNum startFrame endFrame velocity seg/gapType trackLengthFrames]
segs   = vertcat(trackInfo.seg);
fgaps  = vertcat(trackInfo.fgap);
bgaps  = vertcat(trackInfo.bgap);
ugaps  = vertcat(trackInfo.ugap);

compositeMatrix = [];
if ~isempty(segs)
    segs =  [segs   1*ones(size(segs,1),1)];
    compositeMatrix = [compositeMatrix; segs];
end
if ~isempty(fgaps)
    fgaps = [fgaps  2*ones(size(fgaps,1),1)];
    compositeMatrix = [compositeMatrix; fgaps];
end
if ~isempty(bgaps)
    bgaps = [bgaps  3*ones(size(bgaps,1),1)];
    compositeMatrix = [compositeMatrix; bgaps];
end
if ~isempty(ugaps)
    ugaps = [ugaps  4*ones(size(ugaps,1),1)];
    compositeMatrix = [compositeMatrix; ugaps];
end


% put segs/gaps into one matrix and sort to see track profiles in order
aT=sortrows(compositeMatrix,[1 2]);

% lifetime is subtrack length in frames
lifeTimes=aT(:,3)-aT(:,2);
% get total distance traveled over all subtracks
totalDispPix=aT(:,4).*lifeTimes;
% add lifetime and total distplacement to matrix
aT=[aT lifeTimes totalDispPix];


% look at fgaps - if their speed << growth phase just prior, then it's a
% true pause.  it it is greater than 90% growth speed, reclassify it as
% continuation of growth.
fgapIdx=find(aT(:,5)==2);
beforeFgapIdx=fgapIdx-1;
% get speed of growth phase at last 2-3 time pts (depending on how long
% growth phase lasted) just prior to pause
eF=aT(beforeFgapIdx,3);
beforeFgapSpeed=zeros(length(fgapIdx),1);
for iGap=1:length(fgapIdx)
    idxRange=[max(1,eF(iGap)-3):eF(iGap)-1];
    beforeFgapSpeed(iGap)=nanmean(projData.frame2frameVel_micPerMin(aT(beforeFgapIdx(iGap),1),idxRange));
end

% these are the fgaps to consolidate
growthFgapIdx=fgapIdx(aT(fgapIdx,4)>0.9.*beforeFgapSpeed);
% these are the affected track numbers
tracks2check=unique(aT(growthFgapIdx,1));
aTrows2remove=[];
for i=1:length(tracks2check)
    % rows of aT corresponding to i track
    subIdx=find(aT(:,1)==tracks2check(i));
    % rows of aT corresponding to fgaps in track to consolidate
    fgap2remIdx=intersect(growthFgapIdx,subIdx);
    % rows of aT corresonding to bgaps or real pauses in track
    sepIdx=union(subIdx(aT(subIdx,5)==3),setdiff(intersect(subIdx,fgapIdx),fgap2remIdx))';
    % split the track based on bgaps, so that all fgaps that
    % should be consolidated together can be done at the same time
    sIdx=[subIdx(1); sepIdx(:)];
    eIdx=[sepIdx(:); subIdx(end)];
    % loop through groups of subtracks (split by bgaps)
    for j=1:length(sIdx)
        % pTemp contains the ones to consolidate in this section
        pTemp=intersect(fgap2remIdx,sIdx(j):eIdx(j));
        if ~isempty(pTemp)
            fIdx=min(pTemp)-1; % first row - prior to first fgap
            lIdx=max(pTemp)+1; % last row - after final fgap

            aT(fIdx,3)=aT(lIdx,3); % change end frame
            aT(fIdx,6)=sum(aT(fIdx:lIdx,6)); % sum lifetimes
            aT(fIdx,7)=sum(aT(fIdx:lIdx,7)); % sum total displacements
            aT(fIdx,4)=aT(fIdx,7)/aT(fIdx,6); % find new average velocity

            % keep track of which are the extra rows
            aTrows2remove=[aTrows2remove fIdx+1:lIdx];
        end
    end
end
aT(aTrows2remove,:)=[];

% recalculate segment averages to reflect consolidation
projData.segGapAvgVel_micPerMin=zeros(size(projData.frame2frameVel_micPerMin));
for iSub=1:size(aT,1)
    projData.segGapAvgVel_micPerMin(aT(iSub,1),aT(iSub,2):aT(iSub,3)-1)=aT(iSub,4);
end


% fraction of fgaps that were consolidated into growth, since their speeds
% were more than 90% of the growth speed just prior to the gap
projData.percentFgapsReclass=100*(length(aTrows2remove)/2)/length(fgapIdx);

% get track numbers that contain an fgap or bgap
projData.tracksWithFgap = unique(aT(aT(:,5)==2,1));
projData.tracksWithBgap = unique(aT(aT(:,5)==3,1));

% convert pix/frame to micron/min velocities
[aT(:,4)] = pixPerFrame2umPerMin(aT(:,4),secPerFrame,pixSizeNm);

% index for growths which don't start before or end after the movie
gIdx=find(aT(:,5)==1 & aT(:,2)>1 & aT(:,3)<max(aT(:,3)));
projData.stats.nGrowths=length(gIdx);
% median/mean for growth speeds (microns/minute)
gs=aT(gIdx,4);
projData.stats.growth_speed_median = median(gs);
projData.stats.growth_speed_mean_SE = [mean(gs) std(gs)/sqrt(length(gs))];
% median/mean for growth lifetime (min)
gl=(aT(gIdx,6).*projData.secPerFrame)/60;
projData.stats.growth_lifetime_median = median(gl);
projData.stats.growth_lifetime_mean_SE = [mean(gl) std(gl)/sqrt(length(gl))];
% median/mean for growth displacement (microns)
gd=(abs(aT(gIdx,7)).*projData.pixSizeNm)/1000;
projData.stats.growth_length_median = median(gd);
projData.stats.growth_length_mean_SE = [mean(gd) std(gd)/sqrt(length(gd))];


% index for fgaps
fIdx=find(aT(:,5)==2);
projData.stats.nFgaps=length(fIdx);
% median/mean for fgap speeds (microns/minute)
fs=aT(fIdx,4);
projData.stats.fgap_speed_median = median(fs);
projData.stats.fgap_speed_mean_SE = [mean(fs) std(fs)/sqrt(length(fs))];
% median/mean for fgap lifetime (min)
fl=(aT(fIdx,6).*projData.secPerFrame)/60;
projData.stats.fgap_lifetime_median = median(fl);
projData.stats.fgap_lifetime_mean_SE = [mean(fl) std(fl)/sqrt(length(fl))];
% median/mean for fgap displacement (microns)
fd=(abs(aT(fIdx,7)).*projData.pixSizeNm)/1000;
projData.stats.fgap_length_median = median(fd);
projData.stats.fgap_length_mean_SE = [mean(fd) std(fd)/sqrt(length(fd))];
% frequency of fgap is the average of 1 over the total time (min) spent
% growing prior to fgap
beforeFgapIdx=fIdx-1;
% don't count the ones in the first frame, since we don't know how long
% they are, or those which last for only 3 frames
beforeFgapIdx=beforeFgapIdx(aT(beforeFgapIdx,2)~=1 & aT(beforeFgapIdx,6)>2);
if isempty(beforeFgapIdx)
    projData.stats.fgap_freq_time=NaN;
    projData.stats.fgap_freq_length=NaN;
else
    freq=1./(aT(beforeFgapIdx,6).*(secPerFrame/60));
    projData.stats.fgap_freq_time=[mean(freq) std(freq)/sqrt(length(freq))];
    
    freq=1./(aT(beforeFgapIdx,7).*(pixSizeNm/1000));
    projData.stats.fgap_freq_length=[mean(freq) std(freq)/sqrt(length(freq))];
end

% index for bgaps
bIdx=find(aT(:,5)==3);
projData.stats.nBgaps=length(bIdx);
% median/mean for bgap speeds (microns/minute)
bs=abs(aT(bIdx,4));
projData.stats.bgap_speed_median = median(bs);
projData.stats.bgap_speed_mean_SE = [mean(bs) std(bs)/sqrt(length(bs))];
% median/mean for bgap lifetime (min)
bl=(aT(bIdx,6).*projData.secPerFrame)/60;
projData.stats.bgap_lifetime_median = median(bl);
projData.stats.bgap_lifetime_mean_SE = [mean(bl) std(bl)/sqrt(length(bl))];
% median/mean for bgap displacement (microns)
bd=(abs(aT(bIdx,7)).*projData.pixSizeNm)/1000;
projData.stats.bgap_length_median = median(bd);
projData.stats.bgap_length_mean_SE = [mean(bd) std(bd)/sqrt(length(bd))];
% frequency of bgap is the average of 1 over the total time (min) spent
% growing prior to bgap
beforeBgapIdx=bIdx-1;
% don't count the ones in the first frame, since we don't know how long
% they are, or those which last for only 3 frames
beforeBgapIdx=beforeBgapIdx(aT(beforeBgapIdx,2)~=1 & aT(beforeBgapIdx,6)>2); 
if isempty(beforeBgapIdx)
    projData.stats.bgap_freq_time=NaN;
    projData.stats.bgap_freq_length=NaN;
else
    freq=1./(aT(beforeBgapIdx,6).*(secPerFrame/60));
    projData.stats.bgap_freq_time=[mean(freq) std(freq)/sqrt(length(freq))];
    
    freq=1./(aT(beforeBgapIdx,7).*(pixSizeNm/60));
    projData.stats.bgap_freq_length=[mean(freq) std(freq)/sqrt(length(freq))];
end

% percent of time spent in growth, fgap, and bgap
totalTime=sum(gl)+sum(fl)+sum(bl);
projData.stats.percentTimeGrowth=100*(sum(gl)/totalTime);
projData.stats.percentTimeFgap  =100*(sum(fl)/totalTime);
projData.stats.percentTimeBgap  =100*(sum(bl)/totalTime);

% percent nFgaps/nGaps
projData.stats.percentGapsForward = 100*(projData.stats.nFgaps/(projData.stats.nFgaps+projData.stats.nBgaps));
% percent nBgaps/nGaps
projData.stats.percentGapsBackward= 100*(projData.stats.nBgaps/(projData.stats.nFgaps+projData.stats.nBgaps));


% if next one after growth has index 2, it's an fgap; if 3, it's a bgap; if
% 1 or 4, it is unlinked to another subtrack
f=sum(aT(gIdx+1,5)==2);
b=sum(aT(gIdx+1,5)==3);
u=sum(aT(gIdx+1,5)==1 | aT(gIdx+1,5)==4);

% percent of growths ending in fgap, bgap, or nothing
projData.stats.percentGrowthLinkedForward=100*(f/length(gIdx));
projData.stats.percentGrowthLinkedBackward=100*(b/length(gIdx));
projData.stats.percentGrowthTerminal=100*(u/length(gIdx));

% calculate dynamicity
% all track indices where there is either a forward or backward gap
idx=unique([projData.tracksWithFgap; projData.tracksWithBgap]);
idxCell=mat2cell(idx,ones(length(idx),1),1);
% sub track indices of the full tracks
subIdx=cellfun(@(x) find(aT(:,1)==x),idxCell,'uniformoutput',0);
% full track lifetimes and displacements
ltf=cell2mat(cellfun(@(x) sum(aT(x,6)),subIdx,'uniformoutput',0));
disp=cell2mat(cellfun(@(x) sum(abs(aT(x,7))),subIdx,'uniformoutput',0));
% collective displacement of all gap-containing MTs over their collective lifetime
projData.stats.dynamicity=((projData.pixSizeNm/1000)*sum(disp))/((projData.secPerFrame/60)*sum(ltf));


projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix=aT;

% save each projData in its own directory
save([runInfo.metaDir filesep 'projData'],'projData')




