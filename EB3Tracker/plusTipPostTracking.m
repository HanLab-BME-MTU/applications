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
%       vector containing the difference in displacement between any two
%       frame pairs from the initial growth trajectories (before gap
%       closing) - provides a roughly normal distribution, the width of
%       which indicates how much intra-track speed variation there is
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
%       converted to microns/minute
%  .tracksWithFgap
%       track numbers corresponding to tracks that contain forward gaps
%       (pause or out of focus events)
%  .tracksWithBgap
%       track numbers corresponding to tracks that contain backward gaps
%       (shrinkage events)
%  .tracksWithUnclassified
%       track numbers corresponding to tracks that contain unclassified gaps
%  .typeStats
%       structure containing the following fields
%           .type1_median_micPerMin : median of all growth average velocities
%           .type1_mean_micPerMin   : mean of all growth average velocities
%           .type1_std_micPerMin    : std of all growth average velocities
%           .type1_Ppause           : probability of pause (1 over the
%                                     average total time (in seconds) spent
%                                     growing prior to pause event
%           .type1_Pcat             : probability of shrinkage (1 over the
%                                     average total time (in seconds) spent
%                                     growing prior to catastrophe
%           .type2_median_micPerMin : median of all pause average velocities
%           .type2_mean_micPerMin   : mean of all pause average velocities
%           .type2_std_micPerMin    : std of all pause average velocities
%           .type3_median_micPerMin : median of all shrinkage average velocities
%           .type3_mean_micPerMin   : mean of all shrinkage average velocities
%           .type3_std_micPerMin    : std of all shrinkage average velocities
%           .type4_median_micPerMin : median of all unclassified average velocities
%           .type4_mean_micPerMin   : mean of all unclassified average velocities
%           .type4_std_micPerMin    : std of all unclassified average velocities
%  .nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix
%       matrix containing the profile of all the tracks, where the columns
%       represent:
%           1. track number
%           2. start frame
%           3. end frame
%           4. velocity (microns/min)
%           5. subtrack type
%              (1=growth, 2=pause, 3=shrinkage, 4=unclassifed gap)
%           6. length in frames
%           7. total displacement (pixels)



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
[trackedFeatureInfo,trackedFeatureInfoInterp,trackInfo,trackVelocities]=...
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
[projData.segGapAvgVel_micPerMin]=pixPerFrame2umPerMin(trackVelocities.segmentAvgs,secPerFrame,pixSizeNm);

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
pauseIdx=find(aT(:,5)==2);
beforePauseIdx=pauseIdx-1;

% these are the pauses to consolidate
growthPauseIdx=pauseIdx(aT(pauseIdx,4)>0.9.*aT(beforePauseIdx,4));
% these are the affected track numbers
tracks2check=unique(aT(growthPauseIdx,1));
aTrows2remove=[];
for i=1:length(tracks2check)
    % rows of aT corresponding to i track
    subIdx=find(aT(:,1)==tracks2check(i));
    % rows of aT corresponding to pause events i track to consolidate
    pause2remIdx=intersect(growthPauseIdx,subIdx);
    % rows of aT corresonding to shrinkage events or real pauses in i track
    sepIdx=union(subIdx(aT(subIdx,5)==3),setdiff(intersect(subIdx,pauseIdx),pause2remIdx))';
    % split the track based on shrinkage events, so that all pauses that
    % should be consolidated together can be done at the same time
    sIdx=[subIdx(1); sepIdx(:)];
    eIdx=[sepIdx(:); subIdx(end)];
    % loop through groups of subtracks (split by shrinkage)
    for j=1:length(sIdx)
        % pTemp contains the ones to consolidate in this section
        pTemp=intersect(pause2remIdx,sIdx(j):eIdx(j));
        if ~isempty(pTemp)
            fIdx=min(pTemp)-1; % first row - prior to first pause
            lIdx=max(pTemp)+1; % last row - after final pause

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

% figure out which track numbers contain a pause, catastrophe, or
% unclassified gap
projData.tracksWithFgap       = unique(aT(aT(:,5)==2,1));
projData.tracksWithBgap = unique(aT(aT(:,5)==3,1));
projData.tracksWithUnclassifed = unique(aT(aT(:,5)==4,1));

% convert pix/frame to micron/min velocities
[aT(:,4)] = pixPerFrame2umPerMin(aT(:,4),secPerFrame,pixSizeNm);

% median/mean/std for growth speed (microns per minute)
projData.typeStats.type1_median_micPerMin = median(aT(aT(:,5)==1,4));
projData.typeStats.type1_mean_micPerMin = mean(aT(aT(:,5)==1,4));
projData.typeStats.type1_std_micPerMin  =  std(aT(aT(:,5)==1,4));

% probability of pausing is 1 over the average total time (in seconds) spent
% growing prior to pause event
beforePauseIdx=find(aT(:,5)==2)-1;
beforePauseIdx=beforePauseIdx(aT(beforePauseIdx,2)~=1);
if isempty(beforePauseIdx)
    projData.typeStats.type1_Ppause=NaN;
else
    projData.typeStats.type1_Ppause=mean(1./(aT(beforePauseIdx,6).*secPerFrame));
end

% probability of shrinking is 1 over the average total time (in seconds) spent
% growing prior to shrinkage event
beforeShrinkIdx=find(aT(:,5)==3)-1;
beforeShrinkIdx=beforeShrinkIdx(aT(beforeShrinkIdx,2)~=1);
if isempty(beforeShrinkIdx)
    projData.typeStats.type1_Pcat=NaN;
else
    projData.typeStats.type1_Pcat=mean(1./(aT(beforeShrinkIdx,6).*secPerFrame));
end

% median/mean/std for pause speed (microns per minute)
projData.typeStats.type2_median_micPerMin = median(aT(aT(:,5)==2,4));
projData.typeStats.type2_mean_micPerMin = mean(aT(aT(:,5)==2,4));
projData.typeStats.type2_std_micPerMin  =  std(aT(aT(:,5)==2,4));

% median/mean/std for shrinkage speed (microns per minute)
projData.typeStats.type3_median_micPerMin = median(aT(aT(:,5)==3,4));
projData.typeStats.type3_mean_micPerMin = mean(aT(aT(:,5)==3,4));
projData.typeStats.type3_std_micPerMin  =  std(aT(aT(:,5)==3,4));

% median/mean/std for unclassified speed (microns per minute)
projData.typeStats.type4_median_micPerMin = median(aT(aT(:,5)==4,4));
projData.typeStats.type4_mean_micPerMin = mean(aT(aT(:,5)==4,4));
projData.typeStats.type4_std_micPerMin  =  std(aT(aT(:,5)==4,4));


projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix=aT;

% save each projData in its own directory
save([runInfo.metaDir filesep 'projData'],'projData')




