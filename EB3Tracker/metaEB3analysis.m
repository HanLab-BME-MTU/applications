function [projData]=metaEB3analysis(runInfo,secPerFrame,pixSizeNm)
% METAEB3ANALYSIS extracts various statistics from EB3 tracks
%
%SYNOPSIS [projData]=metaEB3analysis(runInfo,secPerFrame,pixSizeNm)
%
%INPUT  runInfo           : structure containing fields .anDir, which gives
%                           the full path to the roi_x directory
%                           and .imDir, which gives the full path to the
%                           folder containing the images for overlay.
%                           if given as [], program will query user for
%                           roi_x directory.
%       secPerFrame       : frame rate (seconds per frame)
%       pixSizeNm         : real-space pixel size in nanometers
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
%       area in pixels of feature from detection (for segs only)
%  .featInt
%       max intensity of feature from detection (for segs only)
%  .frame2frameDispPix
%       vector containing all frame-to-frame displacements from the initial
%       segments (not gaps) of all tracks
%  .pair2pairDiffPix
%       vector containing the difference in displacement between any two 
%       frame pairs from the initial segments (not gaps) - provides a
%       roughly normal distribution, the width of which indicates how much
%       intra-track speed variation there is
%  .medNNdistWithinFramePix
%       median nearest-neighbor distance between all detected features in a
%       frame from movieInfo (though loops through all frames to get
%       composite from full movie, not just first frame)
%  .meanDisp2medianNNDistRatio
%       mean displacement / median NN - 0.3 is upper bound tested in
%       Jaqaman Nature Methods 2008; this statistic is a measure of how
%       challenging the tracking is
%  .frame2frameVel_micPerMin
%       intantaneous velocity based on interpolated positions during gaps,
%       converted to microns/minute
%  .segGapAvgVel_micPerMin
%       average velocity based on interpolated positions during gaps,
%       converted to microns/minute
%  .tracksWithFgaps
%       track numbers corresponding to tracks that contain forward gaps
%  .tracksWithBgaps
%       track numbers corresponding to tracks that contain backward gaps
%       (shrinkage events)
%  .tracksWithUgaps
%       track numbers corresponding to tracks that contain unclassified gaps
%  .nTrack_start_end_velMicPerMin_class_lifetime
%       matrix containing the profile of all the tracks, where the columns
%       represent:
%           1. track number 
%           2. start frame 
%           3. end frame 
%           4. velocity (microns/min)
%           5. segment or gap type
%              (1=growth, 2=forward gap/pause, 3=backward gap/shrinkage,
%              4=unclassifed gap)
%           6. length in frames
%  .segGapMeanStd_micPerMin
%       structure containing the following fields
%           .meanSegVel  : mean of all segment average velocities
%           .stdSegVel   : standard deviation of all segment avg velocities
%           .meanFgapVel : same, but for forward gaps
%           .stdFgapVel  : etc. for backward and unclassified gaps...

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
        error('--metaEB3analysis: first argument should be a structure with fields imDir and anDir');
    else
        [runInfo.anDir] = formatPath(runInfo.anDir);
        [runInfo.imDir] = formatPath(runInfo.imDir);
    end
end


% load movieInfo (detection result)
featDir  = [runInfo.anDir filesep 'feat'];
if ~isdir(featDir)
    error('--metaEB3analysis: feat directory missing')
else
    if exist([featDir filesep 'movieInfo.mat'])
        load([featDir filesep 'movieInfo.mat'])
    else
        error('--metaEB3analysis: movieInfo missing...')
    end
end

% load tracksFinal (tracking result)
trackDir = [runInfo.anDir filesep 'track'];
if ~isdir(trackDir)
    error('--metaEB3analysis: track directory missing')
else
    [listOfFiles]=searchFiles('.mat',[],trackDir,0);
    if ~isempty(listOfFiles)
        load([listOfFiles{1,2} filesep listOfFiles{1,1}])
        if ~exist('tracksFinal','var')
            error('--metaEB3analysis: tracksFinal missing...');
        end
    else
        error('--metaEB3analysis: tracksFinal missing...');
    end
end

% make a new "meta" directory if it does not exist
runInfo.metaDir = [runInfo.anDir filesep 'meta'];
if ~isdir(runInfo.metaDir)
    mkdir(runInfo.metaDir)
end

% convert tracksFinal to matrix
if isstruct(tracksFinal)
    [trackedFeatureInfo,trackedFeatureIndx] = convStruct2MatNoMS(tracksFinal);
    %clear trackedFeatureIndx
end

% get interpolated positions for gaps and calculate velocities
[trackedFeatureInfo,trackedFeatureInfoInterp,trackInfo,trackVelocities] = getVelocitiesFromMat(tracksFinal,3);

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


% save misc info for output
projData.imDir = runInfo.imDir;
projData.anDir = runInfo.anDir;
projData.secPerFrame = secPerFrame;
projData.pixSizeNm = pixSizeNm;
projData.numTracks = numTracks;
projData.numFrames = numTimePoints;
projData.xCoord = trackedFeatureInfoInterp(:,1:8:end); 
projData.yCoord = trackedFeatureInfoInterp(:,2:8:end);
projData.featArea = featArea;
projData.featInt = featInt;


% get frame-to-frame displacement for segments only (not gaps)
frame2frameDispPix=sqrt(diff(x,1,2).^2+diff(y,1,2).^2);
% get rid of NaNs and linearize the vector
projData.frame2frameDispPix=frame2frameDispPix(~isnan(frame2frameDispPix(:)));

% get change in velocity between frame *pairs* for segments only
pair2pairDiffPix=diff(frame2frameDispPix,1,2);
% get rid of NaNs and linearize the vector
projData.pair2pairDiffPix=pair2pairDiffPix(~isnan(pair2pairDiffPix(:)));


% get all feature nearest neighbor distances from all frames in one vector
NNdist=zeros(length(vertcat(movieInfo.xCoord)),1);
count=1;
for iFrame=1:length(movieInfo)

    xCoord = movieInfo(iFrame).xCoord(:,1);
    yCoord = movieInfo(iFrame).yCoord(:,1);

    D=createDistanceMatrix([xCoord yCoord],[xCoord yCoord]);
    [sD,idx]=sort(D,2);

    NNdist(count:count+length(xCoord)-1)=sD(:,2);

    count=count+length(xCoord);
end

% median NN dist
projData.medNNdistWithinFramePix=median(NNdist);

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
    segs =  [segs   1*ones(size(segs,1),1)   segs(:,3)-segs(:,2)];
    
    compositeMatrix = [compositeMatrix; segs];
    % get mean/std speed (microns/min) for whole population of segs (taken from seg averages)
    projData.segGapMeanStd_micPerMin.meanSegVel = mean(pixPerFrame2umPerMin(segs(:,4),secPerFrame,pixSizeNm));
    projData.segGapMeanStd_micPerMin.medianSegVel = median(pixPerFrame2umPerMin(segs(:,4),secPerFrame,pixSizeNm));
    projData.segGapMeanStd_micPerMin.stdSegVel  =  std(pixPerFrame2umPerMin(segs(:,4),secPerFrame,pixSizeNm));
    
end
if ~isempty(fgaps)
    fgaps = [fgaps  2*ones(size(fgaps,1),1) fgaps(:,3)-fgaps(:,2)];
    % get indices of tracks where there are forward gaps
    projData.tracksWithFgaps = unique(fgaps(:,1));
    
    compositeMatrix = [compositeMatrix; fgaps];
    % get mean/std speed (microns/min) for whole population of fgaps
    projData.segGapMeanStd_micPerMin.meanFgapVel = mean(pixPerFrame2umPerMin(fgaps(:,4),secPerFrame,pixSizeNm));
    projData.segGapMeanStd_micPerMin.stdFgapVel  =  std(pixPerFrame2umPerMin(fgaps(:,4),secPerFrame,pixSizeNm));
    
end
if ~isempty(bgaps)
    bgaps = [bgaps  3*ones(size(bgaps,1),1) bgaps(:,3)-bgaps(:,2)];
    % get indices of tracks where there are backward gaps
    projData.tracksWithBgaps = unique(bgaps(:,1));
    
    compositeMatrix = [compositeMatrix; bgaps];
    % get mean/std speed (microns/min) for whole population of bgaps
    projData.segGapMeanStd_micPerMin.meanBgapVel = mean(pixPerFrame2umPerMin(bgaps(:,4),secPerFrame,pixSizeNm));
    projData.segGapMeanStd_micPerMin.stdBgapVel  =  std(pixPerFrame2umPerMin(bgaps(:,4),secPerFrame,pixSizeNm));

end
if ~isempty(ugaps)
    ugaps = [ugaps  4*ones(size(ugaps,1),1) ugaps(:,3)-ugaps(:,2)];
    % get indices of tracks where there are unclassified gaps
    projData.tracksWithUgaps = unique(ugaps(:,1));
    
    compositeMatrix = [compositeMatrix; ugaps];
    % get mean/std speed (microns/min) for whole population of ugaps
    projData.segGapMeanStd_micPerMin.meanUgapVel = mean(pixPerFrame2umPerMin(ugaps(:,4),secPerFrame,pixSizeNm));
    projData.segGapMeanStd_micPerMin.stdUgapVel  =  std(pixPerFrame2umPerMin(ugaps(:,4),secPerFrame,pixSizeNm));

end
    

% put segs/gaps into one matrix and sort to see track profiles in order
aT=sortrows(compositeMatrix,[1 2]);

% convert pix/frame --> micron/min velocities
[aT(:,4)] = pixPerFrame2umPerMin(aT(:,4),secPerFrame,pixSizeNm);
projData.nTrack_start_end_velMicPerMin_class_lifetime=aT;


% get stats for each track in the form:
% [trackNumber trackLength sum(nSegs+nGaps) totalSegFrames totalFgapFrames totalBgapFrames totalUgapFrames]
% tS=zeros(numTracks,7);
% for iTrack=1:numTracks
%     tS(iTrack,1)=iTrack;
%     idx=find(aT(:,1)==iTrack); % index of those rows corresponding to iTrack
%     tS(iTrack,2)=max(aT(idx,3))-min(aT(idx,2)); % total track length in frames
%     tS(iTrack,3)=length(idx); % number of segs and gaps that make up the track
% 
%     tS(iTrack,4)=sum(aT(idx(aT(idx,5)==1),6)); % total number of frames spent in a segment
%     tS(iTrack,5)=sum(aT(idx(aT(idx,5)==2),6)); % total number of frames spent in a forward gap
%     tS(iTrack,6)=sum(aT(idx(aT(idx,5)==3),6)); % total number of frames spent in a backward gap
%     tS(iTrack,7)=sum(aT(idx(aT(idx,5)==4),6)); % total number of frames spent in a unclassified gap
% 
% end

%projData.trackStats=tS;

% save each projData in its own directory
save([runInfo.metaDir filesep 'projData'],'projData')




