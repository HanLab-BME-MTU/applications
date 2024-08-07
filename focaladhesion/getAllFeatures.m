function [tracksFinal, allFeatures, tFirst, lifetime] = getAllFeatures(featuresInfo, tracksFinal, minLifetime)
%This function formats information from the tracker: Tracks starting point,
%ending point and lifetime.
%
%Usage: 
%      [tracksFinal, allFeatures, tFirst, lifetime] = getAllFeatures(featuresInfo, trackPath, minLifetime)
%
% Input:
%       featuresInfo - structure
%       trackFinal   - structure   
%       minLifetime  - tracks shorter than this value are thrown away
%
% Output:
%       tracksFinal
%       allFeatures
%       tFirst
%       lifetime
%
% See also:
%
%This is a slightly modified version of Sylvain's function: getAllFeatures.m
%See vinculinFA project for the original function
%
%Marco Vilela, 2012
ip = inputParser;
ip.addRequired('featuresInfo',@isstruct);
ip.addRequired('trackPath',   @isstruct);
ip.addRequired('minLifetime', @isscalar);

ip.parse(featuresInfo,tracksFinal,minLifetime);


nFrames = length(featuresInfo);

% Load feature and track infos


% Check there is no split and merge. If split and merge is enabled, it
% would complexify the interpolation in gaps.

seqOfEvents = vertcat(tracksFinal(:).seqOfEvents); 
assert(nnz(isnan(seqOfEvents(:,4))) == size(seqOfEvents,1));

SEL = getTrackSEL(tracksFinal);

% Remove any 1-frame long track.
isValid = SEL(:,3) >= minLifetime;
tracksFinal = tracksFinal(isValid);
SEL = SEL(isValid,:);

tFirst = SEL(:,1);                   % first frame of the track
lifetime = SEL(:,3);                 % lifetime of the track

% Linearize frame indices of all the tracks
allTrackFrameIdx = arrayfun(@(a,b) a:b, tFirst, tFirst + lifetime - 1, 'UniformOutput', false);
allTrackFrameIdx = [allTrackFrameIdx{:}];

% allFeatures: feature parameters of all tracks appended together
% x11 y11 A11 Sx11 Sy11 T11 C11
% x12 ... (2nd point of the 1st track)
% ...
% x1n ... (nth point of the 1st track)
% x21 ... (1st point of the 2nd track)
% ...
allFeatures = nan(numel(allTrackFrameIdx),7);

tracksFeatIndx = [tracksFinal(:).tracksFeatIndxCG];

for iFrame = 1:nFrames
    ind = allTrackFrameIdx == iFrame & tracksFeatIndx ~= 0;
    xCoord = featuresInfo(iFrame).xCoord(:,1);
    yCoord = featuresInfo(iFrame).yCoord(:,1);
    amp = featuresInfo(iFrame).amp(:,1);
    sX = featuresInfo(iFrame).sigmaX(:,1);
    sY = featuresInfo(iFrame).sigmaY(:,1);
    theta = featuresInfo(iFrame).theta(:,1);
    bkg = featuresInfo(iFrame).bkg(:,1);
    
    allFeatures(ind,1) = xCoord(tracksFeatIndx(ind));
    allFeatures(ind,2) = yCoord(tracksFeatIndx(ind));
    allFeatures(ind,3) = amp(tracksFeatIndx(ind));
    allFeatures(ind,4) = sX(tracksFeatIndx(ind));
    allFeatures(ind,5) = sY(tracksFeatIndx(ind));
    allFeatures(ind,6) = theta(tracksFeatIndx(ind));
    allFeatures(ind,7) = bkg(tracksFeatIndx(ind));
end

% Interpolate feature parameters in gaps
gacombIdx = diff(isnan(allFeatures(:,1)));
gapStarts = find(gacombIdx==1)+1;
gapEnds = find(gacombIdx==-1);
gapLengths = gapEnds-gapStarts+1;
nGaps = length(gapLengths);

for iGap = 1:nGaps
    borderIdx = [gapStarts(iGap)-1 gapEnds(iGap)+1];
    gacombIdx = gapStarts(iGap):gapEnds(iGap);
    allFeatures(gacombIdx,1) = interp1(borderIdx, allFeatures(borderIdx,1), gacombIdx);
    allFeatures(gacombIdx,2) = interp1(borderIdx, allFeatures(borderIdx,2), gacombIdx);
    allFeatures(gacombIdx,3) = interp1(borderIdx, allFeatures(borderIdx,3), gacombIdx);
    allFeatures(gacombIdx,4) = interp1(borderIdx, allFeatures(borderIdx,4), gacombIdx);
    allFeatures(gacombIdx,5) = interp1(borderIdx, allFeatures(borderIdx,5), gacombIdx);
    allFeatures(gacombIdx,6) = interp1(borderIdx, allFeatures(borderIdx,6), gacombIdx);
    allFeatures(gacombIdx,7) = interp1(borderIdx, allFeatures(borderIdx,7), gacombIdx);
end

assert(all(~isnan(allFeatures(:))));
