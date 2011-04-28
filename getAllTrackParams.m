function [allTrackParams, tFirst, lifetime] = getAllTrackParams(movieData)

nFrames = movieData.nImages(1);

% Load feature and track infos
load(fullfile(movieData.particleDetection.directory, movieData.particleDetection.filename));
load(fullfile(movieData.particleTracking.directory, movieData.particleTracking.filename));

% Check there is no split and merge. If split and merge is enabled, it
% would complexify the interpolation in gaps.

seqOfEvents = vertcat(tracksFinal(:).seqOfEvents); %#ok<NODEF>
assert(nnz(isnan(seqOfEvents(:,4))) == size(seqOfEvents,1));

SEL = getTrackSEL(tracksFinal);

% Remove any 1-frame long track.
isValid = SEL(:,3) ~= 1;
tracksFinal = tracksFinal(isValid);
SEL = SEL(isValid,:);

tFirst = SEL(:,1);                   % first frame of the track
lifetime = SEL(:,3);                 % lifetime of the track

% Linearize frame indices of all the tracks
allTrackFrameIdx = arrayfun(@(a,b) a:b, tFirst, tFirst + lifetime - 1, 'UniformOutput', false);
allTrackFrameIdx = [allTrackFrameIdx{:}];

% allTrackParams: feature parameters of all tracks appended together
% x11 y11 A11 Sx11 Sy11 T11 C11
% x12 ... (2nd point of the 1st track)
% ...
% x1n ... (nth point of the 1st track)
% x21 ... (1st point of the 2nd track)
% ...
allTrackParams = nan(numel(allTrackFrameIdx),7);

tracksFeatIndx = [tracksFinal(:).tracksFeatIndxCG];

for iFrame = 1:nFrames
    ind = allTrackFrameIdx == iFrame & tracksFeatIndx ~= 0;
    xCoord = featuresInfo(iFrame).xCoord(:,1);
    yCoord = featuresInfo(iFrame).yCoord(:,1);
    amp = featuresInfo(iFrame).amp(:,1);
    sX = featuresInfo(iFrame).stdAlong(:,1);
    sY = featuresInfo(iFrame).stdAside(:,1);
    theta = featuresInfo(iFrame).theta(:,1);
    bkg = featuresInfo(iFrame).bkg(:,1);
    
    allTrackParams(ind,1) = xCoord(tracksFeatIndx(ind));
    allTrackParams(ind,2) = yCoord(tracksFeatIndx(ind));
    allTrackParams(ind,3) = amp(tracksFeatIndx(ind));
    allTrackParams(ind,4) = sX(tracksFeatIndx(ind));
    allTrackParams(ind,5) = sY(tracksFeatIndx(ind));
    allTrackParams(ind,6) = theta(tracksFeatIndx(ind));
    allTrackParams(ind,7) = bkg(tracksFeatIndx(ind));
end

% Interpolate feature parameters in gaps
gacombIdx = diff(isnan(allTrackParams(:,1)));
gapStarts = find(gacombIdx==1)+1;
gapEnds = find(gacombIdx==-1);
gapLengths = gapEnds-gapStarts+1;
nGaps = length(gapLengths);

for iGap = 1:nGaps
    borderIdx = [gapStarts(iGap)-1 gapEnds(iGap)+1];
    gacombIdx = gapStarts(iGap):gapEnds(iGap);
    allTrackParams(gacombIdx,1) = interp1(borderIdx, allTrackParams(borderIdx,1), gacombIdx);
    allTrackParams(gacombIdx,2) = interp1(borderIdx, allTrackParams(borderIdx,2), gacombIdx);
    allTrackParams(gacombIdx,3) = interp1(borderIdx, allTrackParams(borderIdx,3), gacombIdx);
    allTrackParams(gacombIdx,4) = interp1(borderIdx, allTrackParams(borderIdx,4), gacombIdx);
    allTrackParams(gacombIdx,5) = interp1(borderIdx, allTrackParams(borderIdx,5), gacombIdx);
    allTrackParams(gacombIdx,6) = interp1(borderIdx, allTrackParams(borderIdx,6), gacombIdx);
    allTrackParams(gacombIdx,7) = interp1(borderIdx, allTrackParams(borderIdx,7), gacombIdx);
end

assert(all(~isnan(allTrackParams(:))));