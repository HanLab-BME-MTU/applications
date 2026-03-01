function sampleMovieOtherChannel(movieData, varargin)
% sampleMovieOtherChannel
% Sample intensity statistics from a specified channel inside a binary mask
% for each frame, and store per-cell time series by tracking objects across frames.
% Adds per-track dF/F0 and progressbar updates.
%
% INPUT
%   movieData : MovieData
%   paramsIn  : struct (optional)
%
% OUTPUT (saved to funParams.OutputDirectory)
%   otherChannelSampling.mat with struct 'S'
%
% Notes
% - Mask is expected to be a binary image per frame (logical or 0/1).
% - Per-object stats are computed using connected components or optional watershed splitting.
% - Track IDs are created by linking objects frame-to-frame using IoU overlap, with centroid fallback.
% - If StageDriftCorrectionProcess exists in TFMPackage and
%   UseStageDriftCorrection=true, masks are warped per-frame using T shifts.
% - progressbar.m must be on the MATLAB path.
%
% 2026 - Sangyoon Han lab

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn', [], @(x) isempty(x) || isstruct(x));
ip.parse(movieData, varargin{:});
paramsIn = ip.Results.paramsIn;

% Ensure process exists
iProc = movieData.getProcessIndex('OtherChannelSamplingProcess', 1, 0);
if isempty(iProc)
    iProc = numel(movieData.processes_) + 1;
    movieData.addProcess(OtherChannelSamplingProcess(movieData, ...
        fullfile(movieData.outputDirectory_, 'otherChannelSampling')));
end
proc = movieData.processes_{iProc};

% Parse params
p = parseProcessParams(proc, paramsIn);

% -------------------------
% Defaults for NEW features
% -------------------------
if ~isfield(p,'UseLabeling'); p.UseLabeling = true; end
if ~isfield(p,'MinAreaPix'); p.MinAreaPix = 200; end

% Instance segmentation from binary mask: 'cc' or 'watershed'
if ~isfield(p,'InstanceSegMethod'); p.InstanceSegMethod = 'cc'; end % stable default
if ~isfield(p,'SmoothSigma'); p.SmoothSigma = 1.5; end
if ~isfield(p,'MinSeedH'); p.MinSeedH = 2.0; end          % for imextendedmax on distance map
if ~isfield(p,'BorderClear'); p.BorderClear = false; end

% Tracking parameters
if ~isfield(p,'TrackByOverlap'); p.TrackByOverlap = true; end
if ~isfield(p,'MinIoU'); p.MinIoU = 0.10; end
if ~isfield(p,'MaxCentroidDist'); p.MaxCentroidDist = 40; end

% dF/F0 parameters (per track and/or whole mask)
if ~isfield(p,'ComputeDFF0'); p.ComputeDFF0 = false; end
if ~isfield(p,'BaselineFrames'); p.BaselineFrames = 1; end
if ~isfield(p,'BaselineFallbackN'); p.BaselineFallbackN = 5; end % if track absent in baseline frames

% Preview
if ~isfield(p,'SavePerFrameTifPreview'); p.SavePerFrameTifPreview = false; end

% progressbar
if ~isfield(p,'UseProgressBar'); p.UseProgressBar = true; end
if ~isfield(p,'ProgressBarPosition'); p.ProgressBarPosition = 4; end % bottom-right

% Output dir
if ~exist(p.OutputDirectory,'dir'); mkdir(p.OutputDirectory); end
outMat = fullfile(p.OutputDirectory, 'otherChannelSampling.mat');
proc.setOutFilePaths({outMat});

nFrames = movieData.nFrames_;
ch = p.ChannelIndex;

% --- Resolve mask process ---
maskProc = [];
if isfield(p,'MaskProcessName') && ~isempty(p.MaskProcessName)
    iMask = movieData.getProcessIndex(p.MaskProcessName, 1, 0);
    if ~isempty(iMask)
        maskProc = movieData.getProcess(iMask);
    end
else
    cand = {'MaskIntersectionProcess','MaskRefinementProcess','ThresholdProcess'};
    for i = 1:numel(cand)
        iMask = movieData.getProcessIndex(cand{i}, 1, 0);
        if ~isempty(iMask)
            maskProc = movieData.getProcess(iMask);
            p.MaskProcessName = cand{i};
            break;
        end
    end
end
if isempty(maskProc)
    error('OtherChannelSamplingProcess:NoMask', ...
        'Could not find a mask process. Set params.MaskProcessName or run a mask/threshold process first.');
end

maskChan = p.MaskChannelIndex;
if isempty(maskChan) || ~isscalar(maskChan) || maskChan < 1
    maskChan = 1;
end

% --- Optional: Stage drift correction shifts (T) ---
useSDC = isfield(p,'UseStageDriftCorrection') && p.UseStageDriftCorrection;
T = [];
existSDC = false;
if useSDC
    try
        iPack = movieData.getPackageIndex('TFMPackage');
        if ~isempty(iPack)
            pack = movieData.getPackage(iPack);
            sdcProc = pack.processes_{1}; % typically stage drift correction in TFMPackage
            if ~isempty(sdcProc) && isprop(sdcProc,'outFilePaths_')
                try
                    iBeadChan = sdcProc.funParams_.iBeadChannel;
                catch
                    iBeadChan = 1;
                end
                if size(sdcProc.outFilePaths_,1) >= 3 && size(sdcProc.outFilePaths_,2) >= iBeadChan
                    s = load(sdcProc.outFilePaths_{3,iBeadChan}, 'T');
                    T = s.T;
                    existSDC = true;
                end
            end
        end
    catch
        existSDC = false;
    end
end

ref_obj = [];
if existSDC
    ref_obj = imref2d(movieData.imSize_);
end

% --- Preallocate output ---
S = struct();
S.params = p;
S.frame = (1:nFrames)';

% Whole-mask summary stats
S.mean = nan(nFrames,1);
S.median = nan(nFrames,1);
S.areaPix = nan(nFrames,1);

% Per-frame per-object stats
S.nObjects = nan(nFrames,1);
S.objectMeans = cell(nFrames,1);
S.objectMedians = cell(nFrames,1);
S.objectAreas = cell(nFrames,1);
S.objectCentroids = cell(nFrames,1);
S.objectPixelIdxList = cell(nFrames,1);
S.frameToTrackIds = cell(nFrames,1);

% Track-level matrices (dynamic sizing)
S.nTracks = 0;
S.trackMean = nan(0, nFrames);
S.trackMedian = nan(0, nFrames);
S.trackArea = nan(0, nFrames);
S.trackCentroidX = nan(0, nFrames);
S.trackCentroidY = nan(0, nFrames);
S.trackPresent = false(0, nFrames);

% dF/F0 outputs (filled later if ComputeDFF0)
S.F0 = [];
S.dFF0 = [];
S.trackF0 = [];
S.trackdFF0 = [];

if p.SavePerFrameTifPreview
    previewDir = fullfile(p.OutputDirectory, 'preview');
    if ~exist(previewDir,'dir'); mkdir(previewDir); end
end

% Tracking state
prevL = [];
prevTrackIds = [];
prevCentroids = [];

% progressbar setup
doPB = false;
if isfield(p,'UseProgressBar') && p.UseProgressBar
    doPB = (exist('progressbar','file') == 2);
end
if doPB
    progressbar(0, p.ProgressBarPosition, 'OtherChannelSampling');
    pbCleanup = onCleanup(@() safeCloseProgressbar());
end

% --- Main loop ---
for ii = 1:nFrames

    if doPB
        progressbar((ii-1)/nFrames, p.ProgressBarPosition, sprintf('Frame %d/%d', ii, nFrames));
    end

    % Load image frame from desired channel
    I = movieData.channels_(ch).loadImage(ii);
    I = double(I);

    % Load mask (binary)
    mask = [];
    try
        mask = maskProc.loadChannelOutput(maskChan, ii);
    catch
        try
            mask = maskProc.loadChannelOutput(ii);
        catch
            mask = [];
        end
    end
    if isempty(mask)
        error('OtherChannelSamplingProcess:MaskLoadFailed', ...
            'Failed to load mask for frame %d from %s.', ii, class(maskProc));
    end
    mask = logical(mask);

    % Apply stage drift correction to mask if requested and available
    if existSDC && size(T,1) >= ii
        Tr = affine2d([1 0 0; 0 1 0; fliplr(T(ii,:)) 1]);
        mask = imwarp(mask, Tr, 'OutputView', ref_obj);
        mask = mask > 0;
    end

    % Whole-mask summary stats within mask
    pix = I(mask);
    if isempty(pix)
        S.mean(ii) = NaN;
        S.median(ii) = NaN;
        S.areaPix(ii) = 0;

        S.nObjects(ii) = 0;
        S.frameToTrackIds{ii} = [];
        continue;
    end
    S.mean(ii) = mean(pix, 'omitnan');
    S.median(ii) = median(pix, 'omitnan');
    S.areaPix(ii) = nnz(mask);

    % ------------------------------
    % Per-object stats + TRACKING
    % ------------------------------
    if isfield(p,'UseLabeling') && p.UseLabeling

        % Build label image from binary mask
        L = labelFromMask(mask, p);

        % Compute regionprops
        rp = regionprops(L, I, 'Area','Centroid','PixelIdxList','MeanIntensity');
        if isempty(rp)
            S.nObjects(ii) = 0;
            S.frameToTrackIds{ii} = [];
            continue;
        end

        % Filter by area
        areas = [rp.Area]';
        keep = areas >= p.MinAreaPix;
        rp = rp(keep);
        if isempty(rp)
            S.nObjects(ii) = 0;
            S.frameToTrackIds{ii} = [];
            continue;
        end

        nObj = numel(rp);
        S.nObjects(ii) = nObj;

        objArea = nan(nObj,1);
        objCent = nan(nObj,2);
        objMean = nan(nObj,1);
        objMed  = nan(nObj,1);
        pixIdx  = cell(nObj,1);

        for k = 1:nObj
            objArea(k) = rp(k).Area;
            objCent(k,:) = rp(k).Centroid;
            objMean(k) = rp(k).MeanIntensity;
            pixIdx{k} = rp(k).PixelIdxList;

            v = I(pixIdx{k});
            objMed(k) = median(v, 'omitnan');
        end

        % Store per-frame per-object
        S.objectAreas{ii} = objArea;
        S.objectCentroids{ii} = objCent;
        S.objectMeans{ii} = objMean;
        S.objectMedians{ii} = objMed;
        S.objectPixelIdxList{ii} = pixIdx;

        % --- Track linking (frame-to-frame) ---
        if ii == 1 || isempty(prevL) || isempty(prevTrackIds) || isempty(prevCentroids)
            currTrackIds = (S.nTracks + 1 : S.nTracks + nObj).';
            S = ensureTrackCapacity(S, max(currTrackIds), nFrames);
            S.nTracks = max(currTrackIds);
        else
            currTrackIds = linkByOverlapOrDistance( ...
                prevL, L, prevTrackIds, prevCentroids, objCent, ...
                p.MinIoU, p.MaxCentroidDist, p.TrackByOverlap);

            % New tracks for unmatched
            newMask = (currTrackIds == 0);
            if any(newMask)
                newIds = (S.nTracks + 1 : S.nTracks + nnz(newMask)).';
                currTrackIds(newMask) = newIds;
                S = ensureTrackCapacity(S, max(currTrackIds), nFrames);
                S.nTracks = max(currTrackIds);
            end
        end

        S.frameToTrackIds{ii} = currTrackIds;

        % Fill track matrices for this frame
        for k = 1:nObj
            tid = currTrackIds(k);
            S.trackMean(tid, ii) = objMean(k);
            S.trackMedian(tid, ii) = objMed(k);
            S.trackArea(tid, ii) = objArea(k);
            S.trackCentroidX(tid, ii) = objCent(k,1);
            S.trackCentroidY(tid, ii) = objCent(k,2);
            S.trackPresent(tid, ii) = true;
        end

        % Update prev state
        prevL = L;
        prevTrackIds = currTrackIds;
        prevCentroids = objCent;
    end

    % Optional preview
    if p.SavePerFrameTifPreview
        rgb = repmat(mat2gray(I),1,1,3);
        b = bwperim(mask);
        rgb(:,:,1) = max(rgb(:,:,1), b);
        imwrite(rgb, fullfile(previewDir, sprintf('frame_%03d.png', ii)));
    end
end

if doPB
    progressbar(1, p.ProgressBarPosition, 'Done');
end

% --- dF/F0 (whole-mask mean) ---
if isfield(p,'ComputeDFF0') && p.ComputeDFF0
    baseFrames = p.BaselineFrames;
    baseFrames = baseFrames(baseFrames>=1 & baseFrames<=nFrames);
    if isempty(baseFrames)
        baseFrames = 1;
    end

    % whole-mask
    F0 = mean(S.mean(baseFrames), 'omitnan');
    if isempty(F0) || ~isfinite(F0) || F0==0
        S.F0 = NaN;
        S.dFF0 = nan(nFrames,1);
    else
        S.F0 = F0;
        S.dFF0 = (S.mean - F0) ./ F0;
    end

    % per-track
    [S.trackF0, S.trackdFF0] = computeTrackDFF0(S.trackMean, S.trackPresent, baseFrames, p.BaselineFallbackN);
end

save(outMat, 'S', '-v7.3');
movieData.save();
end


% =========================
% Helper: safe close progressbar
% =========================
function safeCloseProgressbar()
try
    progressbar(1);
catch
end
end


% =========================
% Helper: label from mask
% =========================
function L = labelFromMask(maskBW, p)
maskBW = logical(maskBW);

% cleanup
maskBW = imfill(maskBW, 'holes');
maskBW = bwareaopen(maskBW, max(10, round(p.MinAreaPix/5)));

if isfield(p,'BorderClear') && p.BorderClear
    maskBW = imclearborder(maskBW);
end

method = 'cc';
if isfield(p,'InstanceSegMethod') && ~isempty(p.InstanceSegMethod)
    method = lower(p.InstanceSegMethod);
end

switch method
    case 'watershed'
        D = bwdist(~maskBW);
        sig = 0;
        if isfield(p,'SmoothSigma'); sig = p.SmoothSigma; end
        if sig > 0
            D = imgaussfilt(D, sig);
        end

        h = 2.0;
        if isfield(p,'MinSeedH'); h = p.MinSeedH; end
        seeds = imextendedmax(D, h);

        if ~any(seeds(:))
            L = bwlabel(maskBW, 8);
        else
            Dneg = -D;
            Dneg(~maskBW) = Inf;
            M = imimposemin(Dneg, seeds | ~maskBW);
            W = watershed(M);
            L = W;
            L(~maskBW) = 0;
            L(W==0) = 0;
            L = relabelDense(L);
        end

    otherwise
        L = bwlabel(maskBW, 8);
end

L = removeSmallLabelRegions(L, p.MinAreaPix);
L = relabelDense(L);
end


function L = removeSmallLabelRegions(L, minArea)
if isempty(L) || max(L(:))==0
    return;
end
stats = regionprops(L, 'Area');
areas = [stats.Area]';
keep = find(areas >= minArea);
if isempty(keep)
    L(:) = 0;
else
    L(~ismember(L, keep)) = 0;
end
end


function Ld = relabelDense(L)
u = unique(L(:));
u(u==0) = [];
Ld = zeros(size(L));
for k = 1:numel(u)
    Ld(L==u(k)) = k;
end
end


% =========================
% Helper: track containers
% =========================
function S = ensureTrackCapacity(S, nTracksNeeded, nFrames)
nNow = size(S.trackMean,1);
if nTracksNeeded <= nNow
    return;
end
addN = nTracksNeeded - nNow;

S.trackMean(end+addN, nFrames) = NaN;
S.trackMedian(end+addN, nFrames) = NaN;
S.trackArea(end+addN, nFrames) = NaN;
S.trackCentroidX(end+addN, nFrames) = NaN;
S.trackCentroidY(end+addN, nFrames) = NaN;

tmp = false(nTracksNeeded, nFrames);
if ~isempty(S.trackPresent)
    tmp(1:nNow,:) = S.trackPresent;
end
S.trackPresent = tmp;
end


% =========================
% Helper: overlap/centroid linking
% =========================
function currTrackIds = linkByOverlapOrDistance(prevL, currL, prevTrackIds, ...
    prevCentroids, currCentroids, minIoU, maxDist, useOverlap)

nCurr = size(currCentroids,1);
currTrackIds = zeros(nCurr,1);

% Ensure dense labels (1..N)
prevL = relabelDense(prevL);
currL = relabelDense(currL);

if useOverlap
    prevVec = prevL(:);
    currVec = currL(:);

    valid = prevVec > 0 & currVec > 0;
    if any(valid)
        pv = prevVec(valid);
        cv = currVec(valid);

        nPrevObj = max(prevVec);
        nCurrObj = max(currVec);

        overlapCounts = accumarray([pv, cv], 1, [nPrevObj, nCurrObj], @sum, 0);

        prevArea = accumarray(prevVec(prevVec>0), 1, [nPrevObj, 1], @sum, 0);
        currArea = accumarray(currVec(currVec>0), 1, [nCurrObj, 1], @sum, 0);

        denom = prevArea + currArea.' - overlapCounts;
        IoU = overlapCounts ./ max(denom, 1);

        [bestIoU, bestPrevIdx] = max(IoU, [], 1);

        usedPrev = false(numel(prevTrackIds),1);
        [~, order] = sort(bestIoU, 'descend');

        for t = 1:numel(order)
            j = order(t);
            if bestIoU(j) < minIoU
                break;
            end
            i = bestPrevIdx(j);
            if i < 1 || i > numel(prevTrackIds)
                continue;
            end
            if usedPrev(i)
                continue;
            end
            currTrackIds(j) = prevTrackIds(i);
            usedPrev(i) = true;
        end
    end
end

% Centroid fallback for any still-unmatched
unmatched = find(currTrackIds == 0);
if ~isempty(unmatched) && ~isempty(prevCentroids)
    D = pdist2(prevCentroids, currCentroids(unmatched,:));
    for uu = 1:numel(unmatched)
        jIdx = unmatched(uu);
        [dMin, iMin] = min(D(:,uu));
        if isfinite(dMin) && dMin <= maxDist
            candTid = prevTrackIds(iMin);
            if ~ismember(candTid, currTrackIds)
                currTrackIds(jIdx) = candTid;
            end
        end
    end
end
end


% =========================
% Helper: per-track dF/F0
% =========================
function [trackF0, trackdFF0] = computeTrackDFF0(trackMean, trackPresent, baseFrames, fallbackN)
% trackMean: [nTracks x nFrames]
% trackPresent: [nTracks x nFrames] logical
% baseFrames: vector of frames used for baseline
% fallbackN: if no baseline points exist for a track, use first fallbackN present frames

nTracks = size(trackMean,1);
nFrames = size(trackMean,2);

trackF0 = nan(nTracks,1);
trackdFF0 = nan(nTracks, nFrames);

for tid = 1:nTracks
    % baseline candidates: present in baseline frames
    bf = baseFrames(baseFrames>=1 & baseFrames<=nFrames);
    okBase = false(1,nFrames);
    okBase(bf) = true;

    idxBase = find(trackPresent(tid,:) & okBase);
    if numel(idxBase) >= 1
        F0 = mean(trackMean(tid, idxBase), 'omitnan');
    else
        % fallback: first fallbackN present frames
        idxAll = find(trackPresent(tid,:));
        if isempty(idxAll)
            F0 = NaN;
        else
            idxUse = idxAll(1:min(fallbackN, numel(idxAll)));
            F0 = mean(trackMean(tid, idxUse), 'omitnan');
        end
    end

    if isempty(F0) || ~isfinite(F0) || F0 == 0
        trackF0(tid) = NaN;
        trackdFF0(tid,:) = nan(1,nFrames);
    else
        trackF0(tid) = F0;
        trackdFF0(tid,:) = (trackMean(tid,:) - F0) ./ F0;
    end
end
end