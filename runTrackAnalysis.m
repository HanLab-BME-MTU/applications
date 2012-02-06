%[data] = runTrackAnalysis(data, varargin)
%
% INPUTS    data        : array of experiment structures
%           {Buffer}    : length of buffer before and after tracks
%           {Overwrite} : 1 to overwrite previous results
%           {FileName}  : name of the output
%
% Outputs
%
%           nSeg : number of segments; 1: regular, 2 or more: merging/splitting
%         status : 1) complete track, 2) incomplete (cut off), 3) persistent
%      gapStatus : 4) valid gap, 5) invalid gap
%          valid : '1' if status == 1 and all gaps have status 4;
%
% Usage example: runTrackAnalysis(data, 'Buffer', 10);
%
% Note: Only tracks with track.status==1 and track.gapStatus==4 are considered.

% Francois Aguet, November 2010 (last modified 02/01/2012)

function runTrackAnalysis(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParamValue('Buffer', 5, @isscalar);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('TrackerOutput', 'trackedFeatures.mat', @ischar);
ip.addParamValue('FileName', 'trackAnalysis.mat', @ischar);
ip.addParamValue('FrameIndexes', arrayfun(@(x) 1:x.movieLength, data, 'UniformOutput', false), @(x) numel(unique(diff(x)))==1); %check that frame rate is constant
ip.parse(data, varargin{:});
filename = ip.Results.FileName;
overwrite = ip.Results.Overwrite;
buffer = ip.Results.Buffer;
trackerOutput = ip.Results.TrackerOutput;
frameIdx = ip.Results.FrameIndexes;

for i = 1:length(data)
    data(i).tracks = [];
    data(i).smTracks = [];
end
parfor i = 1:length(data)
    if ~(exist([data(i).source filesep 'Tracking' filesep filename],'file')==2) || overwrite
        data(i) = main(data(i), buffer, trackerOutput, filename, frameIdx{i});
    else
        fprintf('TrackAnalysis has already been run for: %s\n', getShortPath(data(i)));
    end
end



function [data] = main(data, buffer, trackerOutput, filename, frameIdx)

detection = load([data.source 'Detection' filesep 'detection_v2.mat']);
frameInfo = detection.frameInfo;

ny = data.imagesize(1);
nx = data.imagesize(2);
nFrames = length(frameIdx);

alpha = 0.05;
kLevel = norminv(1-alpha/2.0, 0, 1); % ~2 std above background

%=================================
% Identify master/slave channels
%=================================
nCh = length(data.channels);
mCh = strcmp(data.source, data.channels);

sigmaV = zeros(nCh, 1);
for k = 1:nCh
    sigmaV(k) = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, data.markers{k});
    data.framePaths{k} = data.framePaths{k}(frameIdx);
end

if isempty(data.maskPaths)
    maskPath = [data.source 'Detection' filesep 'Masks' filesep];
    data.maskPaths = arrayfun(@(x) [maskPath x.name], dir([maskPath '*.tif']), 'UniformOutput', false);
end
data.maskPaths = data.maskPaths(frameIdx);

data.framerate = data.framerate*(frameIdx(2)-frameIdx(1));
data.movieLength = length(data.framePaths{1});

sigma = sigmaV(mCh);
w2 = ceil(2*sigma);
w3 = ceil(3*sigma);
w4 = ceil(4*sigma);

[x,y] = meshgrid(-w4:w4);
r = sqrt(x.^2+y.^2);
annularMask = zeros(size(r));
annularMask(r<=w4 & r>=w3) = 1;

%=================================
% Read and convert tracker output
%=================================
tPath = [data.source 'Tracking' filesep trackerOutput];
if exist(tPath, 'file')==2
    trackinfo = load(tPath);
    trackinfo = trackinfo.tracksFinal;
    nTracks = length(trackinfo);  
elseif exist([data.source 'TrackInfoMatrices' filesep 'trackedFeatures.mat'], 'file')==2
    % (for old tracker. oldest version: trackInfo.mat)
    trackinfo = load([data.source 'TrackInfoMatrices' filesep 'trackedFeatures.mat']);
    %trackedFeatureNum = trackinfo.trackedFeatureNum;
    trackinfo = trackinfo.trackedFeatureInfo;
    nTracks = size(trackinfo, 1);
else
    error('No valid tracker output found.');
end

tracks(1:nTracks) = struct('t', [], 'f', [],...
    'x', [], 'y', [], 'A', [], 'c', [],...
    'x_pstd', [], 'y_pstd', [], 'A_pstd', [], 'c_pstd', [],...
    'sigma_r', [], 'SE_sigma_r', [],...
    'pval_Ar', [], 'pval_KS', [], 'isPSF', [],...
    'tracksFeatIndxCG', [], 'gapVect', [], 'gapStatus', [], 'seqOfEvents', [],...
    'nSeg', [], 'visibility', [], 'lifetime_s', [], 'start', [], 'end', [],...
    'startBuffer', [], 'endBuffer', [], 'MotionParameters', []);
    %    'alphaMSD', [], 'MSD', [], 'MSDstd', [], 'totalDisplacement', [], 'D', [], ...

% field names with multiple channels
mcFieldNames = {'x', 'y', 'A', 'c', 'x_pstd', 'y_pstd', 'A_pstd', 'c_pstd', 'sigma_r', 'SE_sigma_r', 'pval_Ar', 'pval_KS', 'isPSF'};
bufferFieldNames = {'t', 'x', 'y', 'A', 'c', 'A_pstd', 'sigma_r', 'SE_sigma_r', 'pval_Ar'};

%==============================
% Loop through tracks
%==============================
fprintf('TrackAnalysis - Converting tracker output:     ');
for k = 1:nTracks

    % convert/assign structure fields
    seqOfEvents = trackinfo(k).seqOfEvents;
    tracksFeatIndxCG = trackinfo(k).tracksFeatIndxCG; % index of the feature in each frame    
    nSeg = size(tracksFeatIndxCG,1);

    segLengths = NaN(1,nSeg);

    % discarding rules: single frame segments w/ merge/split. ADD: longer segments with merge & split
    msIdx = NaN(1,nSeg);
    for s = 1:nSeg
        idx = seqOfEvents(:,3)==s;
        ievents = seqOfEvents(idx, :);
        bounds = ievents(:,1); % beginning & end of this segment
        if ~isnan(ievents(2,4))
            bounds(2) = bounds(2)-1; % correction if end is a merge
        end
        segLengths(s) = bounds(2)-bounds(1)+1;
        
        % determine whether segment has merged and split from same master, if length < 4
        msIdx(s) = segLengths(s)==1 || (diff(ievents(:,4))==0 && segLengths(s)<4);
    end
    if nSeg>1
        segIdx = find(msIdx==0);
        nSeg = numel(segIdx);
        msIdx = find(msIdx);
        if ~isempty(msIdx)
            tracksFeatIndxCG(msIdx,:) = [];
            for i = 1:numel(msIdx)
                seqOfEvents(seqOfEvents(:,3)==msIdx(i),:) = [];
            end
        end
        segLengths = segLengths(segIdx);
    else
        segIdx = 1;
    end
    
    tracks(k).nSeg = nSeg;
    firstIdx = trackinfo(k).seqOfEvents(1,1);
    lastIdx = trackinfo(k).seqOfEvents(end,1);
    
    tracks(k).lifetime_s = (lastIdx-firstIdx+1)*data.framerate;
    tracks(k).start = firstIdx;
    tracks(k).end = lastIdx;
    
    tracks(k).seqOfEvents = seqOfEvents;
    tracks(k).tracksFeatIndxCG = tracksFeatIndxCG; % index of the feature in each frame
    
    if (1<tracks(k).start) && (tracks(k).end<nFrames) % complete tracks
        tracks(k).visibility = 1;
    elseif tracks(k).start==1 && tracks(k).end==nFrames % persistent tracks
        tracks(k).visibility = 3;
    else
        tracks(k).visibility = 2; % incomplete tracks
    end
    
    %==============================================================================
    % Initialize arrays
    %==============================================================================
    
    % Segments are concatenated into single arrays, separated by NaNs.
    fieldLength = sum(segLengths)+nSeg-1;
    for f = 1:length(mcFieldNames)
        tracks(k).(mcFieldNames{f}) = NaN(nCh, fieldLength);
    end
    tracks(k).t = NaN(1, fieldLength);
    tracks(k).f = NaN(1, fieldLength);
    
    if fieldLength>1
        
        % start buffer size for this track
        sb = firstIdx - max(1, firstIdx-buffer);
        eb = min(lastIdx+buffer, data.movieLength)-lastIdx;
        
        if sb>0 && tracks(k).visibility==1
            for f = 1:length(bufferFieldNames)
                tracks(k).startBuffer.(bufferFieldNames{f}) = NaN(nCh, sb);
            end
        end
        if eb>0 && tracks(k).visibility==1
            for f = 1:length(bufferFieldNames)
                tracks(k).endBuffer.(bufferFieldNames{f}) = NaN(nCh, eb);
            end
        end
    end
    
    %==============================================================================
    % Read amplitude & background from detectionResults.mat (localization results)
    %==============================================================================
    delta = [0 cumsum(segLengths(1:end-1))+(1:nSeg-1)];

    for s = 1:nSeg
        ievents = seqOfEvents(seqOfEvents(:,3)==segIdx(s), :);
        bounds = ievents(:,1);
        if ~isnan(ievents(2,4))
            bounds(2) = bounds(2)-1;
        end
        
        nf = bounds(2)-bounds(1)+1;
        frameRange = frameIdx(bounds(1):bounds(2)); % relative to movie (also when movie is subsampled)
        
        for i = 1:length(frameRange)
            idx = tracksFeatIndxCG(s, frameRange(i) - tracks(k).start + 1); % -> relative to IndxCG
            if idx ~= 0 % if not a gap, get detection values
                for f = 1:length(mcFieldNames)
                    tracks(k).(mcFieldNames{f})(:,i+delta(s)) = frameInfo(frameRange(i)).(mcFieldNames{f})(:,idx);
                end
            end
        end
        tracks(k).t(delta(s)+(1:nf)) = (bounds(1)-1:bounds(2)-1)*data.framerate;
        tracks(k).f(delta(s)+(1:nf)) = frameRange;
    end
    tracks(k).pval_Ar = 1-tracks(k).pval_Ar; % bug fix for detection

    fprintf('\b\b\b\b%3d%%', round(100*k/(nTracks)));
end
fprintf('\n');

% remove tracks that fall into image boundary
minx = round(arrayfun(@(t) min(t.x), tracks));
maxx = round(arrayfun(@(t) max(t.x), tracks));
miny = round(arrayfun(@(t) min(t.y), tracks));
maxy = round(arrayfun(@(t) max(t.y), tracks));

idx = minx<=w4 | miny<=w4 | maxx>nx-w4 | maxy>ny-w4;
tracks(idx) = [];

% remove single-frame tracks, store coordinates
singletons = [tracks.start]==[tracks.end];
% singleFrameTracks = tracks(singletons);
tracks(singletons) = [];

nTracks = numel(tracks);

%=======================================
% Interpolate gaps and clean up tracks
%=======================================
fprintf('TrackAnalysis - Classification:     ');
for k = 1:nTracks
    
    % gap locations in 'x' for all segments
    gapVect = isnan(tracks(k).x(mCh,:)) & ~isnan(tracks(k).t);
    tracks(k).gapVect = gapVect;

    %=================================
    % Determine track and gap status
    %=================================
    sepIdx = isnan(tracks(k).t);
    
    gapCombIdx = diff(gapVect | sepIdx);
    gapStarts = find(gapCombIdx==1)+1;
    gapEnds = find(gapCombIdx==-1);
    gapLengths = gapEnds-gapStarts+1;
    
    
    
    segmentIdx = diff([0 ~(gapVect | sepIdx) 0]); % these variables refer to segments between gaps
    segmentStarts = find(segmentIdx==1);
    segmentEnds = find(segmentIdx==-1)-1;
    segmentLengths = segmentEnds-segmentStarts+1;
    
    % loop through track segments with gaps
    %if nanIdx{s}(1)==1 || nanIdx{s}(end)==1 % temporary fix for segments that begin or end with gap
    %    tracks(k).gapStatus{s} = 5;
    %else
        
    % loop over gaps
    nGaps = numel(gapLengths);
    if nGaps>0
        gv = 1:nGaps;
        gapStatus = 5*ones(1,nGaps);
        % gap valid if segments that precede/follow are > 1 frame or if gap is a single frame
        gapStatus(segmentLengths(gv)>1 & segmentLengths(gv+1)>1 | gapLengths(gv)==1) = 4;
        
        sepIdx = sepIdx(gapStarts)==1;
        gapStatus(sepIdx) = [];
        gapStarts(sepIdx) = [];
        gapEnds(sepIdx) = [];
        nGaps = numel(gapStatus);
        
        % fill position information for valid gaps using linear interpolation
        for g = 1:nGaps
            borderIdx = [gapStarts(g)-1 gapEnds(g)+1];
            gacombIdx = gapStarts(g):gapEnds(g);
            for c = 1:nCh
                tracks(k).x(c, gacombIdx) = interp1(borderIdx, tracks(k).x(c, borderIdx), gacombIdx);
                tracks(k).y(c, gacombIdx) = interp1(borderIdx, tracks(k).y(c, borderIdx), gacombIdx);
            end
        end
        tracks(k).gapStatus = gapStatus;
    end
    fprintf('\b\b\b\b%3d%%', round(100*k/(nTracks)));
end
fprintf('\n');

%====================================================================================
% Generate buffers before and after track, estimate gap values
%====================================================================================
%hasValidGaps = arrayfun(@(t) max([t.gapStatus 4])==4, tracks);
%singleTrackIdx = find([tracks.nSeg]==1 & hasValidGaps);
% singleTrackIdx = 1:numel(tracks);

% Gap map for fast indexing
%gapMap = zeros(numel(singleTrackIdx), data.movieLength);
gapMap = zeros(nTracks, data.movieLength);
for k = 1:nTracks
    gapMap(k, tracks(k).f(tracks(k).gapVect==1)) = 1;
end

% for buffers:
trackStarts = [tracks.start];
trackEnds = [tracks.end];
fullTracks = [tracks.visibility]==1;

fprintf('TrackAnalysis - Gap interpolation & generation of track buffers:     ');
for f = 1:data.movieLength
    
    mask = double(imread(data.maskPaths{f}));
    % binarize
    mask(mask~=0) = 1;
    labels = bwlabel(mask);
    
    for ch = 1:nCh
        frame = double(imread(data.framePaths{ch}{f}));
        
        %------------------------
        % Gaps
        %------------------------
        % tracks with valid gaps visible in current frame
        currentGapsIdx = find(gapMap(:,f));
        for ki = 1:numel(currentGapsIdx)
            k = currentGapsIdx(ki);
            
            % index in the track structure (.x etc)
            idxList = find(tracks(k).f==f & tracks(k).gapVect==1);
            
            for l = 1:numel(idxList)
                idx = idxList(l);
                xi = round(tracks(k).x(ch,idx));
                yi = round(tracks(k).y(ch,idx));
                
                % window/masks (see psfLocalization.m for details)
                maskWindow = labels(yi-w4:yi+w4, xi-w4:xi+w4);
                maskWindow(maskWindow==maskWindow(w4+1,w4+1)) = 0;
                
                cmask = annularMask;
                cmask(maskWindow~=0) = 0;
                window = frame(yi-w4:yi+w4, xi-w4:xi+w4);
                
                ci = mean(window(cmask==1));
                window(maskWindow~=0) = NaN;
                
                x0 = tracks(k).x(ch,idx)-xi;
                y0 = tracks(k).y(ch,idx)-yi;
                npx = sum(isfinite(window(:)));
                [prm, prmStd, ~, res] = fitGaussian2D(window, [x0 y0 max(window(:))-ci sigmaV(ch) ci], 'xyAc');
                dx = prm(1);
                dy = prm(2);
                if (dx > -w2 && dx < w2 && dy > -w2 && dy < w2)
                    tracks(k).x(ch,idx) = xi+dx;
                    tracks(k).y(ch,idx) = yi+dy;
                    tracks(k).A_pstd(ch,idx) = prmStd(3);
                    tracks(k).c_pstd(ch,idx) = prmStd(4);
                else
                    [prm, prmStd, ~, res] = fitGaussian2D(window, [x0 y0 max(window(:))-ci sigmaV(ch) ci], 'Ac');
                    tracks(k).A_pstd(ch,idx) = prmStd(1);
                    tracks(k).c_pstd(ch,idx) = prmStd(2);
                end
                tracks(k).A(ch,idx) = prm(3);
                tracks(k).c(ch,idx) = prm(5);
                
                tracks(k).sigma_r(ch,idx) = res.std;
                tracks(k).SE_sigma_r(ch,idx) = res.std/sqrt(2*(npx-1));
                
                SE_r = tracks(k).SE_sigma_r(ch,idx) * kLevel;
                
                tracks(k).pval_KS(idx) = res.pval;
                
                df2 = (npx-1) * (tracks(k).A_pstd(ch,idx).^2 + SE_r.^2).^2 ./...
                    (tracks(k).A_pstd(ch,idx).^4 + SE_r.^4);
                
                scomb = sqrt((tracks(k).A_pstd(ch,idx).^2 + SE_r.^2)/npx);
                T = (tracks(k).A(ch,idx) - res.std*kLevel) ./ scomb;
                tracks(k).pval_Ar(ch,idx) = 1-tcdf(T, df2);
            end            
        end
        
        %------------------------
        % start buffer
        %------------------------
        % tracks with start buffers in this frame
        cand = max(1, trackStarts-buffer)<=f & f<trackStarts;        
        % corresponding tracks, only if status = 1
        currentBufferIdx = find(cand & fullTracks);
     
        for ki = 1:length(currentBufferIdx)
            k = currentBufferIdx(ki);
        
            xi = round(tracks(k).x(ch,1));
            yi = round(tracks(k).y(ch,1));
            
            % window/masks (see psfLocalization.m for details)
            maskWindow = labels(yi-w4:yi+w4, xi-w4:xi+w4);
            maskWindow(maskWindow==maskWindow(w4+1,w4+1)) = 0;
            cmask = annularMask;
            cmask(maskWindow~=0) = 0;
            window = frame(yi-w4:yi+w4, xi-w4:xi+w4);
            ci = mean(window(cmask==1));
            window(maskWindow~=0) = NaN;
            
            x0 = tracks(k).x(1)-xi;
            y0 = tracks(k).y(1)-yi;
            
            [prm,prmStd,~,res] = fitGaussian2D(window, [x0 y0 max(window(:))-ci sigmaV(ch) ci], 'xyAc');
            bi = f - max(1, tracks(k).start-buffer) + 1;
            dx = prm(1);
            dy = prm(2);
            if sqrt((dx-x0)^2+(dy-y0)^2) < w2
                tracks(k).startBuffer.x(ch,bi) = xi+dx;
                tracks(k).startBuffer.y(ch,bi) = yi+dy;
                tracks(k).startBuffer.A_pstd(ch,bi) = prmStd(3);
            else
                [prm,prmStd,~,res] = fitGaussian2D(window, [x0 y0 max(window(:))-ci sigmaV(ch) ci], 'Ac');
                tracks(k).startBuffer.x(ch,bi) = tracks(k).x(ch,1);
                tracks(k).startBuffer.y(ch,bi) = tracks(k).y(ch,1);
                tracks(k).startBuffer.A_pstd(ch,bi) = prmStd(1);
            end
            tracks(k).startBuffer.A(ch,bi) = prm(3);
            tracks(k).startBuffer.c(ch,bi) = prm(5);
            tracks(k).startBuffer.sigma_r(ch,bi) = res.std;
            
            npx = sum(isfinite(window(:)));
            
            tracks(k).startBuffer.SE_sigma_r(ch,bi) = res.std/sqrt(2*(npx-1));
            SE_r = tracks(k).startBuffer.SE_sigma_r(ch,bi) * kLevel;
            
            df2 = (npx-1) * (tracks(k).startBuffer.A_pstd(ch,bi).^2 + SE_r.^2).^2 ./...
                (tracks(k).startBuffer.A_pstd(ch,bi).^4 + SE_r.^4);
            
            scomb = sqrt((tracks(k).startBuffer.A_pstd(ch,bi).^2 + SE_r.^2)/npx);
            T = (tracks(k).startBuffer.A(ch,bi) - res.std*kLevel) ./ scomb;
            tracks(k).startBuffer.pval_Ar(ch,bi) = 1-tcdf(T, df2);
        end
    
        %------------------------
        % end buffer
        %------------------------
        % segments with end buffers in this frame
        cand = trackEnds<f & f<=min(data.movieLength, trackEnds+buffer);
        % corresponding tracks
        currentBufferIdx = find(cand & fullTracks);
        
        for ki = 1:length(currentBufferIdx)
            k = currentBufferIdx(ki);
            xi = round(tracks(k).x(ch,end));
            yi = round(tracks(k).y(ch,end));
            
            % window/masks (see psfLocalization.m for details)
            maskWindow = labels(yi-w4:yi+w4, xi-w4:xi+w4);
            maskWindow(maskWindow==maskWindow(w4+1,w4+1)) = 0;
            cmask = annularMask;
            cmask(maskWindow~=0) = 0;
            window = frame(yi-w4:yi+w4, xi-w4:xi+w4);
            ci = mean(window(cmask==1));
            window(maskWindow~=0) = NaN;
            
            x0 = tracks(k).x(end)-xi;
            y0 = tracks(k).y(end)-yi;
            
            [prm,prmStd,~,res] = fitGaussian2D(window, [x0 y0 max(window(:))-ci sigmaV(ch) ci], 'xyAc');
            bi = f - tracks(k).end;
            dx = prm(1);
            dy = prm(2);
            if sqrt((dx-x0)^2+(dy-y0)^2) < w2
                tracks(k).endBuffer.x(ch,bi) = xi+dx;
                tracks(k).endBuffer.y(ch,bi) = yi+dy;
                tracks(k).endBuffer.A_pstd(ch,bi) = prmStd(3);
            else
                [prm,prmStd,~,res] = fitGaussian2D(window, [x0 y0 max(window(:))-ci sigmaV(ch) ci], 'Ac');
                tracks(k).endBuffer.x(ch,bi) = tracks(k).x(ch,end);
                tracks(k).endBuffer.y(ch,bi) = tracks(k).y(ch,end);
                tracks(k).endBuffer.A_pstd(ch,bi) = prmStd(1);
            end
            tracks(k).endBuffer.A(ch,bi) = prm(3);
            tracks(k).endBuffer.c(ch,bi) = prm(5);
            tracks(k).endBuffer.sigma_r(ch,bi) = res.std;
            
            npx = sum(isfinite(window(:)));
            
            tracks(k).endBuffer.SE_sigma_r(ch,bi) = res.std/sqrt(2*(npx-1));
            SE_r = tracks(k).endBuffer.SE_sigma_r(ch,bi) * kLevel;
            
            df2 = (npx-1) * (tracks(k).endBuffer.A_pstd(ch,bi).^2 + SE_r.^2).^2 ./...
                (tracks(k).endBuffer.A_pstd(ch,bi).^4 + SE_r.^4);
            
            scomb = sqrt((tracks(k).endBuffer.A_pstd(ch,bi).^2 + SE_r.^2)/npx);
            T = (tracks(k).endBuffer.A(ch,bi) - res.std*kLevel) ./ scomb;
            tracks(k).endBuffer.pval_Ar(ch,bi) = 1-tcdf(T, df2);
        end
        fprintf('\b\b\b\b%3d%%', round(100*(ch + (f-1)*nCh)/(nCh*data.movieLength)));
    end
end
fprintf('\n');
% 
% % sort tracks by type
% % idx = find([tracks.type]==1);
% % tracks = tracks([idx setdiff(1:nTracks, idx)]);
% 
% %==========================================
% % Compute displacements
% %==========================================
% % Only on tracks with no/valid gaps
% trackIdx = find(arrayfun(@(t) isempty([t.gapStatus{:}]) || max([t.gapStatus{:} 4])==4, tracks));
% fprintf('TrackAnalysis - Properties:     ');
% for ki = 1:length(trackIdx)
%     k = trackIdx(ki);
%     ns = numel(tracks(k).x);
%     msdVect = cell(1,ns);
%     msdStdVect = cell(1,ns);
%     for s = 1:ns
% 
%         x = tracks(k).x{s}(mCh,:);
%         y = tracks(k).y{s}(mCh,:);
%         tracks(k).totalDisplacement{s} = sqrt((x(end)-x(1))^2 + (y(end)-y(1))^2);        
%         % MSD
%         L = 10;
%         msdVect{s} = NaN(1,L);
%         msdStdVect{s} = NaN(1,L);
%         for l = 1:min(L, numel(x)-1)
%             tmp = (x(1+l:end)-x(1:end-l)).^2 + (y(1+l:end)-y(1:end-l)).^2;
%             msdVect{s}(l) = mean(tmp);
%             msdStdVect{s}(l) = std(tmp);
%         end
%         tracks(k).MSD = msdVect;
%         tracks(k).MSDstd = msdStdVect;
% 
%         %if L > 1 % min 2 points to fit
%         %    [D c alpha] = fitMSD(MSDvect(1:L), [MSDvect(L)/(4*L) 0 1], 'Dc');
%         %    tracks(k).D = D;
%         %    tracks(k).c = c;
%         %    tracks(k).alpha = alpha;
%         %end
%         
%         % add buffer time vectors
%         %b = size(tracks(k).startBuffer.x,2);
%         %tracks(k).startBuffer.t = ((-b:-1) + tracks(k).segmentStarts(s)-1) * data.framerate;
%         %b = size(tracks(k).endBuffer.x,2);
%         %tracks(k).endBuffer.t = (tracks(k).segmentEnds(s) + (1:b)-1) * data.framerate;
%     end
%     fprintf('\b\b\b\b%3d%%', round(100*k/(nTracks)));
% end
% fprintf('\n');

for k = 1:nTracks
    % add buffer time vectors
    if ~isempty(tracks(k).startBuffer)
        b = size(tracks(k).startBuffer.x,2);
        tracks(k).startBuffer.t = ((-b:-1) + tracks(k).start-1) * data.framerate;
    end
    if ~isempty(tracks(k).endBuffer)
        b = size(tracks(k).endBuffer.x,2);
        tracks(k).endBuffer.t = (tracks(k).end + (1:b)-1) * data.framerate;
    end
end

%==========================================
% Save results
%==========================================
if ~(exist([data.source 'Tracking'], 'dir')==7)
    mkdir([data.source 'Tracking']);
end
save([data.source 'Tracking' filesep filename], 'tracks');



% function [D c alpha] = fitMSD(MSD, prmVect, prmSel)
%
% opts = optimset('Jacobian', 'off', ...
%     'MaxFunEvals', 1e4, ...
%     'MaxIter', 1e4, ...
%     'Display', 'off', ...
%     'TolX', 1e-6, ...
%     'Tolfun', 1e-6);
%
%
% estIdx = false(1,3); % [D c alpha]
% estIdx(regexp('Dca', ['[' prmSel ']'])) = true;
%
% x = lsqnonlin(@costMSD, prmVect(estIdx), [], [], opts, MSD, prmVect, estIdx);
% prmVect(estIdx) = deal(abs(x));
%
% D = prmVect(1);
% c = prmVect(2); % constant offset
% alpha = prmVect(3);


% function [v] = costMSD(p, MSD, prmVect, estIdx)
% prmVect(estIdx) = deal(abs(p));
%
% D = prmVect(1);
% c = prmVect(2); % constant offset
% alpha = prmVect(3);
% %d % dimensionality
%
% t = 1:length(MSD);
%
% v = MSD - 4*D*t.^alpha - c;