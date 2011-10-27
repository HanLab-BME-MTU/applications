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

% Francois Aguet, November 2010 (last modified 09/21/2011)

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
tracks(1:nTracks) = struct('t', [],...
    'x', [], 'y', [], 'A', [], 'c', [],...
    'x_pstd', [], 'y_pstd', [], 'A_pstd', [], 'c_pstd', [],...
    'x_init', [], 'y_init', [], 'A_mask', [],...
    'sigma_r', [], 'SE_sigma_r', [],...
    'pval_Ar', [], 'pval_KS', [], 'isPSF', [],...
    'status', [], 'nSeg', [], 'gapStatus', [],...
    'gapStarts', [], 'gapEnds', [], 'gapLengths', [], 'gapVect', [],...
    'segmentStarts', [], 'segmentEnds', [],...
    'alphaMSD', [], 'MSD', [], 'MSDstd', [], 'totalDisplacement', [], 'D', [], ...
    'seqOfEvents', [], 'tracksFeatIndxCG', []);


% field names with multiple channels
mcFieldNames = {'x', 'y', 'A', 'c', 'x_pstd', 'y_pstd', 'A_pstd', 'c_pstd', 'sigma_r', 'SE_sigma_r', 'pval_Ar', 'pval_KS', 'isPSF'};
gapFieldNames = {'gapStatus', 'gapStarts', 'gapEnds', 'gapLengths', 'gapVect'};
bufferFieldNames = {'t', 'x', 'y', 'A', 'c', 'A_pstd', 'sigma_r'};

%==============================
% Loop through tracks
%==============================
fprintf('TrackAnalysis - Converting tracker output:     ');
for k = 1:nTracks
    if ~isstruct(trackinfo)
        x = trackinfo(k,1:8:end); % COM coordinates
        y = trackinfo(k,2:8:end);
        maskI = trackinfo(k,4:8:end);
        
        trackIdx = find(~isnan(x)); % index of detected track points
        
        firstIdx = trackIdx(1);
        lastIdx = trackIdx(end);
        
        tracks(k).xcom = x(firstIdx:lastIdx);
        tracks(k).ycom = y(firstIdx:lastIdx);
        tracks(k).maskI = maskI(firstIdx:lastIdx);
    else % new tracker
        % convert/assign structure fields
        tracks(k).seqOfEvents = trackinfo(k).seqOfEvents;
        tracks(k).tracksFeatIndxCG = trackinfo(k).tracksFeatIndxCG;
        
        firstIdx = trackinfo(k).seqOfEvents(1,1);
        lastIdx = trackinfo(k).seqOfEvents(end,1);
    end
    
    
    tracks(k).lifetime_s = (lastIdx-firstIdx+1)*data.framerate;
    
    tracks(k).start = firstIdx;
    tracks(k).end = lastIdx;

    %==============================================================================
    % Read amplitude & background from detectionResults.mat (localization results)
    %==============================================================================
    
    nSeg = size(tracks(k).tracksFeatIndxCG,1);
    tracks(k).nSeg = nSeg;
    
    % All segments are stored in cells, irrespective of whether they contain merges/splits
    
    for f = 1:length(mcFieldNames)
        tracks(k).(mcFieldNames{f}) = cell(1,nSeg);
    end
    for f = 1:length(gapFieldNames)
        tracks(k).(gapFieldNames{f}) = cell(1,nSeg);
    end
    tracks(k).segmentStarts = NaN(1,nSeg);
    tracks(k).segmentEnds = NaN(1,nSeg);
    for f = 1:length(bufferFieldNames)
        tracks(k).startBuffer.(bufferFieldNames{f}) = cell(1,nSeg);
        tracks(k).endBuffer.(bufferFieldNames{f}) = cell(1,nSeg);
    end
    
    for s = 1:nSeg
        ievents = tracks(k).seqOfEvents(tracks(k).seqOfEvents(:,3)==s, :);
        bounds = ievents(:,1);
        if ~isnan(ievents(2,4))
            bounds(2) = bounds(2)-1;
        end
        
        nf = bounds(2)-bounds(1)+1;
        for f = 1:length(mcFieldNames)
            tracks(k).(mcFieldNames{f}){s} = NaN(nCh,nf);
        end
        
        % for this segment
        %smStatus = ievents(:,4);
        frameRange = frameIdx(bounds(1):bounds(2)); % relative to movie (also when movie is subsampled)
        
        for i = 1:length(frameRange)
            idx = tracks(k).tracksFeatIndxCG(s, frameRange(i) - tracks(k).start + 1); % -> relative to IndxCG
            if idx ~= 0 % if not a gap, get detection values
                for f = 1:length(mcFieldNames)
                    tracks(k).(mcFieldNames{f}){s}(:,i) = frameInfo(frameRange(i)).(mcFieldNames{f})(:,idx);
                end
            end
        end
        tracks(k).segmentStarts(s) = bounds(1);
        tracks(k).segmentEnds(s) = bounds(2);
        tracks(k).t{s} = (bounds(1)-1:bounds(2)-1)*data.framerate;
    end
    fprintf('\b\b\b\b%3d%%', round(100*k/(nTracks)));
end
fprintf('\n');

% remove tracks that fall into image boundary
minx = arrayfun(@(t) min(cellfun(@(i) min(round(i(:))), t.x)), tracks);
maxx = arrayfun(@(t) max(cellfun(@(i) max(round(i(:))), t.x)), tracks);
miny = arrayfun(@(t) min(cellfun(@(i) min(round(i(:))), t.y)), tracks);
maxy = arrayfun(@(t) max(cellfun(@(i) max(round(i(:))), t.y)), tracks);
idx = minx<=w4 | miny<=w4 | maxx>nx-w4 | maxy>ny-w4;
tracks(idx) = [];
nTracks = length(tracks);


%=======================================
% Interpolate gaps and clean up tracks
%=======================================
fprintf('TrackAnalysis - Classification:     ');
for k = 1:nTracks
    
    nanIdx = cellfun(@(s) isnan(s(mCh,:)), tracks(k).x, 'UniformOutput', false);
    tracks(k).gapVect = nanIdx; % binary gap vector
    nanPoints = cellfun(@(s) sum(s), nanIdx);
    
    %=================================
    % Determine track and gap status
    %=================================
    if (1<tracks(k).start) && (tracks(k).end<nFrames) % complete tracks
        tracks(k).status = 1;
    elseif tracks(k).start==1 && tracks(k).end==nFrames % persistent tracks
        tracks(k).status = 3;
    else
        tracks(k).status = 2; % incomplete tracks
    end
    
    % loop through track segments with gaps
    segmentsWithGaps = find(nanPoints~=0);
    for s = segmentsWithGaps
        if nanIdx{s}(1)==1 || nanIdx{s}(end)==1 % temporary fix for segments that begin or end with gap
            tracks(k).gapStatus{s} = 5;
        else
            gacombIdx = diff(nanIdx{s});
            gapStarts = find(gacombIdx==1)+1;
            gapEnds = find(gacombIdx==-1);
            gapLengths = gapEnds-gapStarts+1;
            
            segmentIdx = diff([0 ~nanIdx{s} 0]); % these variables refer to segments between gaps
            segmentStarts = find(segmentIdx==1);
            segmentEnds = find(segmentIdx==-1)-1;
            segmentLengths = segmentEnds-segmentStarts+1;
            
            % loop over gaps
            nGaps = length(gapLengths);
            gv = 1:nGaps;
            gapStatus = 5*ones(1,nGaps);
            % gap valid if segments that precede/follow are > 1 frame or if gap is a single frame
            gapStatus(segmentLengths(gv)>1 & segmentLengths(gv+1)>1 | gapLengths(gv)==1) = 4;
            
            % fill position information for valid gaps using linear interpolation
            for g = 1:nGaps
                borderIdx = [gapStarts(g)-1 gapEnds(g)+1];
                gacombIdx = gapStarts(g):gapEnds(g);
                for c = 1:nCh
                    tracks(k).x{s}(c, gacombIdx) = interp1(borderIdx, tracks(k).x{s}(c, borderIdx), gacombIdx);
                    tracks(k).y{s}(c, gacombIdx) = interp1(borderIdx, tracks(k).y{s}(c, borderIdx), gacombIdx);
                end
            end
            tracks(k).gapStarts{s} = gapStarts;
            tracks(k).gapEnds{s} = gapEnds;
            tracks(k).gapLengths{s} = gapLengths;
            tracks(k).gapStatus{s} = gapStatus;
        end
    end
    
    if isempty(segmentsWithGaps)
        tracks(k).valid = tracks(k).status == 1;
    else
        % status for entire track (if a single gap is invalid, track is considered invalid)
        tracks(k).valid = tracks(k).status == 1 && max([tracks(k).gapStatus{:}]) == 4;
    end
    fprintf('\b\b\b\b%3d%%', round(100*k/(nTracks)));
end
fprintf('\n');

%====================================================================================
% Generate buffers before and after track, only for valid, single-segment tracks
%====================================================================================
% Loop through frames and fill gap and buffer values

hasValidGaps = arrayfun(@(t) ~isempty([t.gapStatus{:}]) && max([t.gapStatus{:} 4])==4, tracks);

[gapMap, segStarts, segEnds, seg2trackIndex, track2segIndex] = catTrackFields(tracks, data.movieLength, 'gapVect', mCh);
gapMap = gapMap==1;
status1tracks = find([tracks.status]==1);

fprintf('TrackAnalysis - Gap interpolation & generation of track buffers:     ');
for f = 1:data.movieLength
    
    mask = double(imread(data.maskPaths{f}));
    % binarize
    mask(mask~=0) = 1;
    labels = bwlabel(mask);
    
    for ch = 1:nCh
        frame = double(imread(data.framePaths{ch}{f}));
        
        % tracks with segments that contain gaps
        % -> tracks corresponding to segments in gapMap(:,f): trackSegIndex(gapMap(:,f))
        trackIdx = unique(seg2trackIndex(gapMap(:,f)));
        if ~isempty(trackIdx) % retain tracks with valid gaps only
            trackIdx = trackIdx(hasValidGaps(trackIdx));
        else
            trackIdx = [];
        end
        
        for ki = 1:length(trackIdx)
            k = trackIdx(ki);
            
            % segments with gaps in current frame
            segIdx = find(gapMap(seg2trackIndex==k,f))';
            
            % loop through track segments
            for si = 1:length(segIdx)
                s = segIdx(si);
                
                idx = f-tracks(k).segmentStarts(s) + 1;
                xi = round(tracks(k).x{s}(ch,idx));
                yi = round(tracks(k).y{s}(ch,idx));
                
                % window/masks (see psfLocalization.m for details)
                maskWindow = labels(yi-w4:yi+w4, xi-w4:xi+w4);
                maskWindow(maskWindow==maskWindow(w4+1,w4+1)) = 0;
                
                cmask = annularMask;
                cmask(maskWindow~=0) = 0;
                window = frame(yi-w4:yi+w4, xi-w4:xi+w4);
                
                ci = mean(window(cmask==1));
                window(maskWindow~=0) = NaN;
                
                x0 = tracks(k).x{s}(ch,idx)-xi;
                y0 = tracks(k).y{s}(ch,idx)-yi;
                npx = sum(isfinite(window(:)));
                [prm, prmStd, ~, res] = fitGaussian2D(window, [x0 y0 max(window(:))-ci sigmaV(ch) ci], 'xyAc');
                dx = prm(1);
                dy = prm(2);
                if (dx > -w2 && dx < w2 && dy > -w2 && dy < w2)
                    tracks(k).x{s}(ch,idx) = xi+dx;
                    tracks(k).y{s}(ch,idx) = yi+dy;
                    tracks(k).A_pstd{s}(ch,idx) = prmStd(3);
                    tracks(k).c_pstd{s}(ch,idx) = prmStd(4);
                else
                    [prm, prmStd, ~, res] = fitGaussian2D(window, [x0 y0 max(window(:))-ci sigmaV(ch) ci], 'Ac');
                    tracks(k).A_pstd{s}(ch,idx) = prmStd(1);
                    tracks(k).c_pstd{s}(ch,idx) = prmStd(2);
                end
                tracks(k).A{s}(ch,idx) = prm(3);
                tracks(k).c{s}(ch,idx) = prm(5);
                
                tracks(k).sigma_r{s}(ch,idx) = res.std;
                tracks(k).SE_sigma_r{s}(ch,idx) = res.std/sqrt(2*(npx-1));
                
                SE_r = tracks(k).SE_sigma_r{s}(ch,idx) * kLevel;
                
                tracks(k).pval_KS{s}(idx) = res.pval;
                
                df2 = (npx-1) * (tracks(k).A_pstd{s}(ch,idx).^2 + SE_r.^2).^2 ./...
                    (tracks(k).A_pstd{s}(ch,idx).^4 + SE_r.^4);
                scomb = sqrt((tracks(k).A_pstd{s}(ch,idx).^2 + SE_r.^2)/npx);
                T = (tracks(k).A{s}(ch,idx) - res.std*kLevel) ./ scomb;
                tracks(k).pval_Ar{s}(ch,idx) = tcdf(T, df2);
            end
        end
        
        %------------------------
        % start buffer
        %------------------------
        % segments with start buffers in this frame
        segCand = max(1, segStarts-buffer)<=f & f<segStarts;
        % corresponding tracks, only if status = 1
        trackIdx = intersect(status1tracks, unique(seg2trackIndex(segCand)));
        
        for ki = 1:length(trackIdx)
            k = trackIdx(ki);
            
            segIdx = find(segCand(seg2trackIndex==k))';
            for si = 1:length(segIdx)
                s = segIdx(si);
                
                xi = round(tracks(k).x{s}(ch,1));
                yi = round(tracks(k).y{s}(ch,1));
            
                % window/masks (see psfLocalization.m for details)
                maskWindow = labels(yi-w4:yi+w4, xi-w4:xi+w4);
                maskWindow(maskWindow==maskWindow(w4+1,w4+1)) = 0;
                cmask = annularMask;
                cmask(maskWindow~=0) = 0;
                window = frame(yi-w4:yi+w4, xi-w4:xi+w4);
                ci = mean(window(cmask==1));
                window(maskWindow~=0) = NaN;
                
                x0 = tracks(k).x{s}(1)-xi;
                y0 = tracks(k).y{s}(1)-yi;

                [prm,prmStd,~,res] = fitGaussian2D(window, [x0 y0 max(window(:))-ci sigmaV(ch) ci], 'xyAc');
                bi = f - max(1, tracks(k).segmentStarts(s)-buffer) + 1;
                dx = prm(1);
                dy = prm(2);
                if sqrt((dx-x0)^2+(dy-y0)^2) < w2
                    tracks(k).startBuffer.x{s}(ch,bi) = xi+dx;
                    tracks(k).startBuffer.y{s}(ch,bi) = yi+dy;
                    tracks(k).startBuffer.A_pstd{s}(ch,bi) = prmStd(3);
                else
                    [prm,prmStd,~,res] = fitGaussian2D(window, [x0 y0 max(window(:))-ci sigmaV(ch) ci], 'Ac');
                    tracks(k).startBuffer.x{s}(ch,bi) = tracks(k).x{s}(ch,1);
                    tracks(k).startBuffer.y{s}(ch,bi) = tracks(k).y{s}(ch,1);
                    tracks(k).startBuffer.A_pstd{s}(ch,bi) = prmStd(1);
                end
                tracks(k).startBuffer.A{s}(ch,bi) = prm(3);
                tracks(k).startBuffer.c{s}(ch,bi) = prm(5);
                tracks(k).startBuffer.sigma_r{s}(ch,bi) = res.std;
            end
        end
        
        %------------------------
        % end buffer
        %------------------------
        % segments with end buffers in this frame
        segCand = segEnds<f & f<=min(data.movieLength, segEnds+buffer);
        % corresponding tracks
        trackIdx = intersect(status1tracks, unique(seg2trackIndex(segCand)));
        
        for ki = 1:length(trackIdx)
            k = trackIdx(ki);
            segIdx = find(segCand(seg2trackIndex==k))';
            for si = 1:length(segIdx)
                s = segIdx(si);
                
                xi = round(tracks(k).x{s}(ch,end));
                yi = round(tracks(k).y{s}(ch,end));
                
                % window/masks (see psfLocalization.m for details)
                maskWindow = labels(yi-w4:yi+w4, xi-w4:xi+w4);
                maskWindow(maskWindow==maskWindow(w4+1,w4+1)) = 0;
                cmask = annularMask;
                cmask(maskWindow~=0) = 0;
                window = frame(yi-w4:yi+w4, xi-w4:xi+w4);
                ci = mean(window(cmask==1));
                window(maskWindow~=0) = NaN;
                
                x0 = tracks(k).x{s}(end)-xi;
                y0 = tracks(k).y{s}(end)-yi;
                
                [prm,prmStd,~,res] = fitGaussian2D(window, [x0 y0 max(window(:))-ci sigmaV(ch) ci], 'xyAc');
                bi = f - tracks(k).segmentEnds(s);
                dx = prm(1);
                dy = prm(2);
                if sqrt((dx-x0)^2+(dy-y0)^2) < w2
                    tracks(k).endBuffer.x{s}(ch,bi) = xi+dx;
                    tracks(k).endBuffer.y{s}(ch,bi) = yi+dy;
                    tracks(k).endBuffer.A_pstd{s}(ch,bi) = prmStd(3);
                else
                    [prm,prmStd,~,res] = fitGaussian2D(window, [x0 y0 max(window(:))-ci sigmaV(ch) ci], 'Ac');
                    tracks(k).endBuffer.x{s}(ch,bi) = tracks(k).x{s}(ch,end);
                    tracks(k).endBuffer.y{s}(ch,bi) = tracks(k).y{s}(ch,end);
                    tracks(k).endBuffer.A_pstd{s}(ch,bi) = prmStd(1);
                end
                tracks(k).endBuffer.A{s}(ch,bi) = prm(3);
                tracks(k).endBuffer.c{s}(ch,bi) = prm(5);
                tracks(k).endBuffer.sigma_r{s}(ch,bi) = res.std;
            end
        end
        fprintf('\b\b\b\b%3d%%', round(100*(ch + (f-1)*nCh)/(nCh*data.movieLength)));
    end
end
fprintf('\n');

% sort tracks by type
% idx = find([tracks.type]==1);
% tracks = tracks([idx setdiff(1:nTracks, idx)]);

%==========================================
% Compute displacements
%==========================================
% Only on tracks with no/valid gaps
trackIdx = find(arrayfun(@(t) isempty([t.gapStatus{:}]) || max([t.gapStatus{:} 4])==4, tracks));
fprintf('TrackAnalysis - Properties:     ');
for ki = 1:length(trackIdx)
    k = trackIdx(ki);
    ns = numel(tracks(k).x);
    msdVect = cell(1,ns);
    msdStdVect = cell(1,ns);
    for s = 1:ns

        x = tracks(k).x{s}(mCh,:);
        y = tracks(k).y{s}(mCh,:);
        tracks(k).totalDisplacement{s} = sqrt((x(end)-x(1))^2 + (y(end)-y(1))^2);        
        % MSD
        L = 10;
        msdVect{s} = NaN(1,L);
        msdStdVect{s} = NaN(1,L);
        for l = 1:min(L, numel(x)-1)
            tmp = (x(1+l:end)-x(1:end-l)).^2 + (y(1+l:end)-y(1:end-l)).^2;
            msdVect{s}(l) = mean(tmp);
            msdStdVect{s}(l) = std(tmp);
        end
        tracks(k).MSD = msdVect;
        tracks(k).MSDstd = msdStdVect;

        %if L > 1 % min 2 points to fit
        %    [D c alpha] = fitMSD(MSDvect(1:L), [MSDvect(L)/(4*L) 0 1], 'Dc');
        %    tracks(k).D = D;
        %    tracks(k).c = c;
        %    tracks(k).alpha = alpha;
        %end
        
        % add buffer time vectors
        %b = size(tracks(k).startBuffer.x{s},2);
        %tracks(k).startBuffer.t{s} = ((-b:-1) + tracks(k).segmentStarts(s)-1) * data.framerate;
        %b = size(tracks(k).endBuffer.x{s},2);
        %tracks(k).endBuffer.t{s} = (tracks(k).segmentEnds(s) + (1:b)-1) * data.framerate;
    end
    fprintf('\b\b\b\b%3d%%', round(100*k/(nTracks)));
end
fprintf('\n');

for k = 1:nTracks
    ns = numel(tracks(k).x);
    for s = 1:ns
        % add buffer time vectors
        b = size(tracks(k).startBuffer.x{s},2);
        tracks(k).startBuffer.t{s} = ((-b:-1) + tracks(k).segmentStarts(s)-1) * data.framerate;
        b = size(tracks(k).endBuffer.x{s},2);
        tracks(k).endBuffer.t{s} = (tracks(k).segmentEnds(s) + (1:b)-1) * data.framerate;
    end
end

trackInfo.x = catTrackFields(tracks, data.movieLength, 'x', mCh);
trackInfo.y = catTrackFields(tracks, data.movieLength, 'y', mCh);
trackInfo.gapMap = gapMap;
trackInfo.segStarts = segStarts;
trackInfo.segEnds = segEnds;
trackInfo.seg2trackIndex = seg2trackIndex;
trackInfo.track2segIndex = track2segIndex;
trackInfo.nSeg = [tracks.nSeg];
trackInfo.status = [tracks.status];
trackInfo.valid = [tracks.valid];

%==========================================
% Save results
%==========================================
if ~(exist([data.source 'Tracking'], 'dir')==7)
    mkdir([data.source 'Tracking']);
end
save([data.source 'Tracking' filesep filename], 'tracks', 'trackInfo');


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