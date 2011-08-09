%[data] = runTrackAnalysis(data, varargin)
%
% INPUTS    data        : array of experiment structures
%           {Buffer}    : length of buffer before and after tracks
%           {Overwrite} : 1 to overwrite previous results
%           {FileName}  : name of the output
%
% Usage example: runTrackAnalysis(data, 'Buffer', 10);
%
% Note: Only tracks with track.status==1 and track.gapStatus==4 are considered.

% Francois Aguet, November 2010 (last modified 08/08/2011)

function [data] = runTrackAnalysis(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParamValue('Buffer', 5, @isscalar);
ip.addParamValue('BufferMode', 'xyAc', @(x) strcmpi(x, 'xyAc') | strcmpi(x, 'Ac'));
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
bufferMode = ip.Results.BufferMode;

for i = 1:length(data)
    data(i).tracks = [];
    data(i).smTracks = [];
end
parfor i = 1:length(data)
    if ~(exist([data(i).source filesep 'Tracking' filesep filename],'file')==2) || overwrite
        data(i) = main(data(i), buffer, trackerOutput, filename, frameIdx{i}, bufferMode);
    else
        fprintf('TrackAnalysis has already been run for: %s\n', getShortPath(data(i)));
    end
end



function [data] = main(data, buffer, trackerOutput, filename, frameIdx, bufferMode)

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
    'status', [], 'gapStatus', [],...
    'gapStarts', [], 'gapEnds', [], 'gapLengths', [], 'gapVect', [],...
    'segmentStarts', [], 'segmentEnds', [],...
    'alphaMSD', [], 'MSD', [], 'MSDstd', [], 'totalDisplacement', [], 'D', [], ...
    'seqOfEvents', [], 'tracksFeatIndxCG', []);


% field names with multiple channels
mcFieldNames = {'x', 'y', 'A', 'c', 'x_pstd', 'y_pstd', 'A_pstd', 'c_pstd', 'sigma_r', 'SE_sigma_r', 'pval_Ar', 'pval_KS', 'isPSF'};
gapFieldNames = {'gapStatus', 'gapStarts', 'gapEnds', 'gapLengths', 'gapVect', 'segmentStarts', 'segmentEnds'};

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
%     if nSeg==1
%         tracks(k).t = (firstIdx-1:lastIdx-1)*data.framerate;
%         nf = length(tracks(k).t);
%         
%         % initialize fields
%         for f = 1:length(mcFieldNames)
%             tracks(k).(mcFieldNames{f}) = NaN(nCh,nf);
%         end
%         tracks(k).pval_KS = NaN(1,nf);
%         tracks(k).isPSF = NaN(1,nf);
%         
%         % index of frames in case movie is subsampled
%         frameRange = frameIdx(tracks(k).start:tracks(k).end);
%         % copy data from frameInfo (detection) to track
%         for i = 1:length(frameRange)
%             idx = tracks(k).tracksFeatIndxCG(i); % feature index in frame
%             %idx = trackedFeatureNum(k, frameRange(i)); % for old tracker
%             if idx ~= 0 % if not a gap
%                 for f = 1:length(mcFieldNames)
%                     tracks(k).(mcFieldNames{f})(:,i) = frameInfo(frameRange(i)).(mcFieldNames{f})(:,idx);
%                 end
%                 tracks(k).pval_KS(i) = frameInfo(frameRange(i)).pval_KS(idx);
%                 tracks(k).isPSF(i) = frameInfo(frameRange(i)).isPSF(idx);
%             end
%         end
    %else % contains merges/splits -> segments stored in cell array
        tracks(k).t = (firstIdx-1:lastIdx-1)*data.framerate;
        for f = 1:length(mcFieldNames)
            tracks(k).(mcFieldNames{f}) = cell(1,nSeg);
        end
        for f = 1:length(gapFieldNames)
            tracks(k).(gapFieldNames{f}) = cell(1,nSeg);
        end

        for s = 1:nSeg
            bounds = tracks(k).seqOfEvents(tracks(k).seqOfEvents(:,3)==s);
            nf = bounds(2)-bounds(1)+1;
            for f = 1:length(mcFieldNames)
                tracks(k).(mcFieldNames{f}){s} = NaN(nCh,nf);
            end
            
            % for this segment
            smStatus = tracks(k).seqOfEvents(tracks(k).seqOfEvents(:,3)==s,4);
            if ~isnan(smStatus(2)) % split or merge
                frameRange = frameIdx(bounds(1):bounds(2)-1);
            else
                frameRange = frameIdx(bounds(1):bounds(2)); % relative to movie (also when movie is subsampled)
            end
            for i = 1:length(frameRange)
                idx = tracks(k).tracksFeatIndxCG(s, frameRange(i) - tracks(k).start + 1); % -> relative to IndxCG
                if idx ~= 0 % if not a gap, get detection values
                    for f = 1:length(mcFieldNames)
                        tracks(k).(mcFieldNames{f}){s}(:,i) = frameInfo(frameRange(i)).(mcFieldNames{f})(:,idx);
                    end
                end
            end
            tracks(k).segmentStarts{s} = bounds(1);
            tracks(k).segmentEnds{s} = bounds(2);
        end
    %end
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
    %x = tracks(k).x{s};
    %y = tracks(k).y{s};
    %nanIdx = isnan(x(mCh,:));
    
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
        %gacombIdx = diff([0 nanIdx{s} 0]);
        %gapStarts = find(gacombIdx==1);
        if nanIdx{s}(1)==1 || nanIdx{s}(end)==1 % temporary fix for segments that begin or end with gap
            tracks(k).gapStatus{s} = 5;
        else
            gacombIdx = diff(nanIdx{s});
            gapStarts = find(gacombIdx==1)+1;
            gapEnds = find(gacombIdx==-1);
            gapLengths = gapEnds-gapStarts+1;
            
            segmentIdx = diff([0 ~nanIdx{s} 0]);
            segmentStarts = find(segmentIdx==1);
            segmentEnds = find(segmentIdx==-1)-1;
            segmentLengths = segmentEnds-segmentStarts+1;
            
            % loop over gaps
            nGaps = length(gapLengths);
            gv = 1:nGaps;
            gapStatus = 5*ones(1,nGaps);
            % gap is valid if segments that precede & follow are > 1 frame
            % of if gap is a single frame
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
    
    
    %==========================================
    % Compute displacements
    %==========================================    
    if tracks(k).valid && numel(tracks(k).x)==1

        tracks(k).totalDisplacement = sqrt((x(end)-x(1))^2 + (y(end)-y(1))^2);        
        % MSD
        L = min(5, length(tracks(k).t));
        msdVect = zeros(1,L);
        msdStdVect = zeros(1,L);
        for l = 1:L
            tmp = sqrt((x(1+l:end)-x(1:end-l)).^2 + (y(1+l:end)-y(1:end-l)).^2);
            msdVect(l) = mean(tmp);
            msdStdVect(l) = std(tmp);
            
        end
        tracks(k).MSD = msdVect;
        tracks(k).MSDstd = msdStdVect;

        %if L > 1 % min 2 points to fit
        %    [D c alpha] = fitMSD(MSDvect(1:L), [MSDvect(L)/(4*L) 0 1], 'Dc');
        %    tracks(k).D = D;
        %    tracks(k).c = c;
        %    tracks(k).alpha = alpha;
        %end
    end
    fprintf('\b\b\b\b%3d%%', round(100*k/(nTracks)));
end
fprintf('\n');

% Re-assign fields for valid regular tracks
rTrackIdx = find(arrayfun(@(t) numel(t.x)==1, tracks));
for k = rTrackIdx
    for f = 1:length(mcFieldNames)
        tracks(k).(mcFieldNames{f}) = tracks(k).(mcFieldNames{f}){1};
    end
    for f = 1:length(gapFieldNames)
        tracks(k).(gapFieldNames{f}) = tracks(k).(gapFieldNames{f}){1};
    end
end
nRegTracks = length(rTrackIdx);
tracks = tracks([rTrackIdx setdiff(1:nTracks, rTrackIdx)]);
rTrackIdx = 1:nRegTracks;

%====================================================================================
% Generate buffers before and after track, only for valid, single-segment tracks
%====================================================================================
% Loop through frames and fill gap and buffer values

% number of buffer frames before and after track
bStart = [tracks(rTrackIdx).start] - max(1, [tracks(rTrackIdx).start]-buffer);
bEnd = min(data.movieLength, [tracks(rTrackIdx).end]+buffer) - [tracks(rTrackIdx).end];

gapMap = catTrackFields(tracks(rTrackIdx), data.movieLength, 'gapVect', 1)==1;

fprintf('TrackAnalysis - Gap interpolation & generation of track buffers:     ');
for f = 1:data.movieLength
    
    mask = double(imread(data.maskPaths{f}));
    % binarize
    mask(mask~=0) = 1;
    labels = bwlabel(mask);
    
    for ch = 1:nCh
        frame = double(imread(data.framePaths{ch}{f}));
        
        % tracks with gaps in this frame
        trackIdx = find(gapMap(:,f) & [tracks(rTrackIdx).valid]')';
        for k = trackIdx
            
            idx = f-tracks(k).start + 1;
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
            
            npx = sum(isfinite(window(:)));
            [prm, prmStd, ~, res] = fitGaussian2D(window, [tracks(k).x(ch,idx)-xi tracks(k).y(ch,idx)-yi max(window(:))-ci sigmaV(ch) ci], 'Ac');
            
            tracks(k).A(ch,idx) = prm(3);
            tracks(k).c(ch,idx) = prm(5);
            
            tracks(k).A_pstd(ch,idx) = prmStd(1);
            tracks(k).c_pstd(ch,idx) = prmStd(2);

            tracks(k).sigma_r(ch,idx) = res.std;
            tracks(k).SE_sigma_r(ch,idx) = res.std/sqrt(2*(npx-1));

            SE_r = tracks(k).SE_sigma_r(ch,idx) * kLevel;
                    
            tracks(k).pval_KS(idx) = res.pval;
                    
            df2 = (npx-1) * (tracks(k).A_pstd(ch,idx).^2 + SE_r.^2).^2 ./...
                (tracks(k).A_pstd(ch,idx).^4 + SE_r.^4);
            scomb = sqrt((tracks(k).A_pstd(ch,idx).^2 + SE_r.^2)/npx);
            T = (tracks(k).A(ch,idx) - res.std*kLevel) ./ scomb;
            tracks(k).pval_Ar(ch,idx) = tcdf(T, df2);
        end
        
        % start buffers in this frame
        %trackIdx = find([tracks.start]-bStart<=f & f<[tracks.start] & [tracks.valid]);
        trackIdx = find([tracks(rTrackIdx).start]-bStart<=f & f<[tracks(rTrackIdx).start] & [tracks(rTrackIdx).valid]);
        %trackIdx = find(trackMap(:,f) & [tracks.valid]')
        for k = trackIdx
            bi = f - tracks(k).start + bStart(k) + 1;
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
            
            x_init = tracks(k).x(1)-xi;
            y_init = tracks(k).y(1)-yi;

            [prm,prmStd,~,res] = fitGaussian2D(window, [x_init y_init max(window(:))-ci sigmaV(ch) ci], bufferMode);
            
            if sqrt((prm(1)-x_init)^2+(prm(2)-y_init)^2) < (tracks(k).MSD(1) + tracks(k).MSDstd(1))
                tracks(k).startBuffer.x(ch,bi) = xi+prm(1);
                tracks(k).startBuffer.y(ch,bi) = yi+prm(2);
            else
                [prm,prmStd,~,res] = fitGaussian2D(window, [x_init y_init max(window(:))-ci sigmaV(ch) ci], 'Ac');
                tracks(k).startBuffer.x(ch,bi) = tracks(k).x(ch,1);
                tracks(k).startBuffer.y(ch,bi) = tracks(k).y(ch,1);
            end
            tracks(k).startBuffer.A(ch,bi) = prm(3);
            tracks(k).startBuffer.c(ch,bi) = prm(5);
            tracks(k).startBuffer.A_pstd(ch,bi) = prmStd(1);
            tracks(k).startBuffer.sigma_r(ch,bi) = res.std;
        end
        
        % end buffers in this frame
        trackIdx = find([tracks(rTrackIdx).end]<f & f<=([tracks(rTrackIdx).end]+bEnd) & [tracks(rTrackIdx).valid]);
        for k = trackIdx
            bi = f-tracks(k).end;
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
            
            x_init = tracks(k).x(end)-xi;
            y_init = tracks(k).y(end)-yi;
            
            [prm,prmStd,~,res] = fitGaussian2D(window, [x_init y_init max(window(:))-ci sigmaV(ch) ci], bufferMode);
            
            if sqrt((prm(1)-x_init)^2+(prm(2)-y_init)^2) < (tracks(k).MSD(1) + tracks(k).MSDstd(1))
                tracks(k).endBuffer.x(ch,bi) = xi+prm(1);
                tracks(k).endBuffer.y(ch,bi) = yi+prm(2);                
            else
                [prm,prmStd,~,res] = fitGaussian2D(window, [x_init y_init max(window(:))-ci sigmaV(ch) ci], 'Ac');
                tracks(k).endBuffer.x(ch,bi) = tracks(k).x(ch,end);
                tracks(k).endBuffer.y(ch,bi) = tracks(k).y(ch,end);
            end
            tracks(k).endBuffer.A(ch,bi) = prm(3);
            tracks(k).endBuffer.c(ch,bi) = prm(5);
            tracks(k).endBuffer.A_pstd(ch,bi) = prmStd(1);
            tracks(k).endBuffer.sigma_r(ch,bi) = res.std;
        end
        fprintf('\b\b\b\b%3d%%', round(100*(ch + (f-1)*nCh)/(nCh*data.movieLength)));
    end
end
fprintf('\n');


%==========================================
% Save results
%==========================================
data.tracks = tracks;

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