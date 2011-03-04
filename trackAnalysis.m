%[tracks nMergeSplit] = trackAnalysis(data, buffer, filename)
%
% INPUTS:     data       : experiment structure from 'loadConditionData()'
%             {buffer}   :
%             {filename} : name of the tracker file to load
%
% The only tracks of interest are those with status 1 or 4.

% Francois Aguet, June 2010 (revised from 'determineLifetimeStatus.m')

function [tracks nMergeSplit] = trackAnalysis(data, buffer, filename)

if nargin<2 || isempty(buffer)
    buffer = 5;
end
if nargin<3 || isempty(filename)
    filename = 'trackedFeatures.mat';
end

ny = data.imagesize(1);
nx = data.imagesize(2);
nFrames = data.movieLength;
nCh = length(data.channels);

frameList = cell(1,nCh);
for k = 1:nCh
    frameList{k} = dir([data.channels{k} '*.tif*']);
end
maskList = dir([data.source 'Detection' filesep 'Masks' filesep '*.tif']);

load([data.source 'Detection' filesep 'detectionResults.mat']);


%=================================
% Determine master/slave channels
%=================================
masterChannel = regexp(data.source, data.channels);
masterChannel = find([masterChannel{:}]);
slaveChannels = setdiff(1:nCh, masterChannel);

sigma = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, data.markers{1});
w = ceil(3*sigma);

%=================================
% Read and convert tracker output
%=================================
tPath = [data.source 'Tracking' filesep filename];
if exist(tPath, 'file')==2
    trackinfo = load(tPath);
    trackinfo = trackinfo.tracksFinal;
    
    % Filter out tracks with merge/split events
    msIdx = arrayfun(@(x) size(x.seqOfEvents, 1)>2, trackinfo);
    nMergeSplit = length(msIdx);
    trackinfo(msIdx) = [];
    
    nTracks = length(trackinfo);
    
elseif exist([data.source 'TrackInfoMatrices' filesep 'trackedFeatures.mat'], 'file')==2
    % (for old tracker. oldest version: trackInfo.mat)
    trackinfo = load([data.source 'TrackInfoMatrices' filesep 'trackedFeatures.mat']);
    trackedFeatureNum = trackinfo.trackedFeatureNum;
    trackinfo = trackinfo.trackedFeatureInfo;
    nTracks = size(trackinfo, 1);
    nMergeSplit = NaN;
else
    error('No valid tracker output found.');
end
tracks(1:nTracks) = struct('t', [],...
    'x', [], 'x_std', [],...
    'y', [], 'y_std', [],...
    'A', [], 'A_std', [],...
    'c', [], 'cStd', [], 'cStd_mask', [],...
    'maskI', [],...
    'status', [], 'gapStatus', [],...
    'gapStarts', [], 'gapEnds', [], 'gapLengths', [],...
    'segmentStarts', [], 'segmentEnds', [], 'segmentLengths', [],...
    'D', [], 'alpha', [], 'MSDvect', [], 'totalD', [],...
    'seqOfEvents', [], 'tracksFeatIndxCG', []);


%==============================
% Loop through tracks
%==============================
fprintf('TrackAnalysis - Converting tracker output:     ');
for k = 1:nTracks
    if ~isstruct(trackinfo)
        x = trackinfo(k,1:8:end); % COM coordinates
        y = trackinfo(k,2:8:end);
        maskI = trackinfo(k,4:8:end);
        
        % index of detected track points
        trackIdx = find(~isnan(x));
        
        firstIdx = trackIdx(1);
        lastIdx = trackIdx(end);
        trackLength = lastIdx-firstIdx+1;
        %trackPoints = length(trackIdx);
        
        tracks(k).xcom = x(firstIdx:lastIdx);
        tracks(k).ycom = y(firstIdx:lastIdx);
        tracks(k).maskI = maskI(firstIdx:lastIdx);
        
        tracks(k).t = (firstIdx-1:lastIdx-1)*data.framerate;
        tracks(k).lifetime_s = trackLength*data.framerate;
        
        tracks(k).start = firstIdx;
        tracks(k).end = lastIdx;
    else
        % convert/assign structure fields
        tracks(k).xcom = trackinfo(k).tracksCoordAmpCG(1:8:end);
        tracks(k).ycom = trackinfo(k).tracksCoordAmpCG(2:8:end);
        tracks(k).maskI = trackinfo(k).tracksCoordAmpCG(4:8:end);
        tracks(k).seqOfEvents = trackinfo(k).seqOfEvents;
        tracks(k).tracksFeatIndxCG = trackinfo(k).tracksFeatIndxCG;
        
        tracks(k).lifetime_s = length(tracks(k).x)*data.framerate;
        tracks(k).t = (trackinfo(k).seqOfEvents(1,1)-1:trackinfo(k).seqOfEvents(2,1)-1)*data.framerate;
        
        firstIdx = trackinfo(k).seqOfEvents(1,1);
        lastIdx = trackinfo(k).seqOfEvents(end,1);
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%
        %trackLength = lastIdx-firstIdx+1;
        %trackPoints = sum(~isnan(x));
        
        tracks(k).start = firstIdx;
        tracks(k).end = lastIdx;
        
    end
    
    %==============================================================================
    % Read amplitude & background from detectionResults.mat (localization results)
    %==============================================================================
    nf = length(tracks(k).xcom);
    tracks(k).x = NaN(1,nf);
    tracks(k).x_std = NaN(1,nf);
    tracks(k).y = NaN(1,nf);
    tracks(k).y_std = NaN(1,nf);
    tracks(k).A = NaN(nCh, nf); % NaNs indicate gaps
    tracks(k).A_std = NaN(nCh,nf);
    tracks(k).c = NaN(nCh, nf);
    tracks(k).cStd = NaN(nCh, nf);
    tracks(k).cStd_mask = NaN(nCh, nf);
    
    frameRange = tracks(k).start:tracks(k).end;
    
    for i = 1:length(frameRange)
        if isstruct(trackinfo)
            idx = tracks(k).tracksFeatIndxCG(i);
        else
            idx = trackedFeatureNum(k, frameRange(i)); % for old tracker
        end
        if idx ~= 0
            tracks(k).x(i) = frameInfo(frameRange(i)).xloc(idx);
            tracks(k).x_std(i) = frameInfo(frameRange(i)).x_std(:,idx);
            tracks(k).y(i) = frameInfo(frameRange(i)).yloc(idx);
            tracks(k).y_std(i) = frameInfo(frameRange(i)).y_std(:,idx);
            tracks(k).A(:,i) = frameInfo(frameRange(i)).A(:,idx);
            tracks(k).A_std(:,i) = frameInfo(frameRange(i)).A_std(:,idx);
            tracks(k).c(:,i) = frameInfo(frameRange(i)).c(:,idx);
            tracks(k).cStd(:,i) = frameInfo(frameRange(i)).cStd(:,idx);
            tracks(k).cStd_mask(:,i) = frameInfo(frameRange(i)).cStd_mask(:,idx);
            
            tracks(k).pval_KS(i) = frameInfo(frameRange(i)).pval_KS(idx);
            tracks(k).valid_KS(i) = frameInfo(frameRange(i)).valid_KS(idx);
            tracks(k).maskI(i) = frameInfo(frameRange(i)).amp(idx);
        end
    end
    fprintf('\b\b\b\b%3d%%', round(100*k/(nTracks)));
end
fprintf('\n');


% Remove all tracks that entirely consist of NaNs or that are at image border
trackLengths = [tracks.end] - [tracks.start] + 1;
nanCount = arrayfun(@(t) sum(isnan(t.x)), tracks);

minx = arrayfun(@(t) min(round(t.xcom)), tracks);
maxx = arrayfun(@(t) max(round(t.xcom)), tracks);
miny = arrayfun(@(t) min(round(t.ycom)), tracks);
maxy = arrayfun(@(t) max(round(t.ycom)), tracks);

idx = trackLengths==nanCount | minx<=w | miny<=w | maxx>nx-w | maxy>ny-w;
tracks(idx) = [];

nTracks = length(tracks);


%=======================================
% Interpolate gaps and clean up tracks
%=======================================
fprintf('TrackAnalysis - Classification and gap interpolation:     ');
for k = 1:nTracks
    
    % prune track if 'gap' at beginning or end (NaNs from localization)
    nanIdx = isnan(tracks(k).x);
    diffIdx = diff(nanIdx);
    rmIdx = [];
    
    trackLength = tracks(k).end-tracks(k).start+1;
    
    if nanIdx(1)
        endNaN = find(diffIdx==-1, 1, 'first');
        rmIdx = 1:endNaN;
        tracks(k).start = tracks(k).start + endNaN;
    end
    if nanIdx(end)
        begNaN = find(diffIdx==1, 1, 'last')+1;
        rmIdx = [rmIdx begNaN:trackLength];
        tracks(k).end = tracks(k).end - (trackLength-begNaN+1);
    end
        
    if ~isempty(rmIdx)
        trackLength = trackLength - length(rmIdx);
        tracks(k).t(rmIdx) = [];
        tracks(k).x(rmIdx) = [];
        tracks(k).x_std(rmIdx) = [];
        tracks(k).y(rmIdx) = [];
        tracks(k).y_std(rmIdx) = [];
        tracks(k).A(:,rmIdx) = [];
        tracks(k).A_std(:,rmIdx) = [];
        tracks(k).c(:,rmIdx) = [];
        tracks(k).cStd(:,rmIdx) = [];
        tracks(k).cStd_mask(:,rmIdx) = [];
        
        tracks(k).xcom(rmIdx) = [];
        tracks(k).ycom(rmIdx) = [];
        
        tracks(k).pval_KS(rmIdx) = [];
        tracks(k).valid_KS(rmIdx) = [];
        tracks(k).maskI(rmIdx) = [];
        
        tracks(k).lifetime_s = trackLength*data.framerate;
    end
    
    
    % use localization coordinates from this point on
    x = tracks(k).x;
    y = tracks(k).y;
    
    trackPoints = sum(~isnan(x));
    
    %=================================
    % Determine track and gap status
    %=================================
    if (trackLength == nFrames)
        tracks(k).status = 3;
        tracks(k).lifetime_s = nFrames*data.framerate;
    else
        if (tracks(k).start>1) && (tracks(k).end<nFrames) % complete tracks
            tracks(k).status = 1;
        else
            tracks(k).status = 2; % incomplete tracks
        end
        tracks(k).lifetime_s = trackLength*data.framerate;
    end
    
    if trackPoints < trackLength % tracks has gaps
        if trackLength==1 || trackPoints==0
            tracks(k).valid = 0;
        else
            
            gacombIdx = diff(isnan(x));%diff(~isnan(x));
            gapStarts = find(gacombIdx==1)+1;
            gapEnds = find(gacombIdx==-1);
            
            % ignore 'gaps' at beginning or end
            if gacombIdx(find(gacombIdx~=0, 1, 'first')) == -1
                gapEnds = gapEnds(2:end);
            end
            if gacombIdx(find(gacombIdx~=0, 1, 'last')) == 1
                gapStarts = gapStarts(1:end-1);
            end
            gapLengths = gapEnds-gapStarts+1;
            
            segmentIdx = diff([0 ~isnan(x) 0]);
            segmentStarts = find(segmentIdx==1);
            segmentEnds = find(segmentIdx==-1)-1;
            segmentLengths = segmentEnds-segmentStarts+1;
            
            % loop over gaps
            nGaps = length(gapLengths);
            gapStatus = zeros(1,nGaps);
            for g = 1:nGaps
                % gap is valid if segments that precede & follow are > 1 frame
                % of if gap is a single frame
                if (segmentLengths(g) > 1 && segmentLengths(g+1) > 1) || gapLengths(g) == 1
                    gapStatus(g) = 4;
                else
                    gapStatus(g) = 5;
                end
            end
            
            % fill position information for valid gaps using linear interpolation
            if max(gapStatus) == 4 % only if all gaps are valid
                for g = 1:nGaps
                    borderIdx = [gapStarts(g)-1 gapEnds(g)+1];
                    gacombIdx = gapStarts(g):gapEnds(g);
                    x(gacombIdx) = interp1(borderIdx, x(borderIdx), gacombIdx);
                    y(gacombIdx) = interp1(borderIdx, y(borderIdx), gacombIdx);
                end
            end
            
            tracks(k).gapStarts = gapStarts;
            tracks(k).gapEnds = gapEnds;
            tracks(k).gapLengths = gapLengths;
            tracks(k).gapStatus = gapStatus;
            tracks(k).segmentStarts = segmentStarts;
            tracks(k).segmentEnds = segmentEnds;
            tracks(k).segmentLengths = segmentLengths;
            
            % if a single gap is invalid, the entire tracks is discarded
            tracks(k).valid = tracks(k).status == 1 && max(gapStatus) == 4;
        end
    else
        tracks(k).valid = tracks(k).status == 1;
    end
    
    %tracks(k).t = (firstIdx-1:lastIdx-1)*data.framerate;
    tracks(k).x = x;
    tracks(k).y = y;
    
    
    
    
    %==========================================
    % Fit intensity for buffer and gap frames
    %==========================================
    if tracks(k).valid
        %===============
        % buffer
        %===============
        % number of buffer frames after beginning/end cutoffs
        %         bStart = tracks(k).start - max(1, tracks(k).start-buffer);
        %         bEnd = min(data.movieLength, tracks(k).end+buffer) - tracks(k).end;
        %
        %         % start buffer
        %         mask = double(imread([data.source 'Detection' filesep 'Masks' filesep maskList(tracks(k).start)]));
        %         frameIdx = tracks(k).start+(-bStart:-1);
        %         for g = 1:length(frameIdx)
        %             xi = round(tracks(k).x(1));
        %             yi = round(tracks(k).y(1));
        %
        %             maskWindow = ~mask(yi-w:yi+w, xi-w:xi+w);
        %             c = mean(window(maskWindow));
        %
        %             frame = double(imread([data.source frameList(frameIdx(g)).name]));
        %
        %             window = frame(yi-w:yi+w, xi-w:xi+w);
        %             [p] = fitGaussian2D(window, [tracks(k).x(1)-xi tracks(k).y(1)-yi max(window(:))-c sigma c], 'A');
        %             %tracks(k).startBuffer(g) = p(3);
        %             tracks(k).startBuffer.A(1,g) = p(3);
        %             tracks(k).startBuffer.c(1,g) = c;
        %             tracks(k).startBuffer.cStd(1,g) = std(window(maskWindow));
        %         end
        %
        %         % end buffer
        %         mask = double(imread([data.source 'Detection' filesep 'Masks' filesep maskList(tracks(k).end)]));
        %         frameIdx = tracks(k).end+(1:bEnd);
        %         for g = 1:length(frameIdx)
        %             frame = double(imread([data.source frameList(frameIdx(g)).name]));
        %
        %             xi = round(tracks(k).x(end));
        %             yi = round(tracks(k).y(end));
        %
        %             maskWindow = ~mask(yi-w:yi+w, xi-w:xi+w);
        %             c = mean(window(maskWindow));
        %
        %             window = frame(yi-w:yi+w, xi-w:xi+w);
        %             [p] = fitGaussian2D(window, [tracks(k).x(end)-xi tracks(k).y(end)-yi max(window(:))-c sigma c], 'A');
        %             tracks(k).endBuffer.A(1,g) = p(3);
        %             tracks(k).endBuffer.c(1,g) = c;
        %             tracks(k).endBuffer.cStd(1,g) = std(window(maskWindow));
        %         end
        
        
        %===============
        % gaps
        %===============
        gapIdx = find(isnan(tracks(k).A(masterChannel,:)));
        if ~isempty(gapIdx)
            frameIdx = gapIdx + tracks(k).start - 1;
            
            % linear interpolation of background values %%%%%%%%%%%% change to fit?
            trackIdx = find(~isnan(tracks(k).A(masterChannel,:)));
%             k
%             trackIdx
%             gapIdx
%             tracks(k).A(masterChannel,:)
%             
            tracks(k).c(:,gapIdx) = interp1(trackIdx, tracks(k).c(:,trackIdx)', gapIdx)';
            tracks(k).cStd(:,gapIdx) = interp1(trackIdx, tracks(k).cStd(:,trackIdx)', gapIdx)';
            tracks(k).cStd_mask(:,gapIdx) = interp1(trackIdx, tracks(k).cStd_mask(:,trackIdx)', gapIdx)';
            
            
            for g = 1:length(gapIdx)
                xi = round(tracks(k).x(gapIdx(g)));
                yi = round(tracks(k).y(gapIdx(g)));
                for ch = 1:nCh
                    frame = double(imread([data.channels{ch} frameList{ch}(frameIdx(g)).name]));
                    
                    window = frame(yi-w:yi+w, xi-w:xi+w);
                    c = tracks(k).c(ch,gapIdx(g));
                    [p] = fitGaussian2D(window, [tracks(k).x(gapIdx(g))-xi tracks(k).y(gapIdx(g))-yi max(window(:))-c sigma c], 'A');
                    tracks(k).A(ch,gapIdx(g)) = p(3);
                end
            end
        end
    end
    
    %     %  for valid tracks, compute MSD and total displacement
    %     if tracks(k).valid
    %
    %         tracks(k).t = firstIdx:lastIdx;
    %
    %         tracks(k).totalD = sqrt((x(lastIdx)-x(firstIdx))^2 + (y(lastIdx)-y(firstIdx))^2);
    %
    %         % MSD
    %         xi = x(firstIdx:lastIdx);
    %         yi = y(firstIdx:lastIdx);
    %
    %         N = trackLength-1;
    %         MSDvect = zeros(1,N);
    %         for k = 1:N
    %            MSDvect(k) = mdist(xi(1+k:end), yi(1+k:end), xi(1:end-k), yi(1:end-k)) ;
    %         end
    %         L = min(N,5); % max. available points
    %         %if ~isnan(sum(MSDvect))
    %         if L > 1 % min 2 points to fit
    %             [D c alpha] = fitMSD(MSDvect(1:L), [MSDvect(L)/(4*L) 0 1], 'Dc');
    %             tracks(k).MSD = MSDvect;
    %             tracks(k).D = D;
    %             tracks(k).c = c;
    %             tracks(k).alpha = alpha;
    %         end
    %     end
    fprintf('\b\b\b\b%3d%%', round(100*k/(nTracks)));
end
fprintf('\n');


%==========================================
% Generate buffers before and after track
%==========================================
% number of buffer frames after beginning/end cutoffs
bStart = [tracks.start] - max(1, [tracks.start]-buffer);
bEnd = min(data.movieLength, [tracks.end]+buffer) - [tracks.end];

fprintf('TrackAnalysis - Generating track buffers:     ');
for ch = 1:nCh
    tifFiles = dir([data.channels{ch} '*.tif*']);
    
    sigma = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, data.markers{ch});
    
    for f = 1:data.movieLength
        frame = double(imread([data.channels{ch} tifFiles(f).name]));
        mask = double(imread([data.source 'Detection' filesep 'Masks' filesep maskList(f).name]));

        % start buffers visible in this frame
        trackIdx = find([tracks.start]-bStart<=f & f<[tracks.start] & [tracks.valid]);
        for k = trackIdx
            bi = f-tracks(k).start + bStart(k) + 1;
            xi = round(tracks(k).x(1));
            yi = round(tracks(k).y(1));
            window = frame(yi-w:yi+w, xi-w:xi+w);
            
            maskWindow = ~mask(yi-w:yi+w, xi-w:xi+w);
            c = mean(window(maskWindow));
            [p] = fitGaussian2D(window, [tracks(k).x(1)-xi tracks(k).y(1)-yi max(window(:))-c sigma c], 'A');
            
            tracks(k).startBuffer.A(ch,bi) = p(3);
            tracks(k).startBuffer.c(ch,bi) = c;
            tracks(k).startBuffer.cStd(ch,bi) = std(window(maskWindow));
        end
        
        % end buffers visible in this frame
        trackIdx = find([tracks.end]<f & f<=([tracks.end]+bEnd) & [tracks.valid]);
        for k = trackIdx
            bi = f-tracks(k).end;
            xi = round(tracks(k).x(end));
            yi = round(tracks(k).y(end));
            window = frame(yi-w:yi+w, xi-w:xi+w);
            
            maskWindow = ~mask(yi-w:yi+w, xi-w:xi+w);
            c = mean(window(maskWindow));
            [p] = fitGaussian2D(window, [tracks(k).x(end)-xi tracks(k).y(end)-yi max(window(:))-c sigma c], 'A');
            
            tracks(k).endBuffer.A(ch,bi) = p(3);
            tracks(k).endBuffer.c(ch,bi) = c;
            tracks(k).endBuffer.cStd(ch,bi) = std(window(maskWindow));
        end
        fprintf('\b\b\b\b%3d%%', round(100*((ch-1)*data.movieLength + f)/(nCh*data.movieLength)));
    end
end
fprintf('\n');





%==========================================
% Save results
%==========================================
if ~(exist([data.source 'Tracking'], 'dir')==7)
    mkdir([data.source 'Tracking']);
end
save([data.source 'Tracking' filesep 'trackAnalysis.mat'], 'tracks', 'nMergeSplit');


% loop through frames and estimate intensity
% startIdx = [tracks.start];
% endIdx = [tracks.end];
% valid = [tracks.valid];
%
% fprintf('Progress:     ');
%
% for k = 2:data.movieLength-1
%     frame = double(imread([data.source frameList(k).name]));
%     mask = double(imread([maskPath maskList(k).name]));
%
%     idx = find(valid & k>=startIdx & k<=endIdx);
%
%     for t = idx
%         relIdx = k - tracks(k).start + 1;
%
%         xi = round(tracks(k).x(relIdx));
%         yi = round(tracks(k).y(relIdx));
%
%         window = frame(yi-w:yi+w, xi-w:xi+w);
%         % binary mask
%         maskWindow = ~mask(yi-w:yi+w, xi-w:xi+w);
%         maskWindow(maskWindow~=0) = 1;
%         % background estimate
%         c = mean(mean(window(maskWindow)));
%         [p] = fitGaussian2D(window, [0 0 max(window(:))-c sigma c], 'xyA');
%         tracks(k).I(relIdx) = p(3)*2*pi*sigma^2;
%         tracks(k).c = c;
%     end
%     fprintf('\b\b\b\b%3d%%', round(100*k/(data.movieLength-2)));
% end
% fprintf('\n');



% function D = mdist(x, y, xDelta, yDelta)
% D = mean(sqrt((xDelta-x).^2 + (yDelta-y).^2));
%
%
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