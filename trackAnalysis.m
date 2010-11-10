
% The only tracks of interest are those with status 1 or 4.

% Francois Aguet, June 2010 (revised from 'determineLifetimeStatus.m')

function [tracks nMergeSplit] = trackAnalysis(data, buffer, filename)

if nargin<2
    buffer = 5;
end
if nargin<3
    filename = 'trackedFeatures.mat';
end

ny = data.imagesize(1);
nx = data.imagesize(2);
nFrames = data.movieLength;
frameList = dir([data.source '*.tif*']);
maskList = dir([data.source 'Detection' filesep 'Masks' filesep '*.tif']);

load([data.source 'Detection' filesep 'detectionResults.mat']);

%=================================
% Determine master/slave channels
%=================================
nChannels = length(data.channels);
% exclude master from list of channels
masterChannel = regexp(data.source, data.channels);
masterChannel = find([masterChannel{:}]);
slaveChannels = setdiff(1:nChannels, masterChannel);
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
     
elseif exist([data.source 'TrackInfoMatrices' filesep 'trackInfo.mat'], 'file')==2
    % (for old tracker. oldest version: trackInfo.mat)
    trackinfo = load([data.source 'TrackInfoMatrices' filesep 'trackedFeatures.mat']);
    trackinfo = trackinfo.trackedFeatureInfo;
    nTracks = size(trackinfo, 1);
else
    error('No valid tracker output found.');
end
tracks(1:nTracks) = struct('t', [], 'x', [], 'y', [], 'A', [], 'maskI', [], 'c', [], 'cStd', [],...
    'status', [], 'gapStatus', [],...
    'gapStarts', [], 'gapEnds', [], 'gapLengths', [],...
    'segmentStarts', [], 'segmentEnds', [], 'segmentLengths', [],...
    'D', [], 'alpha', [], 'MSDvect', [], 'totalD', [],...
    'seqOfEvents', [], 'tracksFeatIndxCG', []);


fprintf('Progress:     ');
for k = 1:nTracks
    if ~isstruct(trackinfo)
        x = trackinfo(k,1:8:end);
        y = trackinfo(k,2:8:end);
        
        % index of detected track points
        trackIdx = find(~isnan(x));
        
        firstIdx = trackIdx(1);
        lastIdx = trackIdx(end);
        trackLength = lastIdx-firstIdx+1;
        trackPoints = length(trackIdx);
        
        tracks(k).start = firstIdx;
        tracks(k).end = lastIdx;
    else
        % convert/assign structure fields
        tracks(k).x = trackinfo(k).tracksCoordAmpCG(1:8:end);
        tracks(k).y = trackinfo(k).tracksCoordAmpCG(2:8:end);
        tracks(k).maskI = trackinfo(k).tracksCoordAmpCG(4:8:end);
        tracks(k).seqOfEvents = trackinfo(k).seqOfEvents;
        tracks(k).tracksFeatIndxCG = trackinfo(k).tracksFeatIndxCG;
        
        tracks(k).lifetime = length(tracks(k).x);
        tracks(k).t = trackinfo(k).seqOfEvents(1,1):trackinfo(k).seqOfEvents(2,1);
        
        x = tracks(k).x;
        y = tracks(k).y;
        
        firstIdx = trackinfo(k).seqOfEvents(1,1);
        lastIdx = trackinfo(k).seqOfEvents(end,1);
        
        trackLength = lastIdx-firstIdx+1;
        trackPoints = sum(~isnan(x));
        
        tracks(k).start = firstIdx;
        tracks(k).end = lastIdx;
    end
    
    
    
    % determine tracks characteristics
    if (trackLength == nFrames)
        tracks(k).status = 3;
        tracks(k).lifetime = nFrames;
    else
        if (firstIdx>1) && (lastIdx<nFrames) % complete tracks
            tracks(k).status = 1;
        else
            tracks(k).status = 2; % incomplete tracks
        end
        tracks(k).lifetime = trackLength;
    end
    
    if trackPoints < trackLength % tracks has gaps
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
    else
        tracks(k).valid = tracks(k).status == 1;
    end
    
    
    tracks(k).t = firstIdx:lastIdx;
    if ~isstruct(trackinfo)
        tracks(k).x = x(tracks(k).t);
        tracks(k).y = y(tracks(k).t);
        tracks(k).xv = x;
        tracks(k).yv = y;
    else
        tracks(k).x = x;
        tracks(k).y = y;
    end
    
    %==============================================================================
    % Read amplitude & background from detectionResults.mat (localization results)
    %==============================================================================
    tracks(k).A = NaN(size(x)); % NaNs indicate gaps
    tracks(k).c = NaN(size(x));
    tracks(k).cStd = NaN(size(x));
    
    frameRange = tracks(k).start:tracks(k).end;
    if isstruct(trackinfo)        
        for i = 1:length(frameRange)
            idx = tracks(k).tracksFeatIndxCG(i);
            if idx ~= 0
                tracks(k).A(i) = frameInfo(frameRange(i)).A(idx);
                tracks(k).maskI(i) = frameInfo(frameRange(i)).amp(idx);
                tracks(k).c(i) = frameInfo(frameRange(i)).c(idx);
                tracks(k).cStd(i) = frameInfo(frameRange(i)).cStd(idx);
            end
        end
    else
        % trackedFeatureNum
        %%%%%%%%%%%%%%%%%%%%%
    end
    
    % invalidate if track is within frame border
    xi = round(x);
    yi = round(y);
    if min(xi)<=w || max(xi)>nx-w || min(yi)<=w || max(yi)>ny-w
        tracks(k).valid = 0;
    end
    
    
    %==========================================
    % Fit intensity for buffer and gap frames
    %==========================================
    if tracks(k).valid
        %===============
        % buffer
        %===============
        % number of buffer frames after beginning/end cutoffs
        bStart = tracks(k).start - max(1, tracks(k).start-buffer);
        bEnd = min(data.movieLength, tracks(k).end+buffer) - tracks(k).end;
        
        % start buffer
        c = tracks(k).c(1);
        frameIdx = tracks(k).start+(-bStart:-1);
        for g = 1:length(frameIdx)
            frame = double(imread([data.source frameList(frameIdx(g)).name]));
            xi = round(tracks(k).x(1));
            yi = round(tracks(k).y(1));
            window = frame(yi-w:yi+w, xi-w:xi+w);
            [p] = fitGaussian2D(window, [tracks(k).x(1)-xi tracks(k).y(1)-yi max(window(:))-c sigma c], 'A');
            tracks(k).startBuffer(g) = p(3);
        end
        
        % end buffer
        c = tracks(k).c(end);
        frameIdx = tracks(k).end+(1:bEnd);
        for g = 1:length(frameIdx)
            frame = double(imread([data.source frameList(frameIdx(g)).name]));
            xi = round(tracks(k).x(end));
            yi = round(tracks(k).y(end));
            window = frame(yi-w:yi+w, xi-w:xi+w);
            [p] = fitGaussian2D(window, [tracks(k).x(end)-xi tracks(k).y(end)-yi max(window(:))-c sigma c], 'A');
            tracks(k).endBuffer(g) = p(3);
        end
        
        
        %===============
        % gaps
        %===============
        gapIdx = find(isnan(tracks(k).A));
        frameIdx = gapIdx + tracks(k).start - 1;
        
        % linear interpolation of background values
        trackIdx = find(~isnan(tracks(k).c));
        tracks(k).c(gapIdx) = interp1(trackIdx, tracks(k).c(trackIdx), gapIdx);
        tracks(k).cStd(gapIdx) = interp1(trackIdx, tracks(k).cStd(trackIdx), gapIdx);
        
        for g = 1:length(gapIdx)
            frame = double(imread([data.source frameList(frameIdx(g)).name]));
            
            xi = round(tracks(k).x(gapIdx(g)));
            yi = round(tracks(k).y(gapIdx(g)));
            window = frame(yi-w:yi+w, xi-w:xi+w);
            [p] = fitGaussian2D(window, [tracks(k).x(gapIdx(g))-xi tracks(k).y(gapIdx(g))-yi max(window(:))-tracks(k).c(gapIdx(g)) sigma tracks(k).c(gapIdx(g))], 'A');
            tracks(k).A(gapIdx(g)) = p(3);
        end
    end
    

    
    
% % % % % %     k = 500;
% % % % % %     idx = tf(k).tracksFeatIndxCG;
% % % % % %     
% % % % % %     frameRange = tf(k).start:tf(k).end;
% % % % % %     
% % % % % %     for i = 1:length(frameRange)
% % % % % %         tmp(i) = frameInfo(frameRange(i)).xcom(idx(i));
% % % % % %     end
% % % % % %     
% % % % % %     figure; plot(tf(k).x)
% % % % % %     hold on;
% % % % % %     plot(tmp, 'r--');
    
    
    
    
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
% Read out intensity in slave channels
%==========================================
for ch = slaveChannels
    cPath = data.channels{ch};
    tifFiles = dir([cPath '*.tif*']);

    sigma = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, data.markers{ch});
    %w = ceil(3*sigma); -> move to start, max w for all sigma

    for f = 1:data.movieLength
        frame = double(imread([cPath tifFiles(f).name]));
        mask = double(imread([data.source 'Detection' filesep 'Masks' filesep maskList(f).name]));
        
        % tracks visible in this frame
        trackIdx = find([tracks.start]<=f & f<=[tracks.end] & [tracks.valid]);
        for k = trackIdx
            fi = f-tracks(k).start+1;
            xi = round(tracks(k).x(fi));
            yi = round(tracks(k).y(fi));
            window = frame(yi-w:yi+w, xi-w:xi+w);
            
            maskWindow = ~mask(yi-w:yi+w, xi-w:xi+w);
            maskWindow(maskWindow~=0) = 1;
            % background estimate
            c = mean(window(maskWindow));
            [p] = fitGaussian2D(window, [tracks(k).x(fi)-xi tracks(k).y(fi)-yi max(window(:))-c sigma c], 'A');

            tracks(k).A(ch,fi) = p(3);
            tracks(k).c(ch,fi) = c;
            tracks(k).cStd(ch,fi) = std(window(maskWindow));
        end
    end
end


%==========================================
% Save results
%==========================================
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



function D = mdist(x, y, xDelta, yDelta)
D = mean(sqrt((xDelta-x).^2 + (yDelta-y).^2));


function [D c alpha] = fitMSD(MSD, prmVect, prmSel)

opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-6, ...
    'Tolfun', 1e-6);


estIdx = false(1,3); % [D c alpha]
estIdx(regexp('Dca', ['[' prmSel ']'])) = true;

x = lsqnonlin(@costMSD, prmVect(estIdx), [], [], opts, MSD, prmVect, estIdx);
prmVect(estIdx) = deal(abs(x));

D = prmVect(1);
c = prmVect(2); % constant offset
alpha = prmVect(3);



function [v] = costMSD(p, MSD, prmVect, estIdx)
prmVect(estIdx) = deal(abs(p));

D = prmVect(1);
c = prmVect(2); % constant offset
alpha = prmVect(3);
%d % dimensionality

t = 1:length(MSD);

v = MSD - 4*D*t.^alpha - c;