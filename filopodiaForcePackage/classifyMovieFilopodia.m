function classifyMovieFilopodia(movieData)
%CLASSIFYMOVIEFILOPODIA  Process 4. Assemble filopodia with temporally
%consistent directions.
%
% Two-pass approach:
%   Pass A (track-level, once): Assign ONE shaft direction per well-tracked
%     tip track, using orientation consensus + shaft adhesion support pooled
%     across ALL frames of that track. Directions are assigned jointly
%     (confidence-ordered greedy) so neighboring filopodia stay consistent
%     and shafts do not cross over the full movie.
%   Pass B (frame-level): For each frame, project each tip's fixed direction
%     to the current body mask -> base(t), L(t), shaft adhesions(t).
%     Base moves slightly as the body mask changes, L(t) = tip-to-base length.
%
% This gives temporally continuous L(t) traces (same filopodium, same
% direction, base only changes as body remodels) suitable for measuring
% extension/retraction. Cross-frame chaotic direction changes are eliminated.
% Sangyoon J. Han / 2026

%% process & params
iProc = movieData.getProcessIndex('FilopodiaClassificationProcess', 1, 0);
assert(~isempty(iProc), 'No FilopodiaClassificationProcess found.');
proc = movieData.processes_{iProc};
p = parseProcessParams(proc);
iChan = p.ChannelIndex;

iSeg = p.SegProcessIndex; if isempty(iSeg), iSeg = movieData.getProcessIndex('FilopodiaSegmentationProcess',1,0); end
iDet = p.DetProcessIndex; if isempty(iDet), iDet = movieData.getProcessIndex('FilopodiaDetectionProcess',1,0); end
iTrk = p.TrkProcessIndex; if isempty(iTrk), iTrk = movieData.getProcessIndex('FilopodiaTrackingProcess',1,0); end
assert(~isempty(iSeg)&&~isempty(iDet)&&~isempty(iTrk), 'Processes 1-3 must be run first.');
segProc = movieData.processes_{iSeg};
detProc = movieData.processes_{iDet};
trkProc = movieData.processes_{iTrk};
segOutDir = segProc.funParams_.OutputDirectory;

detFile = detProc.outFilePaths_{1,iChan};
if isempty(detFile)||exist(detFile,'file')~=2, detFile = fullfile(detProc.funParams_.OutputDirectory,'filoDetection.mat'); end
Sdet = load(detFile,'adhesionInfo'); adhesionInfo = Sdet.adhesionInfo;

trkFile = trkProc.outFilePaths_{1,iChan};
if isempty(trkFile)||exist(trkFile,'file')~=2, trkFile = fullfile(trkProc.funParams_.OutputDirectory,'filoTracks.mat'); end
Strk = load(trkFile,'adhesionTracks','tracksFinal');
adhesionTracks = Strk.adhesionTracks; tracksFinalAll = Strk.tracksFinal;

%% I/O
inFilePaths = cell(1,numel(movieData.channels_)); inFilePaths{1,iChan} = trkFile;
proc.setInFilePaths(inFilePaths);
outDir = p.OutputDirectory; mkClrDir(outDir);
outFile = fullfile(outDir,'filoClassification.mat');
outFilePaths = cell(1,numel(movieData.channels_)); outFilePaths{1,iChan} = outFile;
proc.setOutFilePaths(outFilePaths);

%% preload P1 maps
nF = movieData.nFrames_;
bodyByFrame  = cell(1,nF); thetaByFrame = cell(1,nF);
distByFrame  = cell(1,nF); resByFrame   = cell(1,nF);
shaftByFrame_ = cell(1,nF);   % P1 shaftMask (ridge binary)
for t = 1:nF
    fn = fullfile(segOutDir,sprintf('filoSeg_frame_%04d.mat',t));
    if exist(fn,'file')~=2, continue; end
    s = load(fn,'bodyMask','theta','res','shaftMask');
    bodyByFrame{t}=s.bodyMask; thetaByFrame{t}=s.theta;
    distByFrame{t}=bwdist(s.bodyMask); resByFrame{t}=s.res;
    if isfield(s,'shaftMask')
        shaftByFrame_{t} = s.shaftMask & ~s.bodyMask;
    else
        % fallback: threshold res at median+1*mad
        rv = s.res(~s.bodyMask); thr = median(rv)+1.4826*mad(rv,1);
        shaftByFrame_{t} = s.res > thr & ~s.bodyMask;
    end
end

%% tip eligibility
nTr = numel(adhesionTracks);
minLife   = getfielddef(p,'MinTipLifetime',5);
linFrac   = getfielddef(p,'MinLinearFrac',0.85);
minSpread = getfielddef(p,'MinTrajSpread',3);
minTipDist= getfielddef(p,'MinTipDist',6);
tipEligible = false(1,nTr);
for i = 1:nTr
    tr = adhesionTracks(i);
    if numel(tr.frames)<minLife || max(tr.dist)<minTipDist, continue; end
    P2 = tr.pos-mean(tr.pos,1);
    ev = sort(eig((P2'*P2)/max(1,size(P2,1)-1)),'descend');
    if sqrt(sum(ev))<minSpread || ev(1)/max(sum(ev),eps)>=linFrac
        tipEligible(i) = true;
    end
end
tipIdx = find(tipEligible);
tipTracks = adhesionTracks(tipIdx);

%% trkOfAdh: adhesion k at frame t -> track index (into adhesionTracks)
% adhesionInfo may be shorter than movieData.nFrames_ if detection (P2) did
% not cover every frame; clamp to what is actually available.
nAdhF = numel(adhesionInfo);
if nAdhF < nF
    warning('FilopodiaClassification:frameMismatch', ...
        'adhesionInfo has %d frames but movie has %d; using %d.', nAdhF, nF, nAdhF);
end
nFuse = min(nF, nAdhF);
trkOfAdh = cell(1,nFuse);
for t=1:nFuse
    if isempty(adhesionInfo{t}), trkOfAdh{t}=zeros(1,0);
    else, trkOfAdh{t}=zeros(1,numel(adhesionInfo{t})); end
end
for i=1:nTr
    tr=adhesionTracks(i);
    for j=1:numel(tr.frames)
        t=tr.frames(j); k=tr.featIdx(j);
        if t>=1&&t<=nFuse&&k>=1&&k<=numel(trkOfAdh{t}), trkOfAdh{t}(k)=i; end
    end
end

%% PASS A: assign ONE direction per tip track (all frames jointly)
fprintf('Pass A: assigning directions to %d tip tracks...\n', numel(tipTracks));
progressText(0,'Pass A: direction assignment','Filopodia P4 pass A');
trackDir = assignTrackDirections(tipTracks, adhesionInfo, trkOfAdh, ...
    bodyByFrame, thetaByFrame, distByFrame, resByFrame, shaftByFrame_, p);
progressText(1,'Pass A: direction assignment');
fprintf('  -> %d filopodia assigned directions\n', numel(trackDir));

%% PASS B: per-frame projection -> base(t), L(t), shaft adhesions(t)
dt  = movieData.timeInterval_; if isempty(dt), dt=1; end
pix = movieData.pixelSize_;    if isempty(pix), pix=1; end
sw  = max(1,round(getfielddef(p,'VelSmoothWin',3)));
reachFrac = getfielddef(p,'MinReachFrac',0.5);
band      = getfielddef(p,'ShaftBand',4);

posByFrame.tip   = cell(1,nF);
posByFrame.base  = cell(1,nF);
posByFrame.shaft = cell(1,nF);
shaftByFrame     = cell(1,nF);
for t=1:nF
    posByFrame.tip{t}=zeros(0,2); posByFrame.base{t}=zeros(0,2);
    posByFrame.shaft{t}=zeros(0,2); shaftByFrame{t}=zeros(0,2);
end

nFil = numel(trackDir);
filTipFrames = cell(1,nFil); filTipLen = cell(1,nFil);
filTipPos = cell(1,nFil);    filBase   = cell(1,nFil);
H = size(bodyByFrame{find(~cellfun(@isempty,bodyByFrame),1)},1);
W = size(bodyByFrame{find(~cellfun(@isempty,bodyByFrame),1)},2);
maxLen = getfielddef(p,'MaxShaftLen',200);

progressText(0,'Pass B: per-frame projection','Filopodia P4 pass B');
for f = 1:nFil
    td = trackDir(f);
    tr = tipTracks(td.trackIdx);
    ang = td.ang;
    nreach = 0;
    for j = 1:numel(tr.frames)
        t = tr.frames(j);
        if isempty(bodyByFrame{t}), continue; end
        x0 = tr.pos(j,1); y0 = tr.pos(j,2);
        % shoot ray at fixed angle, find body
        base=[]; Lb=NaN;
        for s = 1:1.5:maxLen
            x=x0+s*cos(ang); y=y0+s*sin(ang);
            ix=round(x); iy=round(y);
            if ix<1||ix>W||iy<1||iy>H, break; end
            if bodyByFrame{t}(iy,ix), base=[x y]; Lb=s; break; end
        end
        if isempty(base), continue; end
        nreach = nreach+1;
        % shaft adhesions: other adhesions on the ray, nearer body
        if t <= numel(adhesionInfo), adh = adhesionInfo{t}; else, adh = []; end
        shIdx = [];
        if ~isempty(adh)
            Af=cat(1,adh.pos); dAf=[adh.dist];
            di=distByFrame{t}(min(max(round(y0),1),H),min(max(round(x0),1),W));
            u=[cos(ang),sin(ang)]; rel=Af-[x0,y0];
            proj=rel(:,1)*u(1)+rel(:,2)*u(2);
            perp=abs(-rel(:,1)*u(2)+rel(:,2)*u(1));
            onl=(perp<band)&(proj>1)&(proj<Lb)&(dAf(:)<di);
            shIdx=find(onl)';
        end
        filTipFrames{f}(end+1)=t; filTipLen{f}(end+1)=Lb;
        filTipPos{f}(end+1,:)=[x0,y0]; filBase{f}(end+1,:)=base;
        posByFrame.tip{t}(end+1,:)=[x0,y0];
        posByFrame.base{t}(end+1,:)=base;
        if ~isempty(shIdx)
            posByFrame.shaft{t}=[posByFrame.shaft{t}; Af(shIdx,:)];
        end
        nstep=max(2,round(Lb)); ss=linspace(0,Lb,nstep)';
        line=[x0+ss*cos(ang), y0+ss*sin(ang)];
        shaftByFrame{t}=[shaftByFrame{t}; line; nan(1,2)];
    end
    trackDir(f).nReach = nreach;
    trackDir(f).reachFrac = nreach/max(1,numel(tr.frames));
    progressText(f/nFil,'Pass B: per-frame projection');
end

%% build filopodia struct (accepted tracks)
roleByTrack = repmat({'background'},1,nTr);
filopodia = struct('tipTrackId',{},'frames',{},'tipPos',{},'basePos',{}, ...
    'L',{},'L_nm',{},'velocity_nmps',{},'lifetime',{},'ang',{});
for f = 1:nFil
    td = trackDir(f);
    if isempty(filTipFrames{f}) || td.reachFrac < reachFrac, continue; end
    globalIdx = tipIdx(td.trackIdx);
    roleByTrack{globalIdx} = 'tip';
    L = filTipLen{f}; Lsm = smoothLocal(L,sw);
    n = numel(filopodia)+1;
    filopodia(n).tipTrackId    = adhesionTracks(globalIdx).trackId;
    filopodia(n).frames        = filTipFrames{f};
    filopodia(n).tipPos        = filTipPos{f};
    filopodia(n).basePos       = filBase{f};
    filopodia(n).L             = L;
    filopodia(n).L_nm          = L*pix;
    filopodia(n).velocity_nmps = gradient(Lsm)/dt*pix;
    filopodia(n).lifetime      = numel(filTipFrames{f});
    filopodia(n).ang            = td.ang;
end

% tip tracks (u-track format) for TracksDisplay
acceptIdx = find(strcmp(roleByTrack,'tip'));
if ~isempty(acceptIdx)&&~isempty(tracksFinalAll)
    safe=acceptIdx(acceptIdx<=numel(tracksFinalAll));
    tipTrks=tracksFinalAll(safe);
else, tipTrks=tracksFinalAll([]); end

tipTracks = tipTrks; %#ok
save(outFile,'filopodia','tipTracks','roleByTrack','posByFrame','shaftByFrame','-v7.3');
fprintf('Classification done: %d filopodia (of %d tip-eligible tracks).\n', numel(filopodia), numel(tipIdx));
end

% ===================================================================
function v = getfielddef(s,name,default)
if isfield(s,name)&&~isempty(s.(name)), v=s.(name); else, v=default; end
end
function y = smoothLocal(x,w)
if w<=1||numel(x)<3, y=x; return; end
k=ones(1,w)/w; y=conv(x,k,'same');
h=floor(w/2);
for i=1:min(h,numel(x))
    y(i)=mean(x(1:min(numel(x),i+h)));
    y(end-i+1)=mean(x(max(1,end-i+1-h):end));
end
end
