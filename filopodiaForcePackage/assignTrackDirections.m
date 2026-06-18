function trackDir = assignTrackDirections(tipTracks, adhesionInfo, trkOfAdh, ...
    bodyByFrame, thetaByFrame, distByFrame, resByFrame, p)
%ASSIGNTRACKDIRECTIONS  Assign ONE shaft direction per tip track, using all
%frames jointly. Called once (not per frame) in P4.
%
% Each tip track gets a single angle (phi) determined from its entire
% lifetime. The angle is chosen by scoring candidate directions using:
%   (a) orientation consensus: mean steerable |cos(theta-phi)| sampled at
%       each tip position across all frames (stable, not per-frame noisy),
%   (b) collinear adhesion support: adhesions from OTHER tracks lying on
%       phi-direction rays across frames (shaft adhesions over time),
%   (c) length penalty: prefer more radial (shorter) shafts,
% then assembled jointly (confidence-ordered greedy) so neighboring
% filopodia keep consistent directions and shafts do not cross.
%
% Returns trackDir struct array (one per tip track that gets a direction):
%   .trackIdx   index into tipTracks
%   .ang        assigned shaft angle (rad)
%   .conf       confidence score
% Sangyoon J. Han / 2026

[H, W] = size(bodyByFrame{find(~cellfun(@isempty,bodyByFrame),1)});
maxLen   = gf(p,'MaxShaftLen',160);
sweepRng = gf(p,'SweepRange',40);
sweepStp = gf(p,'SweepStep',2);
bodyMax  = gf(p,'BodyMaxAngle',55);
band     = gf(p,'ShaftBand',4);
wShaft   = gf(p,'WShaft',0.5);
wLen     = gf(p,'WLen',0.25);
wPrior   = gf(p,'WPrior',0.6);
wOverlap = gf(p,'WOverlap',0.8);
wBaseSep = gf(p,'WBaseSep',0.7);
minSep   = gf(p,'MinBaseSep',8);
neighR   = gf(p,'NeighRadius',60);
nF       = numel(bodyByFrame);

nT = numel(tipTracks);
trackDir = struct('trackIdx',{},'ang',{},'conf',{});

%% pass A: per-track best direction (independent, all frames averaged)
cands_per_track = cell(1,nT);
score_per_track = zeros(1,nT);
ang_per_track   = zeros(1,nT);
conf_per_track  = zeros(1,nT);  % for ordering

for i = 1:nT
    tr = tipTracks(i);
    % collect tip positions across frames + steerable maps
    tipXY = tr.pos;    % Nframes x 2 [x y]
    frames = tr.frames;
    % --- candidate directions: seed from time-averaged theta at tip ---
    thAll = zeros(numel(frames),1);
    resAll = zeros(numel(frames),1);
    gXAll = zeros(numel(frames),1); gYAll = zeros(numel(frames),1);
    for j = 1:numel(frames)
        t = frames(j);
        if t<1||t>nF||isempty(thetaByFrame{t}), continue; end
        ix = min(max(round(tipXY(j,1)),1),W); iy = min(max(round(tipXY(j,2)),1),H);
        thAll(j)  = thetaByFrame{t}(iy,ix);
        resAll(j) = resByFrame{t}(iy,ix);
        [gX,gY] = gradient(distByFrame{t});
        gXAll(j) = gX(iy,ix); gYAll(j) = gY(iy,ix);
    end
    % mean body-ward direction across frames
    bwX = mean(-gXAll); bwY = mean(-gYAll);
    bAng = atan2(bwY,bwX);
    % circular mean of theta (pick bodyward sense)
    th_mean_s = mean(sin(thAll)); th_mean_c = mean(cos(thAll));
    th_mean = atan2(th_mean_s,th_mean_c);
    for ang = (th_mean:pi:th_mean+pi)
        if cos(ang-bAng)>0, th_seed = ang; break; end
    end
    if ~exist('th_seed','var'), th_seed = bAng; end
    seeds = [th_seed, bAng];
    cands = [];
    for sd = seeds
        cands = [cands, sd+deg2rad(-sweepRng:sweepStp:sweepRng)]; %#ok<AGROW>
    end
    cands = cands(cos(cands-bAng)>cos(deg2rad(bodyMax)));
    cands = unique(round(cands*1e3)/1e3);
    cands_per_track{i} = cands;

    % score each candidate: time-averaged orientation consensus + shaft support
    bestSc = -inf; bestAng = bAng;
    for ang = cands
        osc_sum = 0; osc_n = 0; shaft_sum = 0; hit_n = 0;
        for j = 1:numel(frames)
            t = frames(j); if t<1||t>nF||isempty(bodyByFrame{t}), continue; end
            x0 = tipXY(j,1); y0 = tipXY(j,2);
            [osc, hit] = rayConsensus(x0,y0,ang,thetaByFrame{t},bodyByFrame{t},maxLen,H,W);
            if isempty(hit), continue; end
            osc_sum = osc_sum+osc; osc_n = osc_n+1; hit_n = hit_n+1;
            % collinear adhesion support at this frame
            adh = adhesionInfo{t};
            if ~isempty(adh)
                Af = cat(1,adh.pos); dAf = [adh.dist];
                di = distByFrame{t}(min(max(round(y0),1),H),min(max(round(x0),1),W));
                u=[cos(ang),sin(ang)]; rel=Af-[x0,y0];
                proj=rel(:,1)*u(1)+rel(:,2)*u(2);
                perp=abs(-rel(:,1)*u(2)+rel(:,2)*u(1));
                onl=(perp<band)&(proj>1)&(proj<hit(3))&(dAf(:)<di);
                % only count adhesions from OTHER tracks (not this tip)
                for k2=find(onl)'
                    ot=trkOfAdh{t}(k2);
                    if ot~=i, shaft_sum=shaft_sum+1; end
                end
            end
        end
        if hit_n<1, continue; end
        sc = osc_sum/osc_n + wShaft*shaft_sum/max(1,hit_n) - wLen*(1/hit_n)*sum(cellfun(@(t2) ...
            rayLen(tipXY(find(frames==t2,1),:),ang,bodyByFrame{t2},maxLen,H,W), ...
            num2cell(frames(arrayfun(@(jj) ~isempty(bodyByFrame{frames(jj)}),1:numel(frames)))))/maxLen);
        if sc>bestSc, bestSc=sc; bestAng=ang; end
    end
    ang_per_track(i)  = bestAng;
    score_per_track(i)= bestSc;
    % confidence = mean steerable res at tip + score
    conf_per_track(i) = mean(resAll(resAll>0)) + max(bestSc,0);
    clear th_seed
end

%% pass B: joint greedy (confidence order, neighbor prior + no-crossing + base-sep)
[~, ord] = sort(conf_per_track,'descend');
fixedXY=[]; fixedAng=[]; fixedSeg={}; fixedBase=zeros(0,2);

for oo = 1:nT
    i = ord(oo);
    tr = tipTracks(i);
    cands = cands_per_track{i};
    if isempty(cands), continue; end
    % neighbor prior from already-fixed tracks
    tipMean = mean(tr.pos,1);
    prior = neighborPrior(tipMean, fixedXY, fixedAng, neighR);
    % re-score with prior + overlap + base-sep using mean tip position
    bestSc=-inf; bestAng=ang_per_track(i);
    % representative frame: most central
    [~,jrep] = min(abs(tr.frames - median(tr.frames)));
    t = tr.frames(jrep);
    if isempty(bodyByFrame{t}), t=tr.frames(1); end
    x0=tr.pos(jrep,1); y0=tr.pos(jrep,2);
    for ang = cands
        [osc,hit] = rayConsensus(x0,y0,ang,thetaByFrame{t},bodyByFrame{t},maxLen,H,W);
        if isempty(hit), continue; end
        sc = osc - wLen*(hit(3)/maxLen);
        if ~isempty(prior)
            sc = sc - wPrior*abs(atan2(sin(ang-prior),cos(ang-prior)))/pi;
        end
        for q=1:numel(fixedSeg)
            S=fixedSeg{q};
            if segCross([x0 y0],hit(1:2),S(1,:),S(2,:)), sc=sc-wOverlap; break; end
        end
        if ~isempty(fixedBase) && wBaseSep>0
            d2b=(fixedBase(:,1)-hit(1)).^2+(fixedBase(:,2)-hit(2)).^2;
            if min(d2b)<minSep^2, sc=sc-wBaseSep; end
        end
        if sc>bestSc, bestSc=sc; bestAng=ang; bestHit=hit; end
    end
    if ~exist('bestHit','var')||isempty(bestHit), continue; end
    k=numel(trackDir)+1;
    trackDir(k).trackIdx = i;
    trackDir(k).ang      = bestAng;
    trackDir(k).conf     = conf_per_track(i);
    fixedXY(end+1,:)   = tipMean; %#ok<AGROW>
    fixedAng(end+1)    = bestAng; %#ok<AGROW>
    fixedSeg{end+1}    = [x0 y0; bestHit(1:2)]; %#ok<AGROW>
    fixedBase(end+1,:) = bestHit(1:2); %#ok<AGROW>
    clear bestHit
end
end

% ===================================================================
function [osc,hit] = rayConsensus(x0,y0,ang,theta,bodyMask,maxLen,H,W)
oss=0;nO=0;hit=[];
for s=1:1.5:maxLen
    x=x0+s*cos(ang);y=y0+s*sin(ang);ix=round(x);iy=round(y);
    if ix<1||ix>W||iy<1||iy>H,break;end
    if bodyMask(iy,ix),hit=[x y s];break;end
    oss=oss+abs(cos(theta(iy,ix)-ang));nO=nO+1;
end
if isempty(hit)||nO<3,osc=-inf;hit=[];else,osc=oss/nO;end
end

function L = rayLen(tipXY,ang,bodyMask,maxLen,H,W)
L=maxLen; x0=tipXY(1);y0=tipXY(2);
if isempty(bodyMask),return;end
for s=1:1.5:maxLen
    x=x0+s*cos(ang);y=y0+s*sin(ang);ix=round(x);iy=round(y);
    if ix<1||ix>size(bodyMask,2)||iy<1||iy>H,break;end
    if bodyMask(iy,ix),L=s;return;end
end
end

function pr=neighborPrior(tip,fixedXY,fixedAng,R)
pr=[];if isempty(fixedXY),return;end
d2=(fixedXY(:,1)-tip(1)).^2+(fixedXY(:,2)-tip(2)).^2;
w=exp(-d2/(2*R^2));if sum(w)<1e-3,return;end
pr=atan2(sum(w.*sin(fixedAng(:))),sum(w.*cos(fixedAng(:))));
end

function tf=segCross(p1,p2,p3,p4)
tf=(ccw(p1,p3,p4)~=ccw(p2,p3,p4))&&(ccw(p1,p2,p3)~=ccw(p1,p2,p4));
end
function t=ccw(a,b,c),t=(c(2)-a(2))*(b(1)-a(1))>(b(2)-a(2))*(c(1)-a(1));end
function v=gf(s,n,d),if isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end,end
