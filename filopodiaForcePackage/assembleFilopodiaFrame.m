function fil = assembleFilopodiaFrame(tipIdx, A, distA, theta, bodyMask, distOut, res, p)
%ASSEMBLEFILOPODIAFRAME  Assign straight, mutually-consistent shafts to tips.
%
%   fil = assembleFilopodiaFrame(tipIdx, A, distA, theta, bodyMask, distOut, res, p)
%
% For one frame, decide each tip's STRAIGHT shaft direction and base on the
% body, jointly so that neighboring filopodia stay consistent and do not
% overlap. Filopodia are straight (formin-driven), spatially-neighboring ones
% share orientation, and orientation varies smoothly along the cell edge.
%
% Algorithm:
%   pass 1 - each tip's independent best direction -> provisional shaft count.
%            confidence = z(steerable res at tip) + z(provisional shaft count).
%   pass 2 - process tips from HIGH to LOW confidence (most certain first).
%            For each, choose the direction maximizing:
%              orientation consensus (steerable along the ray)
%            + WShaft * (# unclaimed shaft adhesions on the ray)   [deep base]
%            - WLen   * (shaft length / MaxShaftLen)               [small now]
%            - WPrior * angular deviation from already-fixed neighbors
%            - WOverlap if the segment crosses an already-fixed shaft
%            restricted to within BodyMaxAngle of the body-ward direction.
%            Adhesions on the chosen ray (perp < ShaftBand, nearer the body,
%            unclaimed) become that filopodium's shaft adhesions (claimed).
%
% Inputs: tipIdx (indices into A of eligible tips), A (Nx2 [x y] all adhesions
% this frame), distA (signed dist to body), theta/bodyMask/distOut/res (P1
% maps). Returns struct array fil with fields: tipPos, ang, base, shaftIdx,
% len, reached, conf.
% Sangyoon J. Han / 2026

[H, W] = size(bodyMask);
maxLen   = gf(p,'MaxShaftLen',160);
sweepRng = gf(p,'SweepRange',40);
sweepStp = gf(p,'SweepStep',2);
bodyMax  = gf(p,'BodyMaxAngle',55);
band     = gf(p,'ShaftBand',4);
wShaft   = gf(p,'WShaft',0.5);     % deep-base reward (passes many shaft adhesions)
wLen     = gf(p,'WLen',0.25);      % small length penalty
wPrior   = gf(p,'WPrior',0.6);     % neighbor-direction consistency / angle continuity
wOverlap = gf(p,'WOverlap',0.8);   % no-crossing penalty
wBaseSep = gf(p,'WBaseSep',0.7);   % penalty for a base near an already-fixed base
minBaseSep = gf(p,'MinBaseSep',8); % px; bases closer than this are penalized
neighR   = gf(p,'NeighRadius',60);

[gX, gY] = gradient(distOut);
fil = struct('tipPos',{},'ang',{},'base',{},'shaftIdx',{},'len',{},'reached',{},'conf',{});
if isempty(tipIdx), return; end

% ---- pass 1: provisional independent direction + shaft count ----
nT = numel(tipIdx);
ang0 = zeros(nT,1); nsh0 = zeros(nT,1); resT = zeros(nT,1);
claimed = false(size(A,1),1);
for a = 1:nT
    i = tipIdx(a);
    resT(a) = res(min(max(round(A(i,2)),1),H), min(max(round(A(i,1)),1),W));
    [~, ang0(a), hit, onl] = bestDir(A(i,:), distA(i), A, distA, claimed, theta, bodyMask, ...
        distOut, gX, gY, maxLen, sweepRng, sweepStp, bodyMax, band, wShaft, wLen, ...
        [], 0, {}, wOverlap, zeros(0,2), 0, minBaseSep, H, W);
    if ~isempty(hit), nsh0(a) = nnz(onl); end
end
% confidence = z(res) + z(shaft count)
conf = zsc(resT) + zsc(nsh0);
[~, ord] = sort(conf, 'descend');

% ---- pass 2: greedy in confidence order, with neighbor prior + overlap ----
claimed = false(size(A,1),1);
fixedXY = zeros(0,2); fixedAng = []; fixedSeg = {}; fixedBase = zeros(0,2);   % fixed shafts so far
for oo = 1:nT
    a = ord(oo);
    i = tipIdx(a);
    if claimed(i), continue; end       % already absorbed as another tip's shaft
    prior = neighborPrior(A(i,:), fixedXY, fixedAng, neighR);
    [sc, ang, hit, onl] = bestDir(A(i,:), distA(i), A, distA, claimed, theta, bodyMask, ...
        distOut, gX, gY, maxLen, sweepRng, sweepStp, bodyMax, band, wShaft, wLen, ...
        prior, wPrior, fixedSeg, wOverlap, fixedBase, wBaseSep, minBaseSep, H, W);
    if isempty(hit), continue; end
    shIdx = find(onl);
    claimed(shIdx) = true; claimed(i) = true;

    k = numel(fil)+1;
    fil(k).tipPos   = A(i,:);
    fil(k).ang      = ang;
    fil(k).base     = hit(1:2);
    fil(k).shaftIdx = shIdx;
    fil(k).len      = hit(3);
    fil(k).reached  = true;
    fil(k).conf     = conf(a);
    fixedXY(end+1,:) = A(i,:); %#ok<AGROW>
    fixedAng(end+1)  = ang;    %#ok<AGROW>
    fixedSeg{end+1}  = [A(i,:); hit(1:2)]; %#ok<AGROW>
    fixedBase(end+1,:) = hit(1:2); %#ok<AGROW>
end
end

% ===================================================================
function [best, bestAng, bestHit, bestOnl] = bestDir(tip, dtip, A, distA, claimed, ...
    theta, bodyMask, distOut, gX, gY, maxLen, sweepRng, sweepStp, bodyMax, band, ...
    wShaft, wLen, prior, wPrior, fixedSeg, wOverlap, fixedBase, wBaseSep, minBaseSep, H, W)
x0 = tip(1); y0 = tip(2);
ix0 = min(max(round(x0),1),W); iy0 = min(max(round(y0),1),H);
bAng = atan2(-gY(iy0,ix0), -gX(iy0,ix0));
th0  = theta(iy0,ix0);
seeds = [th0, th0+pi, bAng];
cands = [];
for sdir = seeds, cands = [cands, sdir + deg2rad(-sweepRng:sweepStp:sweepRng)]; end %#ok<AGROW>
cands = cands(cos(cands - bAng) > cos(deg2rad(bodyMax)));
cands = unique(round(cands*1e3)/1e3);

best = -inf; bestAng = bAng; bestHit = []; bestOnl = false(size(A,1),1);
for ang = cands
    [osc, hit] = rayConsensus(x0,y0,ang,theta,bodyMask,maxLen,H,W);
    if isempty(hit), continue; end
    u = [cos(ang), sin(ang)];
    rel = A - [x0, y0];
    proj = rel(:,1)*u(1) + rel(:,2)*u(2);
    perp = abs(-rel(:,1)*u(2) + rel(:,2)*u(1));
    onl = (perp<band) & (proj>1) & (proj<hit(3)) & (distA(:)<dtip) & ~claimed(:);
    sc = osc + wShaft*nnz(onl) - wLen*(hit(3)/maxLen);
    if ~isempty(prior)
        dpr = abs(atan2(sin(ang-prior), cos(ang-prior)));
        sc = sc - wPrior*(dpr/pi);
    end
    for q = 1:numel(fixedSeg)
        S = fixedSeg{q};
        if segCross([x0 y0], hit(1:2), S(1,:), S(2,:)), sc = sc - wOverlap; break; end
    end
    % base-separation: penalize a base near an already-fixed base (avoid
    % several distinct filopodia emanating from one body point)
    if ~isempty(fixedBase) && wBaseSep > 0
        d2b = (fixedBase(:,1)-hit(1)).^2 + (fixedBase(:,2)-hit(2)).^2;
        if min(d2b) < minBaseSep^2, sc = sc - wBaseSep; end
    end
    if sc > best, best=sc; bestAng=ang; bestHit=hit; bestOnl=onl; end
end
end

% ===================================================================
function [osc, hit] = rayConsensus(x0,y0,ang,theta,bodyMask,maxLen,H,W)
oss=0; nO=0; hit=[];
for s = 1:1.5:maxLen
    x=x0+s*cos(ang); y=y0+s*sin(ang); ix=round(x); iy=round(y);
    if ix<1||ix>W||iy<1||iy>H, break; end
    if bodyMask(iy,ix), hit=[x y s]; break; end
    oss=oss+abs(cos(theta(iy,ix)-ang)); nO=nO+1;
end
if isempty(hit)||nO<3, osc=-inf; hit=[]; else, osc=oss/nO; end
end

% ===================================================================
function pr = neighborPrior(tip, fixedXY, fixedAng, R)
pr = [];
if isempty(fixedXY), return; end
d2 = (fixedXY(:,1)-tip(1)).^2 + (fixedXY(:,2)-tip(2)).^2;
w = exp(-d2/(2*R^2));
if sum(w) < 1e-3, return; end
pr = atan2(sum(w.*sin(fixedAng(:))), sum(w.*cos(fixedAng(:))));
end

% ===================================================================
function tf = segCross(p1,p2,p3,p4)
tf = (ccw(p1,p3,p4)~=ccw(p2,p3,p4)) && (ccw(p1,p2,p3)~=ccw(p1,p2,p4));
end
function t = ccw(a,b,c), t = (c(2)-a(2))*(b(1)-a(1)) > (b(2)-a(2))*(c(1)-a(1)); end

% ===================================================================
function z = zsc(x)
s = std(x); if s<eps, z = zeros(size(x)); else, z = (x-mean(x))/s; end
end
function v = gf(s,n,d), if isfield(s,n)&&~isempty(s.(n)), v=s.(n); else, v=d; end, end
