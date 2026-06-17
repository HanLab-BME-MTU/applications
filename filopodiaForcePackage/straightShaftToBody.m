function [shaft, base, ang, score, reached] = straightShaftToBody(x0, y0, allPts, theta, bodyMask, distOut, p)
%STRAIGHTSHAFTTOBODY  Best STRAIGHT filopodium shaft from a tip to the body.
%
%   [shaft, base, ang, score, reached] = ...
%       straightShaftToBody(x0, y0, allPts, theta, bodyMask, distOut, p)
%
% Filopodia are straight (formin-driven linear actin polymerization), so the
% shaft is modeled as a straight segment from the tip (x0,y0) to the body.
% The direction is chosen robustly (the per-pixel steerable orientation is
% noisy) by scoring candidate directions around the tip's ridge orientation
% and the body-ward direction by:
%   (a) orientation consensus : mean |cos(theta_along_ray - phi)|
%   (b) collinear support     : # of other tip adhesions (allPts) lying along
%                               the ray within a perpendicular band
%   (c) length penalty        : shorter (more radial) shafts preferred
% Only directions within BodyMaxAngle of the body-ward direction are allowed,
% which removes shafts that graze tangentially along the cell edge.
%
% Returns shaft as a densified Mx2 [x y] straight polyline (tip..base), the
% base on the body, the chosen angle, the score, and whether the body was
% reached. p fields (with defaults): MaxShaftLen, SweepRange, SweepStep,
% BodyMaxAngle, AlignBand, AlignWeight, LenPenalty.
% Sangyoon J. Han / 2026

maxLen   = gf(p,'MaxShaftLen',160);
sweepRng = gf(p,'SweepRange',35);     % deg around each seed direction
sweepStp = gf(p,'SweepStep',3);       % deg
bodyMax  = gf(p,'BodyMaxAngle',60);   % deg; shaft must be within this of body-ward
band     = gf(p,'AlignBand',3.5);     % px; perpendicular band for collinear support
wAlign   = gf(p,'AlignWeight',0.12);
wLen     = gf(p,'LenPenalty',0.6);
[H, W] = size(bodyMask);

% body-ward direction at the tip (toward decreasing distance-to-body)
[gX, gY] = gradient(distOut);
ix0 = min(max(round(x0),1),W); iy0 = min(max(round(y0),1),H);
bAng = atan2(-gY(iy0,ix0), -gX(iy0,ix0));
th0  = theta(iy0, ix0);

% candidate directions: sweeps around ridge orientation (both senses) and body-ward
seeds = [th0, th0+pi, bAng];
cands = [];
for sdir = seeds
    cands = [cands, sdir + deg2rad(-sweepRng:sweepStp:sweepRng)]; %#ok<AGROW>
end
% keep only directions within BodyMaxAngle of body-ward
cands = cands(cos(cands - bAng) > cos(deg2rad(bodyMax)));
cands = unique(round(cands*1e3)/1e3);

best = -inf; bestHit = []; bestAng = bAng;
for ang = cands
    [sc, hit] = scoreDir(x0, y0, ang, allPts, theta, bodyMask, distOut, ...
                         maxLen, band, wAlign, wLen, H, W);
    if sc > best, best = sc; bestHit = hit; bestAng = ang; end
end

ang = bestAng; score = best;
if isempty(bestHit)
    reached = false;
    shaft = [x0, y0]; base = [x0, y0];
    return;
end
reached = true;
base = bestHit(1:2);
L = bestHit(3);
nstep = max(2, round(L));
ss = linspace(0, L, nstep)';
shaft = [x0 + ss*cos(ang), y0 + ss*sin(ang)];
end

% -------------------------------------------------------------------
function [sc, hit] = scoreDir(x0,y0,ang,allPts,theta,bodyMask,distOut,maxLen,band,wAlign,wLen,H,W)
oss = 0; nO = 0; hit = [];
for s = 1:1.5:maxLen
    x = x0 + s*cos(ang); y = y0 + s*sin(ang);
    ix = round(x); iy = round(y);
    if ix<1 || ix>W || iy<1 || iy>H, break; end
    if bodyMask(iy,ix), hit = [x, y, s]; break; end
    oss = oss + abs(cos(theta(iy,ix) - ang)); nO = nO + 1;
end
if isempty(hit) || nO < 3, sc = -inf; return; end
osc = oss / nO;
% collinear support: other tips along the ray within perpendicular band
aligned = 0;
if ~isempty(allPts)
    dx = cos(ang); dy = sin(ang);
    rx = allPts(:,1)-x0; ry = allPts(:,2)-y0;
    proj = rx*dx + ry*dy; perp = abs(-rx*dy + ry*dx);
    aligned = sum(proj>3 & proj<hit(3) & perp<band);
end
sc = osc + wAlign*aligned - wLen*(hit(3)/maxLen);
end

% -------------------------------------------------------------------
function v = gf(s,n,d), if isfield(s,n)&&~isempty(s.(n)), v=s.(n); else, v=d; end, end
