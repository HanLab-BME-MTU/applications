function [centerline, arc, len, ok] = traceFilopodiaShaft(res, theta, tipPos, basePos, p)
%TRACEFILOPODIASHAFT  Orientation-guided minimum-cost trace of a dim shaft.
%
%   [centerline, arc, len, ok] = traceFilopodiaShaft(res, theta, tipPos, basePos, p)
%
% Both endpoints are known (talin tip punctum and base focal adhesion), so the
% shaft is recovered as the minimum-cost path between them on an 8-connected
% pixel graph. Node-to-node cost favors high steerable response (cost ~ 1/res)
% and agreement between the step direction and the local ridge orientation
% theta. Because the endpoints anchor the path, dim or broken shafts still
% resolve; with no signal the path falls back to a near-geodesic line, a sane
% prior for a filopodium.
%
% INPUT
%   res      : HxW steerable ridge response (from multiscaleSteerableDetector)
%   theta    : HxW ridge orientation in radians (undirected, mod pi)
%   tipPos   : [x y] sub-pixel tip
%   basePos  : [x y] sub-pixel base
%   p        : struct with fields
%                .PathCostFloor    floor on res to avoid div-by-zero (e.g 1e-3)
%                .OrientLambda     orientation penalty weight (e.g. 2)
%                .OrientTolerance  deg; mismatch beyond this is fully penalized
%                .MaxTipBaseDist   px; reject if straight-line distance exceeds
%
% OUTPUT
%   centerline : Nx2 [x y], ordered base -> tip
%   arc        : Nx1 cumulative arclength (px), arc(1)=0
%   len        : total length (px) = arc(end)
%   ok         : true if a valid path was found
%
% Sangyoon J. Han / 2026

centerline = []; arc = []; len = NaN; ok = false;
if nargin < 5, p = struct(); end
if ~isfield(p, 'PathCostFloor'),   p.PathCostFloor   = 1e-3; end
if ~isfield(p, 'OrientLambda'),    p.OrientLambda    = 2;    end
if ~isfield(p, 'OrientTolerance'), p.OrientTolerance = 30;   end
if ~isfield(p, 'MaxTipBaseDist'),  p.MaxTipBaseDist  = Inf;  end

[H, W] = size(res);
bx = round(basePos(1)); by = round(basePos(2));
tx = round(tipPos(1));  ty = round(tipPos(2));
if any([bx by tx ty] < 1) || bx > W || tx > W || by > H || ty > H, return; end
if hypot(tx - bx, ty - by) > p.MaxTipBaseDist, return; end

% --- restrict to a bounding box around the two endpoints (+ margin) ---
margin = max(10, round(0.4 * hypot(tx - bx, ty - by)));
r0 = max(1, min(by, ty) - margin); r1 = min(H, max(by, ty) + margin);
c0 = max(1, min(bx, tx) - margin); c1 = min(W, max(bx, tx) + margin);
R = res(r0:r1, c0:c1);
T = theta(r0:r1, c0:c1);
[h, w] = size(R);
nodeCost = 1 ./ max(R, p.PathCostFloor);   % low where response is high

% endpoints in local coords
bl = [by - r0 + 1, bx - c0 + 1];   % [row col]
tl = [ty - r0 + 1, tx - c0 + 1];
sInd = sub2ind([h w], bl(1), bl(2));
eInd = sub2ind([h w], tl(1), tl(2));

% --- build 8-connected sparse graph ---
[rr, cc] = ndgrid(1:h, 1:w);
idx = (1:h*w)';
nb = [-1 -1; -1 0; -1 1; 0 -1; 0 1; 1 -1; 1 0; 1 1];
src = []; dst = []; wgt = [];
tolRad = deg2rad(p.OrientTolerance);
for k = 1:size(nb, 1)
    dr = nb(k, 1); dc = nb(k, 2);
    nr = rr + dr; ncl = cc + dc;
    valid = nr >= 1 & nr <= h & ncl >= 1 & ncl <= w;
    a = idx(valid);
    b = sub2ind([h w], nr(valid), ncl(valid));
    stepLen = hypot(dr, dc);
    stepAng = atan2(dr, dc);                 % direction of this edge
    % orientation mismatch vs local ridge (undirected -> fold to [0 pi/2])
    dAng = abs(angleDiffMod(stepAng, T(a)));
    pen = min(dAng / max(tolRad, eps), 1);   % 0 (aligned) .. 1 (orthogonal)
    base = 0.5 * (nodeCost(a) + nodeCost(b)) * stepLen;
    src = [src; a];                          %#ok<AGROW>
    dst = [dst; b];                          %#ok<AGROW>
    wgt = [wgt; base .* (1 + p.OrientLambda * pen)]; %#ok<AGROW>
end
G = graph(src, dst, wgt, h * w);

% --- shortest path ---
[pathNodes, d] = shortestpath(G, sInd, eInd);
if isempty(pathNodes) || ~isfinite(d), return; end

[pr, pc] = ind2sub([h w], pathNodes(:));
xy = [pc + c0 - 1, pr + r0 - 1];             % back to global [x y]

% snap exact sub-pixel endpoints
xy(1, :)   = basePos(:)';
xy(end, :) = tipPos(:)';

% arclength
steps = sqrt(sum(diff(xy, 1, 1).^2, 2));
arc = [0; cumsum(steps)];
centerline = xy;
len = arc(end);
ok = true;
end

% ----- helper: smallest angle between an edge dir and an undirected ridge -----
function d = angleDiffMod(edgeAng, ridgeAng)
% ridgeAng is undirected (period pi). Return value in [0, pi/2].
d = mod(edgeAng - ridgeAng, pi);     % [0, pi)
d = min(d, pi - d);                  % fold to [0, pi/2]
end
