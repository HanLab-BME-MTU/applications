function [centerline, arc, len, basePos, ok] = traceFilopodiaToBody(res, theta, bodyMask, tipPos, p)
%TRACEFILOPODIATOBODY  Trace a filopodium from its tip to the cell body.
%
%   [centerline, arc, len, basePos, ok] = traceFilopodiaToBody(res, theta, bodyMask, tipPos, p)
%
% Starting from a known talin tip punctum (outside the body), the shaft is
% recovered as the minimum-cost path on an 8-connected pixel graph
% (cost ~ 1/res, with an orientation penalty against the local ridge theta)
% to the *nearest-in-cost* body-mask boundary pixel. The base is therefore
% defined geometrically as where the trace first reaches the body edge, not
% as a separate punctum, so bases land on the body boundary by construction
% and shafts never run into the body interior. Tips whose trace cannot reach
% the body within MaxTipBaseDist are rejected (debris / spurious puncta).
%
% A virtual sink connected to all in-box body-boundary pixels lets one
% shortest-path call find the least-cost contact point in a single pass.
%
% INPUT
%   res       : HxW steerable ridge response
%   theta     : HxW ridge orientation (rad, undirected, period pi)
%   bodyMask  : HxW logical cell-body mask
%   tipPos    : [x y] sub-pixel tip
%   p         : struct with .PathCostFloor .OrientLambda .OrientTolerance
%               .MaxTipBaseDist
%
% OUTPUT
%   centerline : Nx2 [x y], ordered base -> tip
%   arc        : Nx1 cumulative arclength (px), arc(1)=0 at base
%   len        : total length (px)
%   basePos    : [x y] body-edge contact point
%   ok         : true if a valid trace to the body was found
%
% Sangyoon J. Han / 2026

centerline = []; arc = []; len = NaN; basePos = []; ok = false;
if ~isfield(p,'PathCostFloor'),   p.PathCostFloor   = 1e-3; end
if ~isfield(p,'OrientLambda'),    p.OrientLambda    = 2;    end
if ~isfield(p,'OrientTolerance'), p.OrientTolerance = 30;   end
if ~isfield(p,'MaxTipBaseDist'),  p.MaxTipBaseDist  = 100;  end

[H, W] = size(res);
tx = round(tipPos(1)); ty = round(tipPos(2));
if tx < 1 || tx > W || ty < 1 || ty > H, return; end
if bodyMask(ty, tx), return; end            % tip must be outside the body

% --- bounding box around the tip, large enough to reach the body edge ---
R  = ceil(p.MaxTipBaseDist) + 5;
r0 = max(1, ty - R); r1 = min(H, ty + R);
c0 = max(1, tx - R); c1 = min(W, tx + R);
Rmap = res(r0:r1, c0:c1);
Tmap = theta(r0:r1, c0:c1);
perim = bwperim(bodyMask);                  % body boundary on full image
Pcrop = perim(r0:r1, c0:c1);
[h, w] = size(Rmap);
if ~any(Pcrop(:)), return; end              % no body edge within reach

nodeCost = 1 ./ max(Rmap, p.PathCostFloor);

% --- 8-connected grid graph (each undirected edge added once) ---
[rr, cc] = ndgrid(1:h, 1:w);
idx = (1:h*w)';
nb  = [-1 -1; -1 0; -1 1; 0 -1; 0 1; 1 -1; 1 0; 1 1];
src = []; dst = []; wgt = [];
tolRad = deg2rad(p.OrientTolerance);
for k = 1:size(nb, 1)
    dr = nb(k, 1); dc = nb(k, 2);
    nr = rr + dr; ncl = cc + dc;
    valid = nr >= 1 & nr <= h & ncl >= 1 & ncl <= w;
    a = idx(valid);
    b = sub2ind([h w], nr(valid), ncl(valid));
    kp = b > a; a = a(kp); b = b(kp);
    if isempty(a), continue; end
    stepLen = hypot(dr, dc);
    dAng = angleDiffMod(atan2(dr, dc), Tmap(a));
    pen  = min(dAng / max(tolRad, eps), 1);
    base = 0.5 * (nodeCost(a) + nodeCost(b)) * stepLen;
    src = [src; a]; dst = [dst; b]; wgt = [wgt; base .* (1 + p.OrientLambda * pen)]; %#ok<AGROW>
end

% --- virtual sink connected (0 cost) to every body-boundary pixel ---
sink = h * w + 1;
perimNodes = find(Pcrop);
src = [src; perimNodes];
dst = [dst; sink * ones(numel(perimNodes), 1)];
wgt = [wgt; zeros(numel(perimNodes), 1)];
G = graph(src, dst, wgt, sink);

sInd = sub2ind([h w], ty - r0 + 1, tx - c0 + 1);
[pathNodes, d] = shortestpath(G, sInd, sink);
if isempty(pathNodes) || ~isfinite(d), return; end
pathNodes(end) = [];                        % drop the sink

[pr, pc] = ind2sub([h w], pathNodes(:));
xy = [pc + c0 - 1, pr + r0 - 1];            % global [x y], tip -> base
xy(1, :) = tipPos(:)';                      % exact sub-pixel tip
xy = flipud(xy);                            % reorder base -> tip

basePos = xy(1, :);
steps = sqrt(sum(diff(xy, 1, 1).^2, 2));
arc = [0; cumsum(steps)];
len = arc(end);
if len > p.MaxTipBaseDist, return; end      % too long -> reject
centerline = xy;
ok = true;
end

% ----- smallest angle between an edge direction and an undirected ridge -----
function d = angleDiffMod(edgeAng, ridgeAng)
d = mod(edgeAng - ridgeAng, pi);
d = min(d, pi - d);                          % fold to [0, pi/2]
end
