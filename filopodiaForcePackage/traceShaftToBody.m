function path = traceShaftToBody(x0, y0, theta, bodyMask, distOut, maxLen, stepLen)
%TRACESHAFTTOBODY  Trace a filopodium shaft from a tip to the cell body.
%
%   path = traceShaftToBody(x0, y0, theta, bodyMask, distOut, maxLen, stepLen)
%
% Starting at the tip (x0,y0), step along the local steerable ridge
% orientation theta toward the body (decreasing distance-to-body), producing
% a near-straight shaft polyline that ends on the body edge. Returns an Mx2
% [x y] polyline; the last row is the base (on the body) if the body was
% reached. If the body is not reached within maxLen steps, the trace is
% returned as-is (caller can reject tips whose trace never reaches the body).
% Sangyoon J. Han / 2026

if nargin < 6 || isempty(maxLen),  maxLen  = 160; end
if nargin < 7 || isempty(stepLen), stepLen = 1.5; end
[H, W] = size(bodyMask);
cx = x0; cy = y0;
path = [x0, y0];

for it = 1:maxLen
    iy = round(cy); ix = round(cx);
    if iy < 1 || iy > H || ix < 1 || ix > W, break; end
    if bodyMask(iy, ix), break; end          % reached body
    th = theta(iy, ix);
    d0 = distOut(iy, ix);
    best = []; bd = d0;
    % follow the ridge orientation (and its perpendicular) toward the body
    for ang = [th, th+pi, th+pi/2, th-pi/2]
        ny = cy + stepLen*sin(ang);
        nx = cx + stepLen*cos(ang);
        jy = round(ny); jx = round(nx);
        if jy >= 1 && jy <= H && jx >= 1 && jx <= W && distOut(jy,jx) < bd
            bd = distOut(jy,jx);
            best = [nx, ny];
        end
    end
    if isempty(best), break; end             % stuck (no bodyward step)
    cx = best(1); cy = best(2);
    path(end+1, :) = [cx, cy]; %#ok<AGROW>
end
end