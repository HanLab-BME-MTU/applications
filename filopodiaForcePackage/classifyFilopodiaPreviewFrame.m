function filo = classifyFilopodiaPreviewFrame(tipXY, bodyMask, p)
%CLASSIFYFILOPODIAPREVIEWFRAME  Single-frame filopodia shaft tracing for the
%P4 settings preview. For each tip position, shoot a ray body-ward and find
%the base where it enters the cell body, exactly like Pass B of
%classifyMovieFilopodia (same MaxShaftLen / BodyMaxAngle gating). This is an
%approximation of the full classifier: the per-track direction that the real
%Pass A assigns from the whole movie is here approximated by the body-ward
%direction at the tip, which is enough to tune the shaft-tracing parameters
%(MaxShaftLen, BodyMaxAngle, ShaftBand) interactively.
%
%   filo = classifyFilopodiaPreviewFrame(tipXY, bodyMask, p)
%
% tipXY : nTip x 2 [x y] candidate tip positions in this frame
% bodyMask : logical cell-body mask for this frame
% p : funParams (MaxShaftLen, BodyMaxAngle, MinTipDist, ShaftBand)
%
% Returns filo struct array with fields tipPos, basePos, L (px), accepted.
% Sangyoon J. Han / 2026

filo = struct('tipPos',{},'basePos',{},'L',{},'accepted',{});
if isempty(tipXY) || isempty(bodyMask) || ~any(bodyMask(:)), return; end

[H,W] = size(bodyMask);
maxLen   = gf(p,'MaxShaftLen',160);
bodyMaxA = gf(p,'BodyMaxAngle',75) * pi/180;
minTipD  = gf(p,'MinTipDist',6);

dOut = bwdist(bodyMask);          % distance outside body (0 inside)
% body-ward direction field: gradient of the inside-distance transform points
% inward; outside, the direction toward nearest body pixel is -grad(dOut).
[gy, gx] = gradient(dOut);

cen = regionprops(bodyMask,'Centroid');
if ~isempty(cen), bc = cen(1).Centroid; else, bc = [W/2 H/2]; end

n = 0;
for i = 1:size(tipXY,1)
    x0 = tipXY(i,1); y0 = tipXY(i,2);
    ix = min(max(round(x0),1),W); iy = min(max(round(y0),1),H);
    if bodyMask(iy,ix), continue; end             % tip must be outside body
    if dOut(iy,ix) < minTipD, continue; end       % and at least MinTipDist out

    % body-ward direction: prefer -gradient(dOut) (toward nearest body edge);
    % fall back to direction toward body centroid.
    gdir = -[gx(iy,ix), gy(iy,ix)];
    if norm(gdir) < 1e-6
        gdir = [bc(1)-x0, bc(2)-y0];
    end
    ang0 = atan2(gdir(2), gdir(1));

    % small sweep around body-ward to find the shortest ray that hits body
    best = []; bestL = inf; bestAng = ang0;
    for da = -bodyMaxA:(5*pi/180):bodyMaxA
        ang = ang0 + da;
        for s = 1:1.5:maxLen
            x = x0 + s*cos(ang); y = y0 + s*sin(ang);
            xi = round(x); yi = round(y);
            if xi<1||xi>W||yi<1||yi>H, break; end
            if bodyMask(yi,xi)
                if s < bestL, bestL = s; best = [x y]; bestAng = ang; end
                break;
            end
        end
    end

    accepted = ~isempty(best);
    n = n + 1;
    filo(n).tipPos   = [x0 y0];
    if accepted
        filo(n).basePos = best; filo(n).L = bestL;
    else
        filo(n).basePos = [NaN NaN]; filo(n).L = NaN;
    end
    filo(n).accepted = accepted; %#ok<STRNU>
end
end

function v = gf(s,n,d), if isfield(s,n)&&~isempty(s.(n)), v=s.(n); else, v=d; end, end
