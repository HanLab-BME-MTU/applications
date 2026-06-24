function [tipXY, baseXY, allXY] = detectFilopodiaPointsPreview(img, bodyMask, p)
%DETECTFILOPODIAPOINTSPREVIEW  Light per-frame point detection for the P2
%settings preview. Runs pointSourceDetection on the talin frame, then splits
%puncta into tips (outside body, within TipMaxDistFromBody) and bases (within
%BaseSearchBand of the body edge). This mirrors the gating in
%detectMovieFilopodia's 'all' mode closely enough for parameter tuning; the
%batch run still does the full ridge-aware detection + shaft tracing.
%
%   [tipXY, baseXY, allXY] = detectFilopodiaPointsPreview(img, bodyMask, p)
% Sangyoon J. Han / 2026

sigma = gf(p,'PSFsigma',2.1);
alpha = gf(p,'Alpha',0.05);
args = {'Alpha',alpha};
cr = gf(p,'ConfRadius',[]); if ~isempty(cr), args=[args {'ConfRadius',cr}]; end
ws = gf(p,'WindowSize',[]);  if ~isempty(ws), args=[args {'WindowSize',ws}]; end

tipXY=zeros(0,2); baseXY=zeros(0,2); allXY=zeros(0,2);
try
    ps = pointSourceDetection(img, sigma, args{:});
catch
    return;
end
if isempty(ps) || ~isfield(ps,'x') || isempty(ps.x), return; end
allXY = [ps.x(:), ps.y(:)];

[H,W] = size(img);
if isempty(bodyMask), bodyMask = false(H,W); end
% distance outside the body (0 inside)
dOut = bwdist(bodyMask);
% distance inside from the edge (for base band)
dIn  = bwdist(~bodyMask);

tipMax   = gf(p,'TipMaxDistFromBody',130);
baseBand = gf(p,'BaseSearchBand',5);
inBand   = gf(p,'BaseInsideBand',4);

for i = 1:size(allXY,1)
    x = min(max(round(allXY(i,1)),1),W);
    y = min(max(round(allXY(i,2)),1),H);
    inside = bodyMask(y,x);
    if ~inside && dOut(y,x) <= tipMax
        tipXY(end+1,:) = allXY(i,:); %#ok<AGROW>
    end
    % base: just outside within band, or just inside within inBand
    if (~inside && dOut(y,x) <= baseBand) || (inside && dIn(y,x) <= inBand)
        baseXY(end+1,:) = allXY(i,:); %#ok<AGROW>
    end
end
end

function v = gf(s,n,d), if isfield(s,n)&&~isempty(s.(n)), v=s.(n); else, v=d; end, end