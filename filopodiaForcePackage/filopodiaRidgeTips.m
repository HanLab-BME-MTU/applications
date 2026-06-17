function tips = filopodiaRidgeTips(bodyMask, shaftMask, minRidgeArea, minBranchLen, tipMinReach, gapBridge, mergeRadius)
%FILOPODIARIDGETIPS  Distal endpoints of the body-connected steerable ridge.
%
%   tips = filopodiaRidgeTips(bodyMask, shaftMask, minRidgeArea, minBranchLen, ...
%                             tipMinReach, gapBridge, mergeRadius)
%
% Returns Nx2 [x y] tip candidates: distal ends of every branch of the
% body-connected shaftMask ridge (multiple tips per component, one per
% filopodium branch). Free-floating debris islands are ignored.
%
%   gapBridge   : px the body is dilated before the connectivity test, so a
%                 filopodium whose ridge breaks within this gap of the body
%                 edge is still kept (recovers near-root disconnections).
%   mergeRadius : px; tips closer than this are merged (one kept) to avoid
%                 clusters of endpoints from a single blobby ridge.
%
% Per-frame detection is intentionally permissive; transient false tips are
% removed downstream by track persistence (P4 MinTipLifetime).
% Shared by P2 (detection augmentation) and P4 (classification).
% Sangyoon J. Han / 2026

if nargin < 6 || isempty(gapBridge),   gapBridge = 2;   end
if nargin < 7 || isempty(mergeRadius), mergeRadius = 0; end
[H, W] = size(bodyMask);
tips = zeros(0, 2);

ridge = shaftMask & ~bodyMask;
ridge = bwareaopen(ridge, minRidgeArea);
if ~any(ridge(:)), return; end

% keep components touching the body (dilated by gapBridge to bridge root gaps);
% drop free-floating debris islands
seedBody = imdilate(bodyMask, strel('disk', max(1, round(gapBridge))));
ridge = imreconstruct(ridge & seedBody, ridge);
if ~any(ridge(:)), return; end

% pruned skeleton: bwskel removes spur branches shorter than minBranchLen,
% so only real filopodia branches keep an endpoint
try
    sk = bwskel(logical(ridge), 'MinBranchLength', round(minBranchLen));
catch
    sk = bwmorph(ridge, 'thin', Inf);
    sk = bwmorph(sk, 'spur', round(minBranchLen));
end

ep = bwmorph(sk, 'endpoints');
if ~any(ep(:)), return; end

distOut = bwdist(bodyMask);
[ye, xe] = find(ep);
keep = distOut(sub2ind([H W], ye, xe)) > tipMinReach;
xe = xe(keep); ye = ye(keep);
if isempty(xe), return; end

% merge nearby endpoints (greedy: keep the more distal one first)
if mergeRadius > 0 && numel(xe) > 1
    dd = distOut(sub2ind([H W], ye, xe));
    [~, ord] = sort(dd, 'descend');
    xe = xe(ord); ye = ye(ord);
    keepM = true(numel(xe),1);
    for i = 1:numel(xe)
        if ~keepM(i), continue; end
        dup = (xe - xe(i)).^2 + (ye - ye(i)).^2 < mergeRadius^2;
        dup(1:i) = false;
        keepM(dup) = false;
    end
    xe = xe(keepM); ye = ye(keepM);
end

tips = [xe, ye];
end