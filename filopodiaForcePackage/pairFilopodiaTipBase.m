function baseForTip = pairFilopodiaTipBase(tipPos, basePos, maxDist)
%PAIRFILOPODIATIPBASE  Global one-to-one tip<->base assignment via LAP.
%
%   baseForTip = pairFilopodiaTipBase(tipPos, basePos, maxDist)
%
% Uses uTrack's lap.m. Each tip is matched to at most one base (and vice
% versa); pairs farther apart than maxDist are forbidden, so a tip with no
% admissible base stays unmatched. This resolves the ambiguity that nearest-
% base assignment creates for crossing or densely packed filopodia.
%
% INPUT
%   tipPos   : Nt x 2 [x y]
%   basePos  : Nb x 2 [x y]
%   maxDist  : px, maximum admissible tip-base distance
%
% OUTPUT
%   baseForTip : Nt x 1, base index for each tip (NaN if unmatched)
%
% Sangyoon J. Han / 2026

Nt = size(tipPos, 1); Nb = size(basePos, 1);
baseForTip = nan(Nt, 1);
if Nt == 0 || Nb == 0, return; end

% Euclidean cost, gated by maxDist
D = createDistanceMatrix(tipPos, basePos);   % Nt x Nb (uTrack helper)
NONLINK = -1;
cc = D;
cc(cc > maxDist) = NONLINK;                  % forbid links beyond range

% augmentCC=1 lets lap add births/deaths so unequal/ungated rows are allowed
[link12, ~] = lap(cc, NONLINK, [], 1);
link12 = double(link12(:));

for it = 1:Nt
    j = link12(it);
    if j >= 1 && j <= Nb
        baseForTip(it) = j;                  % real tip->base link
    end                                       % else linked to a "death" -> NaN
end
end
