function [ pairdist, r, c ] = histpdist( nms , thresh )
%histpdist Display histogram of pairwise distances
    if(nargin < 2)
        thresh = 0.5;
    end
    bps = bwmorph(bwmorph(nms > thresh,'skel'),'branchpoints');
    [r,c] = find(bps);
    pairdist = pdist([r c]);
    optimalHistogram(pairdist);
end

