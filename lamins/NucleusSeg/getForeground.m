function BW = getForeground(BL, conn)
% This function finds the largest connected foreground object for a
% provided background/border
%
%   BL = label or binary matrix of the background
%   BW = mask of the largest connected foreground object

    dim = size(BL);
    BWf = zeros(dim);
    BWf(BL==0) = 1;
    CC = bwconncomp(BWf, conn);
    [dummy, maxIdx] = max(cellfun('length',CC.PixelIdxList));
    BW = zeros(dim); % create new mask
    BW(CC.PixelIdxList{maxIdx}) = 1; % only include biggest component
end