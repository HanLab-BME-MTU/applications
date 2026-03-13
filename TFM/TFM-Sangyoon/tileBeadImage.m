function [tiles, tileOrigins, validMask] = tileBeadImage(img, tileSize, minMeanFrac)
% tileBeadImage  Divide an image into non-overlapping ROI tiles and flag
%                dark/empty tiles so they can be skipped during analysis.
%
%   [tiles, tileOrigins, validMask] = tileBeadImage(img, tileSize, minMeanFrac)
%
%   Inputs
%     img          - 2-D grayscale image (any numeric class)
%     tileSize     - scalar tile side length in pixels (default 500)
%     minMeanFrac  - a tile is 'valid' if its mean intensity >= this fraction
%                   of the whole-image mean (default 0.3).
%                   Set to 0 to keep all tiles.
%
%   Outputs
%     tiles        - cell array of cropped image patches
%     tileOrigins  - Nx2 array of [row col] top-left corners (1-based)
%     validMask    - logical Nx1 vector; true = tile passes intensity check
%
% Tiles that would extend beyond the image border are silently dropped
% (only full tileSize x tileSize patches are returned).
%
% Sangyoon Han, 2026.

if nargin < 2 || isempty(tileSize),    tileSize    = 500; end
if nargin < 3 || isempty(minMeanFrac), minMeanFrac = 0.3; end

[imgH, imgW] = size(img);
imgMean      = mean(double(img(:)));
threshold    = minMeanFrac * imgMean;

% Number of full tiles that fit
nTilesR = floor(imgH / tileSize);
nTilesC = floor(imgW / tileSize);

nTotal     = nTilesR * nTilesC;
tiles      = cell(nTotal, 1);
tileOrigins= zeros(nTotal, 2);
validMask  = false(nTotal, 1);

idx = 0;
for r = 1:nTilesR
    rowStart = (r-1)*tileSize + 1;
    rowEnd   = rowStart + tileSize - 1;
    for c = 1:nTilesC
        colStart = (c-1)*tileSize + 1;
        colEnd   = colStart + tileSize - 1;

        idx = idx + 1;
        patch = img(rowStart:rowEnd, colStart:colEnd);
        tiles{idx}       = patch;
        tileOrigins(idx,:) = [rowStart, colStart];
        validMask(idx)   = mean(double(patch(:))) >= threshold;
    end
end
end