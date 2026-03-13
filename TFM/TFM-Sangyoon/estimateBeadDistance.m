function [beadDistAll, beadDensity, aggregationLevel, meanSegAreaUm2] = estimateBeadDistance(img, pixSize, sigma, outputFolder, h)
%function [beadDistAll,beadDensity,aggregationLevel,meanSegAreaUm2] = estimateBeadDistance(img,pixSize,sigma,outputFolder,h)
% Identifies beads and quantifies statistics about bead-to-bead distance,
% bead density, aggregation level, and segmentation area.
% Saves two overlay images into outputFolder (if provided):
%   - beadDetection.tif  : detected beads (red circles) on raw image
%   - beadSegmentation.tif : mask regions color-coded by area (um^2)
%
%   input:
%           img:            bead image
%           pixSize:        pixel size in nm (output from MovieData)
%           sigma:          point source sigma (from MD)
%           outputFolder:   path to save overlay images ('' or [] to skip)
%           h:              figure handle for live display (optional)
%
%   output:
%           beadDistAll:        closest distance (um) for each bead
%           beadDensity:        number of beads per area (#/um2)
%           aggregationLevel:   mean number of pstruct detections per
%                               mask segmentation region. 1 = non-aggregated,
%                               >1 = aggregated beads per region.
%           meanSegAreaUm2:     mean area of mask segmentation regions (um^2)
%
% May 2020, Sangyoon Han
% Updated: aggregationLevel, meanSegAreaUm2, and image saving added.
% Fixed: empty pstruct guard; correct 3-channel sub2ind indexing. 2026.

if nargin < 4, outputFolder = []; end
if nargin < 5, h = []; end

%% Settings
pixSizeUm = pixSize / 1000;   % nm -> um
maxDistUm = 1.5;
maxDist   = maxDistUm / pixSizeUm;

%% Detection
[pstruct, mask] = pointSourceDetection(img, sigma, 'FitMixtures', true, 'MaxMixtures', 10);

%% Nearest-neighbour distances
% Guard: need at least 2 detections for a meaningful nearest-neighbour query
if isempty(pstruct) || isempty(pstruct.x) || numel(pstruct.x) < 2
    beadDistAll = [];
else
    [~, dist] = KDTreeBallQuery([pstruct.x' pstruct.y'], [pstruct.x' pstruct.y'], maxDist);
    numDist     = cellfun(@length, dist);
    dist2       = dist(numDist > 1);
    distAll     = cellfun(@(x) x(2), dist2);
    beadDistAll = distAll * pixSizeUm;   % um
end

%% Connected-component analysis of mask
CC       = bwconncomp(mask);
nRegions = CC.NumObjects;
regionAreas_pix = cellfun(@numel, CC.PixelIdxList);
regionAreas_um2 = regionAreas_pix * pixSizeUm^2;

if nRegions > 0 && ~isempty(pstruct) && ~isempty(pstruct.x)
    labelImg         = labelmatrix(CC);
    beadRegionLabels = labelImg(sub2ind(size(img), pstruct.y_init, pstruct.x_init));

    beadsPerRegion = zeros(1, nRegions);
    for r = 1:nRegions
        beadsPerRegion(r) = sum(beadRegionLabels == r);
    end

    occupiedRegions  = beadsPerRegion(beadsPerRegion > 0);
    aggregationLevel = mean(occupiedRegions);
    meanSegAreaUm2   = mean(regionAreas_um2);
else
    aggregationLevel = NaN;
    meanSegAreaUm2   = NaN;
end

%% Bead density
curAreaPix  = (size(img,1) - ceil(4*sigma)*2) * (size(img,2) - ceil(4*sigma)*2);
curArea     = curAreaPix * pixSizeUm^2;
if ~isempty(pstruct) && ~isempty(pstruct.x)
    beadDensity = numel(pstruct.x) / curArea;
else
    beadDensity = 0;
end

%% Save overlay images
if ~isempty(outputFolder) && (ischar(outputFolder) || isstring(outputFolder))
    % Ensure the folder exists (create if needed)
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    % --- Normalise image to uint8 for display ---
    imgNorm = double(img);
    imgNorm = imgNorm - min(imgNorm(:));
    denom   = max(imgNorm(:));
    if denom == 0, denom = 1; end
    imgNorm = uint8(255 * imgNorm / denom);
    imgRGB  = repmat(imgNorm, [1 1 3]);   % grayscale RGB base

    % ---- 1) Bead detection overlay ----------------------------------------
    % Draw a circle of radius 2*sigma around each detected bead in red.
    % Uses correct 3-channel linear indexing into the RGB array.
    detImg = imgRGB;
    if ~isempty(pstruct) && ~isempty(pstruct.x)
        theta  = linspace(0, 2*pi, 60);
        r_circ = 2 * sigma;
        imgH   = size(detImg, 1);
        imgW   = size(detImg, 2);
        for k = 1:numel(pstruct.x)
            cx = round(pstruct.x(k) + r_circ * cos(theta));
            cy = round(pstruct.y(k) + r_circ * sin(theta));
            valid = cx >= 1 & cx <= imgW & cy >= 1 & cy <= imgH;
            cy_v  = cy(valid);
            cx_v  = cx(valid);
            if isempty(cy_v), continue; end
            % Correct 3-channel linear indexing
            nv   = numel(cy_v);
            linR = sub2ind([imgH imgW 3], cy_v, cx_v, ones(1, nv, 'int32'));
            linG = sub2ind([imgH imgW 3], cy_v, cx_v, 2*ones(1, nv, 'int32'));
            linB = sub2ind([imgH imgW 3], cy_v, cx_v, 3*ones(1, nv, 'int32'));
            detImg(linR) = 255;
            detImg(linG) = 0;
            detImg(linB) = 0;
        end
    end
    imwrite(detImg, fullfile(outputFolder, 'beadDetection.tif'));

    % ---- 2) Colorised segmentation overlay --------------------------------
    % Each connected component is colored by its area (um^2) using 'jet'.
    % Background remains as the grayscale image.
    segImg = imgRGB;
    if nRegions > 0
        cmap       = jet(256);
        areaMin    = min(regionAreas_um2);
        areaMax    = max(regionAreas_um2);
        areaRange  = max(areaMax - areaMin, eps);
        imgH       = size(segImg, 1);
        imgW       = size(segImg, 2);

        for r = 1:nRegions
            cmapIdx = round(1 + 255 * (regionAreas_um2(r) - areaMin) / areaRange);
            color   = uint8(cmap(cmapIdx, :) * 255);   % [R G B]

            pixIdx      = CC.PixelIdxList{r};   % linear indices into 2D image
            [pyy, pxx]  = ind2sub([imgH imgW], pixIdx);
            np          = numel(pyy);

            linR = sub2ind([imgH imgW 3], pyy, pxx, ones(np, 1, 'int32'));
            linG = sub2ind([imgH imgW 3], pyy, pxx, 2*ones(np, 1, 'int32'));
            linB = sub2ind([imgH imgW 3], pyy, pxx, 3*ones(np, 1, 'int32'));

            % Blend 60% jet color + 40% original gray for readability
            segImg(linR) = uint8(0.6*double(color(1)) + 0.4*double(segImg(linR)));
            segImg(linG) = uint8(0.6*double(color(2)) + 0.4*double(segImg(linG)));
            segImg(linB) = uint8(0.6*double(color(3)) + 0.4*double(segImg(linB)));
        end

        % Colorbar strip on the right edge (~2.5% image width)
        barW   = max(10, round(imgW * 0.025));
        barImg = zeros(imgH, barW, 3, 'uint8');
        for row = 1:imgH
            ci = round(1 + 255*(1 - (row-1)/max(imgH-1, 1)));   % top = max area
            barImg(row, :, 1) = uint8(cmap(ci, 1) * 255);
            barImg(row, :, 2) = uint8(cmap(ci, 2) * 255);
            barImg(row, :, 3) = uint8(cmap(ci, 3) * 255);
        end
        segImg = [segImg, barImg];
    end
    imwrite(segImg, fullfile(outputFolder, 'beadSegmentation.tif'));
end

%% Optional live figure
if ~isempty(h) && ~isempty(pstruct) && ~isempty(pstruct.x)
    figure(h);
    imshow(img, []), hold on
    plot(pstruct.x, pstruct.y, 'ro')
end

%% Display summary
disp(['Mean distance=       ' num2str(mean(beadDistAll),   3) ' um'])
disp(['Bead density=        ' num2str(beadDensity,         3) ' /um^2'])
disp(['Aggregation level=   ' num2str(aggregationLevel,    3) ' beads/region'])
disp(['Mean segment area=   ' num2str(meanSegAreaUm2,      3) ' um^2'])
end