function [beadDistAll, beadDensity, aggregationLevel, meanSegAreaUm2] = estimateBeadDistance(img, pixSize, sigma, outputFolder, h)
%function [beadDistAll,beadDensity,aggregationLevel,meanSegAreaUm2] = estimateBeadDistance(img,pixSize,sigma,outputFolder,h)
% Identifies beads and quantifies statistics about bead-to-bead distance,
% bead density, aggregation level, and segmentation area.
% Saves overlay figures into outputFolder (if provided):
%   - beadDetection.fig/.tif/.eps : detected beads (red vector circles) on contrast-enhanced image
%   - beadSegmentation.fig/.tif   : mask regions color-coded by area (um^2) with labeled colorbar
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
% Updated: vector-graphics circles (.fig + .eps), contrast stretch with
%          0.5% saturation, labeled colorbar on segmentation map. 2026.

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

%% Illumination mask ? continuous area where beads are present
% NOTE: the 'mask' returned by pointSourceDetection is a tight per-bead
% segmentation mask (tiny regions around each PSF).  Using that as the
% density denominator would grossly overestimate density because the area
% is far smaller than the true illuminated field.
%
% Instead we build an independent illumination mask:
%   1. Strongly blur the image (large Gaussian, radius ~ 10*sigma) to obtain
%      a smooth background intensity map that is high wherever beads exist
%      in reasonable numbers and low in truly dark regions.
%   2. Threshold at illumThreshFrac * (median of non-zero blurred values)
%      to separate the bright bead-containing region from dark borders.
%   3. Morphologically close with a disk of radius closingRadPx to bridge
%      small gaps between beads and produce one solid illuminated blob.
%   4. Remove very small isolated specks (< minRegionFracOfMax * max region)
%      that are noise rather than real illuminated field.
%
illumSmoothSigma  = 10 * sigma;      % heavy blur to capture field envelope
illumThreshFrac   = 0.15;            % threshold = 15% of median bright value
closingRadPx      = round(5 * sigma);% morphological closing radius (px)
minRegionFracOfMax = 0.1;            % drop regions < 10% of largest region

% Step 1: Gaussian blur
blurred = imgaussfilt(double(img), illumSmoothSigma);

% Step 2: Threshold
blurVals   = blurred(:);
medBright  = median(blurVals(blurVals > 0));
illumMask  = blurred > illumThreshFrac * medBright;

% Step 3: Morphological closing to fill gaps
se        = strel('disk', closingRadPx);
illumMask = imclose(illumMask, se);

% Step 4: Remove small spurious regions
illumCC   = bwconncomp(illumMask);
illumSizes = cellfun(@numel, illumCC.PixelIdxList);
if ~isempty(illumSizes)
    maxSz = max(illumSizes);
    keepIdx = find(illumSizes >= minRegionFracOfMax * maxSz);
    illumMask = false(size(illumMask));
    for ki = keepIdx
        illumMask(illumCC.PixelIdxList{ki}) = true;
    end
end

illumAreaPix = sum(illumMask(:));
if illumAreaPix > 0
    curArea = illumAreaPix * pixSizeUm^2;
else
    % Ultimate fallback: trimmed tile area (border of 4*sigma on each side)
    curAreaPix = (size(img,1) - ceil(4*sigma)*2) * (size(img,2) - ceil(4*sigma)*2);
    curArea    = curAreaPix * pixSizeUm^2;
end

%% Bead density
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

    % --- FIX 2: Contrast-stretch with 0.5% saturation ---
    % imadjust + stretchlim clips the bottom and top 0.5% of pixel values,
    % producing a brighter image with gentle highlight saturation.
    imgDouble = double(img);
    imgNorm   = imgDouble - min(imgDouble(:));
    denom     = max(imgNorm(:));
    if denom == 0, denom = 1; end
    imgNorm   = imgNorm / denom;   % [0,1] double

    satFrac   = 0.005;             % 0.5% saturation at each tail
    loVal     = quantile(imgNorm(:), satFrac);
    hiVal     = quantile(imgNorm(:), 1 - satFrac);
    if hiVal <= loVal, hiVal = loVal + eps; end
    imgStretch = (imgNorm - loVal) / (hiVal - loVal);
    imgStretch = uint8(255 * min(max(imgStretch, 0), 1));   % clamp & convert
    imgRGB     = repmat(imgStretch, [1 1 3]);   % grayscale RGB base

    imgH = size(imgRGB, 1);
    imgW = size(imgRGB, 2);

    % ---- FIX 1) Bead detection ? MATLAB figure with vector circles --------
    % Circles are drawn as line objects (vector graphics) so they remain
    % editable in Illustrator after export. Saved as:
    %   beadDetection.fig  ? editable MATLAB figure
    %   beadDetection.eps  ? vector export for Illustrator
    %   beadDetection.tif  ? raster for quick preview
    hDet = figure('Visible', 'off', 'Color', 'k');
    imshow(imgStretch, [0 255], 'Border', 'tight');
    hold on;
    if ~isempty(pstruct) && ~isempty(pstruct.x)
        r_circ  = 2 * sigma;   % circle radius in pixels
        theta   = linspace(0, 2*pi, 61);
        cosT    = cos(theta);
        sinT    = sin(theta);
        for k = 1:numel(pstruct.x)
            % plot() draws a closed vector circle ? stays scalable in .fig/.eps
            plot(pstruct.x(k) + r_circ * cosT, ...
                 pstruct.y(k) + r_circ * sinT, ...
                 'r-', 'LineWidth', 0.8);
        end
    end
    hold off;
    % Save all three formats
    hgsave(hDet,  fullfile(outputFolder, 'beadDetection'),         '-v7.3');
    print(hDet,   fullfile(outputFolder, 'beadDetection'),         '-depsc2', '-painters');
    print(hDet,   fullfile(outputFolder, 'beadDetection.tif'),     '-dtiff',  '-r150');
    close(hDet);

    % ---- Segmentation map ? black background, parula colormap, labeled colorbar ----
    % Background (non-mask) pixels are pure black; only segmented regions are
    % colored by area using parula.  No grayscale underlay, so the background
    % is unambiguously black rather than brownish-gray.
    % Saved as:
    %   beadSegmentation.fig  ? editable MATLAB figure
    %   beadSegmentation.tif  ? raster for quick preview
    hSeg = figure('Visible', 'off', 'Color', 'k');
    ax = axes('Parent', hSeg, 'Color', 'k');
    if nRegions > 0
        areaMin   = min(regionAreas_um2);
        areaMax   = max(regionAreas_um2);
        areaRange = max(areaMax - areaMin, eps);

        % Build label image: NaN = background (shows as black), [0,1] = parula color
        labelAreaImg = nan(imgH, imgW);
        for r = 1:nRegions
            pixIdx = CC.PixelIdxList{r};
            labelAreaImg(pixIdx) = (regionAreas_um2(r) - areaMin) / areaRange;
        end

        % Render with parula; NaN pixels are fully transparent over the black axes
        hOver = imagesc(ax, labelAreaImg, [0 1]);
        set(hOver, 'AlphaData', double(~isnan(labelAreaImg)));
        colormap(ax, parula(256));
        axis(ax, 'image', 'off');

        % Labeled colorbar with actual area values (um^2), white text on black
        cb = colorbar(ax);
        cb.Color          = 'w';
        cb.Label.String   = 'Segment area (\mum^2)';
        cb.Label.FontSize = 9;
        cb.Label.Color    = 'w';
        midVal = (areaMin + areaMax) / 2;
        cb.Ticks      = [0, 0.5, 1];
        cb.TickLabels = { sprintf('%.3g', areaMin), ...
                          sprintf('%.3g', midVal),  ...
                          sprintf('%.3g', areaMax) };
    else
        % No regions: show a plain black image
        imagesc(ax, zeros(imgH, imgW, 'uint8'));
        colormap(ax, gray(2));
        axis(ax, 'image', 'off');
    end

    % Save
    hgsave(hSeg, fullfile(outputFolder, 'beadSegmentation'),     '-v7.3');
    print(hSeg,  fullfile(outputFolder, 'beadSegmentation.tif'), '-dtiff', '-r150');
    close(hSeg);

    % ---- Illumination mask QC overlay -------------------------------------
    % Semi-transparent green overlay on the contrast-stretched image shows
    % exactly which area is used as the density denominator.  Inspect this
    % to verify the mask parameters (illumThreshFrac, closingRadPx) are
    % appropriate for your data.
    hIllum = figure('Visible', 'off', 'Color', 'k');
    imshow(imgStretch, [0 255], 'Border', 'tight'); hold on;
    hM = imagesc(double(illumMask));
    colormap(gca, [0 0 0; 0 0.8 0]);   % index 0 = black, index 1 = green
    set(hM, 'AlphaData', illumMask * 0.35);   % 35% opacity green overlay
    title(sprintf('Illumination mask  (area = %.1f \\mum^2)', curArea), ...
          'Color', 'w', 'FontSize', 8);
    hold off;
    print(hIllum, fullfile(outputFolder, 'illumMask.tif'), '-dtiff', '-r150');
    close(hIllum);
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