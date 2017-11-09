function mask = segmentNucleiLocalOtsu(imInput, dataProperties, varargin)
% segmentNucleiLocalOtsu uses layer-by-layer local otsu thresholding to
% segment nuclei in 3D
%
% INPUT
%   imInput
% OUTPUT
%   mask:

% 05/2016 Ning Zhang

p = inputParser;
p.addRequired('imInput', @(x) (isnumeric(x) && ~isempty(x)));
p.addRequired('pixelSizeXY', @(x) (isnumeric(x) && numel(x) == 1));
p.addRequired('pixelSizeZ', @(x) (isnumeric(x) && numel(x) == 1));
p.addRequired('imSize', @(x) (isnumeric(x) && numel(x) == 2));
p.addRequired('nDepth', @(x) (isnumeric(x) && numel(x) == 1));
p.addOptional('flagDebugMode', 1, @isnumeric);
p.parse(imInput, dataProperties.PIXELSIZE_XY, dataProperties.PIXELSIZE_Z, ...
    dataProperties.imSize, dataProperties.nDepth, varargin{:})

% Check dataProperties and make sure none of those parameters are zero
imInput = p.Results.imInput;
pixelSizeXY = p.Results.pixelSizeXY;
pixelSizeZ = p.Results.pixelSizeZ;
imSize = p.Results.imSize;
nDepth = p.Results.nDepth;
flagDebugMode = p.Results.flagDebugMode;

% Default parameters
% localThresholdWindowRadius in um (Refer to thresholdLocalCalculate.m by Liya Ding)
% Carefully adjust parameters!!!
localThresholdWindowRadius = 30;
localWindowPaceFraction = 1/3;
minSliceLocalGlobalThresholdRatio = 0.6;
minSliceToStackThresholdRatio = 0.2;

% Convert um to pixel and calculate pace
localWindowRadius = round(localThresholdWindowRadius ./ pixelSizeXY);
localWindowPace = round(localWindowRadius * localWindowPaceFraction);


% local otsu thresholding in each slice
otsuMask = zeros(size(imInput));

% Calculate the threshold for whole 3D stack
globalStackThresh = thresholdOtsu(imInput);

% Debug module
if flagDebugMode
    imLocalThresholdVals = zeros(size(imInput));
    imGlobalSliceThresholdVals = zeros(size(imInput));
    totalTimer = tic;
    fprintf('\nPerforming Local Otsu Thresholding in each of the %d slices ... ', size(imInput, 3));
end

for sliceId = 1:size(imInput, 3)
    imSlice = imInput(:,:,sliceId);
    [sliceGlobalThresh, imMask, imLocalThreshVal] = thresholdLocalSeg(imSlice, 'Otsu', ...
        localWindowRadius, localWindowPace, ...
        minSliceLocalGlobalThresholdRatio * 100);
    
    if sliceGlobalThresh < minSliceToStackThresholdRatio * globalStackThresh
        imMask = imSlice > minSliceToStackThresholdRatio * globalStackThresh;
    end
    otsuMask(:,:,sliceId) = imMask;
    
    if flagDebugMode
        imLocalThresholdVals(:,:,sliceId) = imLocalThreshVal;
        imGlobalSliceThresholdVals(:,:,sliceId) = sliceGlobalThresh;
    end
    
end

if flagDebugMode
    timeElapsed = toc(totalTimer);
    fprintf('took %f seconds\n', timeElapsed);
end

% display stuff in debug mode
if flagDebugMode
    
    if ndims(imInput) > 2
        sliceThreshVals = squeeze(imGlobalSliceThresholdVals(1,1,:));
        figure, plot(1:size(imInput,3), sliceThreshVals, 'b-', 'LineWidth', 2.0 );
        hold on;
        plot(1:numel(sliceThreshVals), globalStackThresh * ones(1,numel(sliceThreshVals)), ...
            'g-', 'LineWidth', 2.0 );
        plot(1:numel(sliceThreshVals), minSliceToStackThresholdRatio * globalStackThresh * ones(1,numel(sliceThreshVals)), ...
            'r-', 'LineWidth', 2.0 );
        
        hold off;
        xlabel( 'Z-slice' );
        ylabel( 'Otsu threshold' );
        title( 'Variation of otsu threshold from slice to slice' );
        legend( { 'slice threshold', 'global threshold', 'slice threshold lower-bnd' } );
    end
    
    imseriesmaskshow(imInput, otsuMask);
    set(gcf, 'Name', sprintf('Otsu mask before post processing'));
    
end


% post-processing from /home2/nzhang/matlab/applications/FISHprobe
% /InvivoCytometer_2.0_source_code/code_package/segmentCellForegroundUsingLocalOtsu.m
diskRad = ones(1,ndims(imInput));
diskRad(1:2) = 3;
otsuMask = imopen(otsuMask, streldisknd(diskRad));

% post-processing from /home2/nzhang/matlab/applications/FISHprobe
% /InvivoCytometer_2.0_source_code/code_package/segmentCellsInIntravitalData.m
otsuMask = imfill( double(otsuMask) );
otsuMask = imclose(otsuMask, streldisknd(2*ones(1,ndims(imInput))));

% Remove unqualified regions
threshL = bwlabeln(otsuMask);
regionStats = regionprops( threshL, {'Area', 'BoundingBox'} );

% Remove regions with unreasonable small/large size
% Find smarter methods for size selection!!!
minVolIndex = 0.5;
maxVolIndex = 0.8;

spacing = [pixelSizeXY, pixelSizeXY, pixelSizeZ];
regArea = [regionStats.Area] * prod(spacing);
medianArea = median(regArea);
minCellVolume = medianArea * (1 - minVolIndex);
maxCellVolume = medianArea * (1 + maxVolIndex);

midianDiameter = (medianArea/pi/4*3)^(1/3)*2;
fprintf('The diameter of median size nucleus is %.2f um.\n', midianDiameter)

dumpRegionIndex = find(regArea < minCellVolume | regArea > maxCellVolume);
threshL(ismember(threshL, dumpRegionIndex)) = 0;

% Remove regions that touch the boundary
flagCrossBoundary = false(1, numel(regionStats));

for i = 1:numel(regionStats)
    if any(regionStats(i).BoundingBox(1:2) <1)
        
        flagCrossBoundary(i) = 1;
        
    else if ((regionStats(i).BoundingBox(1) + regionStats(i).BoundingBox(4)) >= imSize(2)) ...
                || ((regionStats(i).BoundingBox(2) + regionStats(i).BoundingBox(5)) >= imSize(1))
            
            flagCrossBoundary(i) = 1;
            
        end
    end
end

% for i = 1:numel(regionStats)
%     if any(regionStats(i).BoundingBox(1:(ndims(imInput))) <1)
%
%         flagCrossBoundary(i) = 1;
%
%     else if (any((regionStats(i).BoundingBox(1) + regionStats(i).BoundingBox(4)) >= imSize(1))) ...
%             || (any((regionStats(i).BoundingBox(2) + regionStats(i).BoundingBox(5)) >= imSize(2))) ...
%             || ((regionStats(i).BoundingBox(3) + regionStats(i).BoundingBox(6)) >= nDepth)
%
%             flagCrossBoundary(i) = 1;
%
%         end
%     end
% end

CrossBoundaryIndex = find(flagCrossBoundary);
threshL(ismember(threshL, CrossBoundaryIndex)) = 0;

% Remove regions with unclosed holes in 3D that cannot be filled with
% morphological processing
% xyMaxProj = max(threshL,[],3);
% xyMaxProjFilled = imfill(xyMaxProj);
% xyDiff = xyMaxProjFilled - xyMaxProj;
% % How to keep masks with only small holes at the boundary ???
% threshL(ismember(threshL, xyDiff)) = 0;

% yzMaxProj = max(threshL,[],1);
% yzMaxProjFilled = imfill(yzMaxProj);
% yzDiff = yzMaxProjFilled - yzMaxProj;
% threshL(ismember(threshL, yzDiff)) = 0;
%
% xzMaxProj = max(threshL,[],2);
% xzMaxProjFilled = imfill(xzMaxProj);
% xzDiff = xzMaxProjFilled - xzMaxProj;
% threshL(ismember(threshL, xzDiff)) = 0;

mask = bwlabeln(double(threshL > 0));

maskZeroOne = mask > 0;
if flagDebugMode
    imseriesmaskshow(imInput, maskZeroOne);
    set(gcf, 'Name', sprintf('Otsu mask after post processing'));
end

% cleanDiskRad = max( [round(0.25 * min(cellDiameterRange) ./ spacing); [2,2,1]] );
% mask = imopen(mask, streldisknd(cleanDiskRad) ); % removes thin vessel like structures

end