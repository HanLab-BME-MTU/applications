function mask = segmentNucleiLocalOtsu(imInput, dataProperties)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% 
% 03/2016 Ning

% Default parameters
% localThresholdWindowRadius in um (see /matlab/applications/FISHprobe/
% InvivoCytometer_2.0_source_code/code_package/thresholdLocalCalculate.m)
parameters.localThresholdWindowRadius = 30;
parameters.localWindowPaceFraction = 1/3;
parameters.minSliceLocalGlobalThresholdRatio = 0.6;
parameters.minSliceToStackThresholdRatio = 0.4;

% Convert um to pixel and calculate pace
localWindowRadius = round(parameters.localThresholdWindowRadius ./ dataProperties.PIXELSIZE_XY);
localWindowPace = round(localWindowRadius * parameters.localWindowPaceFraction);


% local otsu thresholding in each slice
otsuMask = zeros(size(imInput));

% Calculate the threshold for whole 3D stack
globalStackThresh = thresholdOtsu( imInput );

for sliceId = 1:size(imInput, 3)
    imSlice = imInput(:,:,sliceId);
    [sliceGlobalThresh, imMask, imLocalThreshVal] = thresholdLocalSeg(imSlice, 'Otsu', ...
                                                    localWindowRadius, localWindowPace, ...
                                                    parameters.minSliceLocalGlobalThresholdRatio * 100);

    if sliceGlobalThresh < parameters.minSliceToStackThresholdRatio * globalStackThresh
        imMask = imMask > globalStackThresh;
    end
    otsuMask(:,:,sliceId) = imMask;
end

% post-processing from /home2/nzhang/matlab/applications/FISHprobe
% /InvivoCytometer_2.0_source_code/code_package/segmentCellForegroundUsingLocalOtsu.m
diskRad = ones(1,ndims(imInput));
diskRad(1:2) = 3;
otsuMask = imopen(otsuMask, streldisknd(diskRad));

% post-processing from /home2/nzhang/matlab/applications/FISHprobe
% /InvivoCytometer_2.0_source_code/code_package/segmentCellsInIntravitalData.m
otsuMask = imfill( double(otsuMask) ); 
otsuMask = imclose(otsuMask, streldisknd(2*ones(1,ndims(imInput))) );

% Remove unqualified regions
threshL = bwlabeln(otsuMask);
regionStats = regionprops( threshL, {'Area', 'BoundingBox'} );

% Remove regions with unreasonable small/large size
minCellDiameter = input('Enter the minimum cell diameter in um for selection > ');
maxCellDiameter = input('Enter the maximum cell diameter in um for selection > ');
minCellVolume = (4/3) * pi * (0.5 * minCellDiameter)^3;
maxCellVolume = (4/3) * pi * (0.5 * maxCellDiameter)^3;
spacing = [dataProperties.PIXELSIZE_XY, dataProperties.PIXELSIZE_XY, dataProperties.PIXELSIZE_Z];

regArea = [regionStats.Area] * prod(spacing);
smallRegionIndex = find(regArea < minCellVolume | regArea > maxCellVolume);
threshL(ismember(threshL, smallRegionIndex)) = 0;

% Remove regions that touch the boundary
flagCrossBoundary = false(1, numel(regionStats));
if dataProperties.imSize(1) ~= dataProperties.imSize(2)
    error('SizeX and SizeY are not identical.')
end
for i = 1:numel(regionStats)
    if any(regionStats(i).BoundingBox(1:(ndims(imInput))) <1)
        flagCrossBoundary(i) = 1;
    else if any((regionStats(i).BoundingBox(1:2) + regionStats(i).BoundingBox(4:5))...
                > dataProperties.imSize(1)) || ((regionStats(i).BoundingBox(3) + ...
                regionStats(i).BoundingBox(6)) > dataProperties.nDepth)
            flagCrossBoundary(i) = 1;
        end
    end
end
CrossBoundaryIndex = find(flagCrossBoundary);
threshL(ismember(threshL, CrossBoundaryIndex)) = 0;

% Remove regions with unclosed holes in 3D that cannot be filled with
% morphological processing
xyMaxProj = max(threshL,[],3);
xyMaxProjFilled = imfill(xyMaxProj);
xyDiff = xyMaxProjFilled - xyMaxProj;
threshL(ismember(threshL, xyDiff)) = 0;

yzMaxProj = max(threshL,[],1);
yzMaxProjFilled = imfill(yzMaxProj);
yzDiff = yzMaxProjFilled - yzMaxProj;
threshL(ismember(threshL, yzDiff)) = 0;

xzMaxProj = max(threshL,[],2);
xzMaxProjFilled = imfill(xzMaxProj);
xzDiff = xzMaxProjFilled - xzMaxProj;
threshL(ismember(threshL, xzDiff)) = 0;



mask = bwlabeln(double(threshL > 0));

% cleanDiskRad = max( [round(0.25 * min(cellDiameterRange) ./ spacing); [2,2,1]] );
% mask = imopen(mask, streldisknd(cleanDiskRad) ); % removes thin vessel like structures

end