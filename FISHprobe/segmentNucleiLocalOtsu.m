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
mask = zeros(size(imInput));

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
    mask(:,:,sliceId) = imMask;
end

% post-processing
diskRad = ones(1,ndims(imInput));
diskRad(1:2) = 3;
mask = imopen(mask, streldisknd(diskRad));
% mask = bwlabeln( mask > 0 );


% Possible post-processing from /home2/nzhang/matlab/applications/FISHprobe
% /InvivoCytometer_2.0_source_code/code_package/segmentCellsInIntravitalData.m
% 
mask = imfill( double(mask) ); 
mask = imclose(mask, streldisknd(2*ones(1,ndims(imInput))) );

% remove regions with small and invalid/unusual sizes
threshL = bwlabeln( mask );
thRegStats = regionprops( threshL, {'Area', 'BoundingBox'} );

flagIsBBoxSizeValid = false(1, numel(thRegStats));
imsize = size(imInput);

for i = 1:numel( thRegStats )

    bboxSideLength = thRegStats(i).BoundingBox((ndims(imInput)+1):end);
    flagIsCurBBoxBigEnough = all( bboxSideLength >= minCellBBox );
    if ~flagIsCurBBoxBigEnough
        continue;
    end

    % check if this is the region containing the grid using a heuristic
    bboxSideLengthXY = bboxSideLength(1:2);
    flagIsGridRegion = any(bboxSideLengthXY >= 0.5 * imsize([2,1])) && min(bboxSideLengthXY) / max(bboxSideLengthXY) < 0.25;
    flagIsBBoxSizeValid(i) = ~flagIsGridRegion;

end

regArea = [thRegStats.Area] * prod(spacing);

smallRegInd = find( ~flagIsBBoxSizeValid | regArea < minCellVolume );
threshL( ismember(threshL, smallRegInd) ) = 0;
mask = double(threshL > 0);

cleanDiskRad = max( [round(0.25 * min(cellDiameterRange) ./ spacing); [2,2,1]] );
mask = imopen(mask, streldisknd(cleanDiskRad) ); % removes thin vessel like structures

imAdjusted( ~mask ) = min( imAdjusted(:) );

    
end