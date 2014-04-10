function [regProps] = ComputeRegionProperties( imLabel, labelId, varargin )
%% [regProps] = ComputeRegionProperties( imLabel, labelId )
%
%  Computes some basic region properties:
% 
%   - PixelIdxList
%   - ptPixel
%   - Bounding Box
%   - indCropBox
%   - BoundaryPixelIdxList
%   - Area/Volume
%   - Perimeter/Surface-Area
%   - ConvexArea
%
    p = inputParser;
    p.addRequired( 'imLabel', @(x) (ismember( ndims(x), [2,3])) ); 
    p.addRequired( 'labelId', @(x) (isscalar(x)) );
    p.addOptional( 'spacing', ones(1, ndims(imLabel)), @(x) (isnumeric(x) && numel(x) == ndims(imLabel)) );
    p.addParamValue( 'bboxPadding', 5 * ones(1, ndims(imLabel)), @(x) (isscalar(x) || numel(x) == ndims(imLabel)) );
    p.parse( imLabel, labelId, varargin{:});
    
    bboxPadding = p.Results.bboxPadding;
    
    % pixel indices
    regProps.PixelIdxList = find(imLabel == labelId); 

    % Area
    regProps.Area = numel( regProps.PixelIdxList );
    
    % region centroid
    ptPixelLocations = cell(1, ndims(imLabel));
    [ptPixelLocations{:}] = ind2sub( size(imLabel), regProps.PixelIdxList );
    ptPixelLocations = cell2mat( ptPixelLocations );    
    
    regProps.ptCentroid = round( mean( ptPixelLocations ) );
    
    % region bounding box 
    regProps.BoundingBoxTight = ([ min(ptPixelLocations) ; max(ptPixelLocations) ])';
    
    % construct indices to crop region bounding box from any image  
    regProps.BoundingBox = regProps.BoundingBoxTight;
    
    indCropBBox = cell(1,ndims(imLabel));
    for j = 1:ndims(imLabel)
        indStart = max(1, regProps.BoundingBox(j,1) - bboxPadding(j));
        indEnd = min( size(imLabel,j), regProps.BoundingBox(j,2) + bboxPadding(j) );
        indCropBBox{j} = indStart:indEnd;
        regProps.BoundingBox(j, :) = [indStart, indEnd];        
    end
    regProps.indCropBBox = indCropBBox;
    
    % region boundary pixels 
    imLabelCropped = imLabel( indCropBBox{:} );
    ptBoundaryPixelLocCropped = ind2submat( size(imLabelCropped), find(bwperim(imLabelCropped == labelId)) );
    ptBoundaryPixelLoc = bsxfun(@plus, ptBoundaryPixelLocCropped, (regProps.BoundingBox(:,1))' - 1);
    regProps.BoundaryPixelIdxList = submat2ind(size(imLabel), ptBoundaryPixelLoc); 
    
    % compute perimeter/surface-area
    regProps.Perimeter = numel( regProps.BoundaryPixelIdxList );
    
    % get list of coordinates of pixels on cell convex hull
    pts = ptBoundaryPixelLoc(:, [2, 1, 3:ndims(imLabel)]);

    try
        [K, v] = convhulln( pts );
        regProps.ConvexArea = v;          
    catch
        regProps.ConvexArea = 0;          
    end
    
end