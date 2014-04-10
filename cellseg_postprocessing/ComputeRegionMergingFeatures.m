function [featureStruct, varargout] = ComputeRegionMergingFeatures(imInput, imLabelCellSeg, cellStats, c1, c2, varargin)
    
    p = inputParser;
    p.addRequired( 'imInput', @(x) (ismember( ndims(x), [2,3])) ); 
    p.addRequired( 'imLabelCellSeg', @(x) ( all(size(x) == size(imInput)) ) );
    p.addRequired( 'cellStats', @(x) ( isstruct(x) && all(isfield(x, {'PixelIdxList', 'BoundingBox'})) ));
    p.addRequired( 'c1', @(x) ( isnumeric(x) && isscalar(x) && ~(x - floor(x) > 0) && x <= max(imLabelCellSeg(:)) ));
    p.addRequired( 'c2', @(x) ( isnumeric(x) && isscalar(x) && ~(x - floor(x) > 0) && x <= max(imLabelCellSeg(:)) ));

    p.addParamValue( 'spacing', ones( 1, ndims(imInput) ), @(x) (isnumeric(x) && numel(x) == ndims(imInput)) );
    p.addParamValue( 'numIntensityLevels', 64, @(x) (isnumeric(x) && isscalar(x)) );
    p.addParamValue( 'glcmOffsetDistanceVec', [1, 2, 4], @(x) (isnumeric(x)) );
    p.addParamValue( 'junctionCurvatureStelRadPhysp', 2, @(x) (isnumeric(x)) );
    
    p.parse(imInput, imLabelCellSeg, cellStats, c1, c2, varargin{:});

    PARAMETERS = p.Results;
    
    imdims = ndims(imInput);
    
    % Crop portion of image enclosing the two cells for faster processing
    ptBBoxCorners = [ cellStats(c1).BoundingBox, cellStats(c2).BoundingBox ];
    indCropBox = cell(1, ndims(imInput));
    for i = 1:imdims
        indCropBox{i} = min(ptBBoxCorners(i,:)):max(ptBBoxCorners(i,:));        
    end
    imInputCropped = imInput( indCropBox{:} );
    imLabelCellSegCropped = imLabelCellSeg( indCropBox{:} );
    
    
    %% Appearance similarity
    imageIntensityRange = ComputeImageIntensityRange( imInput );        
    
%     imCellUnionMaskCropped = ismember(imLabelCellSegCropped, [c1, c2]);
%     imCellUnionMaskDilatedCropped = imdilate(imCellUnionMaskCropped, ones(5 * ones(1,ndims(imInput))));
%     imageIntensityRange = ComputeImageIntensityRange( imInputCropped(imCellUnionMaskDilatedCropped > 0) );
    
        % difference in overall brightness
        median_c1 = median( imInput( cellStats(c1).PixelIdxList ) );
        median_c2 = median( imInput( cellStats(c2).PixelIdxList ) );
        
        appearance.brightnessSimilarity = min([median_c1, median_c2]) / max([median_c1, median_c2]);
        
        % similarity between intensity histograms
        binLocs = linspace( imageIntensityRange(1), imageIntensityRange(2), PARAMETERS.numIntensityLevels );
        hist_c1 = hist( imInput( cellStats(c1).PixelIdxList ), binLocs );
        hist_c2 = hist( imInput( cellStats(c2).PixelIdxList ), binLocs );

        hist_c1 = hist_c1 / sum(hist_c1);
        hist_c2 = hist_c2 / sum(hist_c2);

        appearance.intensityDistSimilarity = dot( sqrt( hist_c1 ), sqrt( hist_c2 ) ); % bhattacharya-coefficient

        % texture similarity
        neighRasterVisitTimes = reshape( 1:3^imdims, 3 * ones(1,imdims) );
        neighOffSetInd = find( neighRasterVisitTimes > (3^imdims + 1)/2 );
        neighOffSetSubInd = ind2submat( size(neighRasterVisitTimes), neighOffSetInd ) - 2;
        OffsetDictionary = neighOffSetSubInd;

        textureCompFunc = @(roiMask, offsets) ( graycopropsext( sum( ...
                                                    graycomatrixnd( imInputCropped, ...
                                                                    'ROIMask', roiMask, ...  
                                                                    'Offset', offsets, ...
                                                                    'NumLevels', PARAMETERS.numIntensityLevels, ...
                                                                    'GrayLimits', imageIntensityRange, ...
                                                                    'Symmetric', true ), 3 ) ...
                                                             ) );

                                                 
        disp = PARAMETERS.glcmOffsetDistanceVec;

        glcmTextureVec_c1 = [];
        glcmTextureVec_c2 = [];

        for i = 1:numel(disp)

            curOffsets = OffsetDictionary * disp(i);

            if imdims == 3
                anisotropyFactor = PARAMETERS.spacing(3) / PARAMETERS.spacing(1);
                curOffsets(:,3) = ceil( curOffsets(:,3) / anisotropyFactor );
            end
            
            curGlcmTextureStruct_c1 = textureCompFunc( imLabelCellSegCropped == c1, curOffsets );
            curGlcmTextureVec_c1 = cell2mat( ConvertFeatureStructToFeatureVec( curGlcmTextureStruct_c1 ) );                
            glcmTextureVec_c1 = cat(2, glcmTextureVec_c1, curGlcmTextureVec_c1);

            curGlcmTextureStruct_c2 = textureCompFunc( imLabelCellSegCropped == c2, curOffsets );
            curGlcmTextureVec_c2 = cell2mat( ConvertFeatureStructToFeatureVec( curGlcmTextureStruct_c2 ) );                
            glcmTextureVec_c2 = cat(2, glcmTextureVec_c2, curGlcmTextureVec_c2);
            
        end
        
        glcmTextureVec_c1 = glcmTextureVec_c1 / norm(glcmTextureVec_c1);
        glcmTextureVec_c2 = glcmTextureVec_c2 / norm(glcmTextureVec_c2);
        
        appearance.textureSimilarity = dot( glcmTextureVec_c1 , glcmTextureVec_c2 ); % cosine similarity measure
    
    featureStruct.appearance = appearance;
    
    %% Boundary Saliency    
%     neighStrel = ones(3 * ones(1,ndims(imInputCropped)));
%     imMask_c1 = imclose( imfill( imLabelCellSegCropped == c1, 'holes' ), neighStrel );
%     imMask_c2 = imclose( imfill( imLabelCellSegCropped == c2, 'holes' ), neighStrel );
% 
%     imComponentCellMaskUnion = imMask_c1 | imMask_c2;
%     imMergedCellMask = imclose(imComponentCellMaskUnion, neighStrel);
% 
%     L = bwlabeln(imMergedCellMask & ~imComponentCellMaskUnion); 
%     stats = regionprops( L, 'Area' );
%     [~, maxind] = max( [stats.Area] );
%     imWatershed = (L == maxind);
%     imWatershedDilated = imdilate(imWatershed, neighStrel) & imMergedCellMask;

    neighStrel = ones(5 * ones(1,ndims(imInput)));

    imMask_c1 = imclose( imfill( imLabelCellSegCropped == c1, 'holes' ), neighStrel );
    imMask_c2 = imclose( imfill( imLabelCellSegCropped == c2, 'holes' ), neighStrel );
    imComponentCellMaskUnion = imMask_c1 | imMask_c2;
    
    imMergedCellMask = imclose(imComponentCellMaskUnion, neighStrel);
    imWatershed = imMergedCellMask & ~imComponentCellMaskUnion & imdilate(imMask_c1, neighStrel) & imdilate(imMask_c2, neighStrel);
    
    if ~any( imWatershed(:) ) % degenerate cases at the boundary
        imMask_c1 = imfill( imLabelCellSegCropped == c1, 'holes' );
        imMask_c2 = imfill( imLabelCellSegCropped == c2, 'holes' );
        imComponentCellMaskUnion = imMask_c1 | imMask_c2;
        imWatershed = imMergedCellMask & ~imComponentCellMaskUnion & imdilate(imMask_c1, neighStrel) & imdilate(imMask_c2, neighStrel);
    end
    
    dilStrel = ones(3 * ones(1,ndims(imInput)));
    imWatershedDilated = imdilate(imWatershed, dilStrel) & imMergedCellMask;

    median_c1 = median( imInput( cellStats(c1).PixelIdxList ) );
    median_c2 = median( imInput( cellStats(c2).PixelIdxList ) );
    median_wshed = median( imInputCropped(imWatershedDilated > 0) );
    
    watershedSaliency(1) = min( min(median_c1, median_c2) / (eps + median_wshed), 4096 );
    watershedSaliency(2) = min( max(median_c1, median_c2) / (eps + median_wshed), 4096 );
    
    featureStruct.boundarySaliency.watershedSaliency = watershedSaliency;
    
    iqr_wshed = iqr( imInputCropped(imWatershedDilated > 0) );
    iqr_c1 = iqr( imInput( cellStats(c1).PixelIdxList ) );
    iqr_c2 = iqr( imInput( cellStats(c2).PixelIdxList ) );
    featureStruct.boundarySaliency.watershedIQR_Ratio(1) = min( iqr_wshed / (eps + min(iqr_c1, iqr_c2)), 4096 );
    featureStruct.boundarySaliency.watershedIQR_Ratio(2) = min( iqr_wshed / (eps + max(iqr_c1, iqr_c2)), 4096 );
    
    %% Boundary Continuity
        
        % amount of boundary overlap
        imPerim_c1 = bwperim( imMask_c1 );        
        imPerim_c2 = bwperim( imMask_c2 );
        
        boundaryOverlap_c1 = numel( find( imWatershedDilated & imPerim_c1 ) );
        boundaryOverlap_c1 = boundaryOverlap_c1 / numel( find( imPerim_c1 ) );
        
        boundaryOverlap_c2 = numel( find( imWatershedDilated & imPerim_c2 ) );
        boundaryOverlap_c2 = boundaryOverlap_c2 / numel( find( imPerim_c2 ) );
        
        boundaryContinuity.overlapAmount = sort( [boundaryOverlap_c1, boundaryOverlap_c2] );
        
        % junction curvature
        bndContSrelRadImsp = max( [ceil(PARAMETERS.junctionCurvatureStelRadPhysp ./ PARAMETERS.spacing); 3 * ones(1, ndims(imInput))] );
        bndContStrel = streldisknd( bndContSrelRadImsp ); 
        bndContStrel = bndContStrel / sum(bndContStrel(:));

        imLocalCurvature = imfilter( double(imMergedCellMask), bndContStrel);
        imLocalCurvature(imLocalCurvature > 0.99) = 0;    
        imJunction = bwperim(imMergedCellMask) & imWatershed;
        junctionCurvatureVals = abs(imLocalCurvature(imJunction > 0) - 0.5);
        
        junctionCurvature.mean = mean(junctionCurvatureVals);
        junctionCurvature.stddev = std(junctionCurvatureVals);
        junctionCurvature.min = min(junctionCurvatureVals); 
        junctionCurvature.max = max(junctionCurvatureVals);

        boundaryContinuity.junctionCurvature = junctionCurvature;

        % gap/distance between the boundaries of the two surfaces
        if ndims( imInput ) > 2
            
           imBndInterface_c1 = imWatershedDilated & imPerim_c1; 
           imBndInterface_c2 = imWatershedDilated & imPerim_c2; 
           
           imBndInterfaceDistMap_c1 = bwdistsc(imBndInterface_c1, PARAMETERS.spacing);              
           distValsToBndInterface_c2 = imBndInterfaceDistMap_c1(imBndInterface_c2 > 0);
           
           boundaryContinuity.interBoundaryDist.mean = mean( distValsToBndInterface_c2 );
           boundaryContinuity.interBoundaryDist.stddev = std( distValsToBndInterface_c2 );
           boundaryContinuity.interBoundaryDist.min = min( distValsToBndInterface_c2 );
           boundaryContinuity.interBoundaryDist.max = max( distValsToBndInterface_c2 );
           
        end
        
    featureStruct.boundaryContinuity =  boundaryContinuity;
    
    % shape
    ptMergedRegionBoundary = ind2submat( size(imMergedCellMask), find( bwperim(imMergedCellMask) > 0 ) );
    [~, mergedRegionConvexArea] = convhulln( ptMergedRegionBoundary(:, [2, 1, 3:ndims(imInput)]) );
    mergedRegionArea = numel( find( imMergedCellMask > 0 ) );
    mergedRegionConvexity =  mergedRegionArea / (eps + mergedRegionConvexArea);
    
    convexity_c1 = cellStats(c1).Area / (eps + cellStats(c1).ConvexArea);
    convexity_c2 = cellStats(c2).Area / (eps + cellStats(c2).ConvexArea);
    
    convexityChange = sort( [ mergedRegionConvexity / convexity_c1, mergedRegionConvexity / convexity_c2] );
    
    featureStruct.shape.mergedRegionConvexity = mergedRegionConvexity;
    featureStruct.shape.convexityChange = convexityChange;
    
    % size
%     featureStruct.size.volumeBeforeMerging = prod( PARAMETERS.spacing ) * sort( [numel(cellStats(c1).PixelIdxList), numel(cellStats(c2).PixelIdxList)] );
%     featureStruct.size.volumeAfterMerging = mergedRegionArea * prod( PARAMETERS.spacing );
%     
%     radii_c1 = ComputeEllipsoidRadii(imLabelCellSegCropped == c1, PARAMETERS.spacing);
%     radii_c2 = ComputeEllipsoidRadii(imLabelCellSegCropped == c2, PARAMETERS.spacing);
%     radii_merged = ComputeEllipsoidRadii(imMergedCellMask, PARAMETERS.spacing);
%     
%     featureStruct.size.majorAxisLengthBeforeMerging = sort( [max(radii_c1), max(radii_c2)] );
%     featureStruct.size.majorAxisLengthAfterMerging = max(radii_merged);
    
    if nargout > 1
       varargout{1} = PARAMETERS; 
    end
end