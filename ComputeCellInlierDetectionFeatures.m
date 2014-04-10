function [featureStruct, varargout] = ComputeCellInlierDetectionFeatures( imInput, imLabelCellSeg, cellId, spacing, varargin )

    imdims = ndims( imInput );
    imsize = size( imInput );
    
    p = inputParser;
    p.addRequired( 'imInput', @(x) (isnumeric(x) && ndims(x) == 3) ); 
    p.addRequired( 'imLabelCellSeg', @(x) ( (isnumeric(x) || islogical(x)) && all(size(x) == imsize) ) );
    p.addRequired( 'cellId', @(x) ( isnumeric(x) && isscalar(x) ) );
    p.addRequired( 'spacing', @(x) (isnumeric(x) && numel(x) == 3) ); 
    p.addParamValue( 'imageIntensityRange', [], @(x) (isnumeric(x) && all(size(x) == [3,2])) );
    p.addParamValue( 'glcmIntensityLevels', 32, @(x) (isnumeric(x) && isscalar(x)) );
    p.addParamValue( 'glcmOffsetDistanceVec', [1, 2, 4], @(x) (isnumeric(x)) );
    p.parse(imInput, imLabelCellSeg, cellId, spacing, varargin{:} );
    
    PARAMETERS = p.Results;
    
    % compute some basic cell props
    cellProps = ComputeRegionProperties(imLabelCellSeg, cellId, spacing, 'bboxPadding', [16,16,2] );
    
    % Crop portion of image tightly enclosing the cell with some buffer
    imCellDataCropped = imInput( cellProps.indCropBBox{:} );        
    imLabelCellSegCropped = imLabelCellSeg( cellProps.indCropBBox{:} );
    imCellMaskCropped = (imLabelCellSegCropped == cellId);
    imCellMaskDilatedCropped = imdilate( imCellMaskCropped, ones(5,5,3) );    
    
    % Appearance features
    
        % intensity statistics
        cellPixelIntensities = imCellDataCropped(imCellMaskCropped > 0);                    
             
            % histogram
            imageIntensityRange = ComputeImageIntensityRange( imCellDataCropped(imCellMaskDilatedCropped > 0) ); 
            binLocs = linspace( imageIntensityRange(1), imageIntensityRange(2), 8 );
            curChannelHistogram = histc( cellPixelIntensities, binLocs );
            curChannelHistogram = curChannelHistogram / sum(curChannelHistogram);
            curIntensityStats.entropy = sum(-curChannelHistogram .* log2(curChannelHistogram + eps));
            
            % global intensity
            curIntensityStats.median = median( cellPixelIntensities );
            curIntensityStats.mad = mad( cellPixelIntensities );
            curIntensityStats.iqr = iqr( cellPixelIntensities );
            curIntensityStats.skewness = skewnessRobustHinkley( cellPixelIntensities );
            curIntensityStats.kurtosis = kurtosisRobustCrow( cellPixelIntensities );

            % local intensity stddev        
            imLocalStd = stdfilt( imCellDataCropped );
            curIntensityStats.LocalStddev = median( imLocalStd(imCellMaskCropped > 0) );

            % local range
            imLocalRange = rangefilt( imCellDataCropped );
            curIntensityStats.LocalRange = median( imLocalRange(imCellMaskCropped > 0) );
            
            % contrast between foreground and background
            curIntensityStats.BgndToFgndContrast = median( imLocalStd(~imLabelCellSegCropped) ) / curIntensityStats.median;
            
            appearance.intensity = curIntensityStats;
        
        % texture features
        neighRasterVisitTimes = reshape( 1:3^imdims, 3 * ones(1,imdims) );
        neighOffSetInd = find( neighRasterVisitTimes > (3^imdims + 1)/2 );
        neighOffSetSubInd = ind2submat( size(neighRasterVisitTimes), neighOffSetInd ) - 2;
        OffsetDictionary = neighOffSetSubInd;

        textureCompFunc = @(roiMask, offsets, intensityRange) ( graycopropsext( sum( ...
                                                                    graycomatrixnd( imCellDataCropped, ...
                                                                                    'ROIMask', roiMask, ...  
                                                                                    'Offset', offsets, ...
                                                                                    'NumLevels', PARAMETERS.glcmIntensityLevels, ...
                                                                                    'GrayLimits', intensityRange, ...
                                                                                    'Symmetric', true ), 3 ) ...
                                                                             ) );

        disp = PARAMETERS.glcmOffsetDistanceVec;
        anisotropyFactor = PARAMETERS.spacing(3) / PARAMETERS.spacing(1);

        appearance.texture.haralick_fgnd = [];
        appearance.texture.haralick_bgnd = [];
        
        bboxIntensityRange = [min(imCellDataCropped(:)) max(imCellDataCropped(:))]; 
        
        for i = 1:numel(disp)

            curOffsets = OffsetDictionary * disp(i);
            curOffsets(:,3) = ceil( curOffsets(:,3) / anisotropyFactor );

            % foreground texture
            cellPixelIntensities = imCellDataCropped(imCellMaskDilatedCropped > 0);                    
            cellPixelIntensityRange = [min(cellPixelIntensities) max(cellPixelIntensities)]; 
            curGlcmTextureStruct_fgnd = textureCompFunc(imCellMaskDilatedCropped, curOffsets, cellPixelIntensityRange );
            appearance.texture.haralick_fgnd = [ appearance.texture.haralick_fgnd; curGlcmTextureStruct_fgnd ];

            % background texture
            curGlcmTextureStruct_bgnd = textureCompFunc(~imLabelCellSegCropped, curOffsets, bboxIntensityRange );
            appearance.texture.haralick_bgnd = [ appearance.texture.haralick_bgnd; curGlcmTextureStruct_bgnd ];
            
        end
        
    featureStruct.appearance = appearance;
    
    % geometry
    
        % shape of Z-MIP
        imCurCellMaskMIP = max( imCellMaskCropped, [], 3 );
        imCurCellMaskMIPBnd = bwperim( imCurCellMaskMIP );
        
        mipRegionProps = regionprops( imCurCellMaskMIP, { 'Solidity', 'Eccentricity' } );

        geometry.MaskMIP_ZDir.Convexity = mipRegionProps.Solidity;
        geometry.MaskMIP_ZDir.Eccentricity = mipRegionProps.Eccentricity;
        
        ptCellPixels = ind2submat( size(imCurCellMaskMIP), find(imCurCellMaskMIP > 0) );
        ptCellBoundaryPixels = ind2submat( size(imCurCellMaskMIP), find(imCurCellMaskMIPBnd > 0) );        
        ptCellCentroid = mean(ptCellPixels);
        bndDistances = sqrt( sum( (bsxfun(@minus, ptCellBoundaryPixels, ptCellCentroid)).^2, 2) );
        
        geometry.MaskMIP_ZDir.NonCircularity = std(bndDistances) / mean(bndDistances);
        
        % 3D shape
        if cellProps.ConvexArea > 0
            geometry.Convexity = cellProps.Area / cellProps.ConvexArea;
        else            
            geometry.Convexity =  geometry.MaskMIP_ZDir.Convexity; % single-slice cell
        end
    
        % size
        geometry.volume = cellProps.Area * prod(spacing);
        geometry.bbox_size = cellProps.BoundingBoxTight(:,2) - cellProps.BoundingBoxTight(:,1) + 1;
        geometry.bbox_size = geometry.bbox_size' .* spacing;
        
    featureStruct.geometry = geometry;
    
    if nargout > 1
        varargout{1} = PARAMETERS;
    end
    
end
