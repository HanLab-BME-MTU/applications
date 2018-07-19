function [featureStruct, varargout] = ComputeCellStateClassificationFeatures( imageData, imValidROIMask, imLabelCellSeg, cellId, spacing, varargin )

    imdims = ndims( imageData{1} );
    imsize = size( imageData{1} );
    
    p = inputParser;
    p.addRequired( 'imageData', @(x) (iscell(x) && numel(x) == 3) ); 
    p.addRequired( 'imValidROIMask', @(x) ( all(size(x) == imsize) ) );
    p.addRequired( 'imLabelCellSeg', @(x) ( (isnumeric(x) || islogical(x)) && all(size(x) == imsize) ) );
    p.addRequired( 'cellId', @(x) ( isnumeric(x) && isscalar(x) ) );
    p.addRequired( 'spacing', @(x) (isnumeric(x) && numel(x) == 3) ); 
    p.addParamValue( 'imageIntensityRange', [], @(x) (isnumeric(x) && all(size(x) == [3,2])) );
    p.addParamValue( 'glcmIntensityLevels', 64, @(x) (isnumeric(x) && isscalar(x)) );
    p.addParamValue( 'glcmOffsetDistanceVec', [1, 2, 4, 8], @(x) (isnumeric(x)) );
    p.parse(imageData, imValidROIMask, imLabelCellSeg, cellId, spacing, varargin{:} );
    
    PARAMETERS = p.Results;
    
    imCellMask = (imLabelCellSeg == cellId);
    
    % compute some basic cell props
    cellProps = ComputeRegionProperties( imCellMask, 1, spacing, 'bboxPadding', [5,5,2] );
    
    % Crop portion of image tightly enclosing the cell with some buffer
    imCellDataCropped = imageData;
    for i = 1:3
        imCellDataCropped{i} = imageData{i}( cellProps.indCropBBox{:} );        
    end
    imCellMaskCropped = imCellMask( cellProps.indCropBBox{:} );
    imCellMaskDilatedCropped = imdilate( imCellMaskCropped, ones(3,3,3) );    
    imValidROIMaskCropped = imValidROIMask( cellProps.indCropBBox{:} );
    
    % Appearance features
    if ~isempty(PARAMETERS.imageIntensityRange)
        imageIntensityRange = PARAMETERS.imageIntensityRange;
    else

        imageIntensityRange = zeros(3,2);
        for i = 1:3
             imageIntensityRange(i,:) = ComputeImageIntensityRange(imageData{i}(imValidROIMask > 0));
        end        

    end
    
        % intensity statistics
        for i = 1:3        

            if i > 1
                cellPixelIntensities = imCellDataCropped{i}(imCellMaskCropped & imValidROIMaskCropped);                    
            else
                cellPixelIntensities = imCellDataCropped{i}(imCellMaskCropped);                    
            end
             
            % histogram
            binLocs = linspace( imageIntensityRange(i, 1), imageIntensityRange(i, 2), 8 );
            curChannelHistogram = histc( cellPixelIntensities, binLocs );
            curChannelHistogram = curChannelHistogram / sum(curChannelHistogram);
            curIntensityStats.histogram_entropy = sum(-curChannelHistogram .* log2(curChannelHistogram + eps));
            
            % global intensity
            curIntensityStats.median = median( cellPixelIntensities );
            curIntensityStats.mad = mad( cellPixelIntensities );
            curIntensityStats.iqr = iqr( cellPixelIntensities );
            curIntensityStats.skewness = skewnessRobustHinkley( cellPixelIntensities );
            curIntensityStats.kurtosis = kurtosisRobustCrow( cellPixelIntensities );

            % local intensity        
            imLocalStd = stdfilt( imCellDataCropped{i} );
            curIntensityStats.LocalStddev = median( imLocalStd(imCellMaskCropped > 0) );

            chStr = ['ch_' num2str(i)];
            
            appearance.intensity.(chStr) = curIntensityStats;
            
        end
        
        % meta intensity features
        infVal = 4096;
        chcombs = combnk(1:numel(imageData),2);
        metaFeatureNameList = { 'median', 'mad', 'iqr' };

        for i = 1:numel(metaFeatureNameList)                   
            for j = 1:size(chcombs,1)  
                
                curComb = chcombs(j,:);        
                curFeatureName = sprintf( '%sRatio_ch%d_by_%d', ...
                                          metaFeatureNameList{i}, ...
                                          curComb(2), curComb(1) );    
                                      
                numerVal = appearance.intensity.(['ch_' num2str(curComb(2))]).(metaFeatureNameList{i});
                denomVal = appearance.intensity.(['ch_' num2str(curComb(1))]).(metaFeatureNameList{i});
                curFeatureVal = numerVal / denomVal;    
                curFeatureVal = FixInvalidFeatureVal( curFeatureVal, infVal );
                
                appearance.metaIntensity.(curFeatureName) = curFeatureVal;                
            end
        end
        
        % inter-channel colocalization features
        chcombs = combnk(1:numel(imageData),2);
        
        for i = 1:size(chcombs,1)  
            
            imChannel_1 = mat2gray( imCellDataCropped{chcombs(i,1)}, imageIntensityRange(chcombs(i,1), :) );
            imChannel_2 = mat2gray( imCellDataCropped{chcombs(i,2)}, imageIntensityRange(chcombs(i,2), :) );
            
            curColoc = corrPearson( imChannel_1(imCellMaskCropped & imValidROIMaskCropped), ...
                                    imChannel_2(imCellMaskCropped & imValidROIMaskCropped) );
            
            appearance.colocalization.( sprintf('pearsonCoeff_ch_%d_%d', chcombs(i,1), chcombs(i,2)) ) = curColoc;
            
        end
        
        % texture features
        neighRasterVisitTimes = reshape( 1:3^imdims, 3 * ones(1,imdims) );
        neighOffSetInd = find( neighRasterVisitTimes > (3^imdims + 1)/2 );
        neighOffSetSubInd = ind2submat( size(neighRasterVisitTimes), neighOffSetInd ) - 2;
        OffsetDictionary = neighOffSetSubInd;

        textureCompFunc = @(roiMask, offsets, intensityRange) ( graycopropsext( sum( ...
                                                                    graycomatrixnd( imCellDataCropped{1}, ...
                                                                                    'ROIMask', roiMask, ...  
                                                                                    'Offset', offsets, ...
                                                                                    'NumLevels', PARAMETERS.glcmIntensityLevels, ...
                                                                                    'GrayLimits', intensityRange, ...
                                                                                    'Symmetric', true ), 3 ) ...
                                                                             ) );

        disp = PARAMETERS.glcmOffsetDistanceVec;
        anisotropyFactor = PARAMETERS.spacing(3) / PARAMETERS.spacing(1);

        appearance.texture.haralick_nmzd = [];
        
        for i = 1:numel(disp)

            curOffsets = OffsetDictionary * disp(i);
            curOffsets(:,3) = ceil( curOffsets(:,3) / anisotropyFactor );

            % normalized to the intensity range of the cell pixels
            cellPixelIntensities = imCellDataCropped{1}(imCellMaskDilatedCropped > 0);                    
            cellPixelIntensityRange = [min(cellPixelIntensities) max(cellPixelIntensities)]; 
            curGlcmTextureStruct_nmzd = textureCompFunc(imCellMaskDilatedCropped, curOffsets, cellPixelIntensityRange );
            appearance.texture.haralick_nmzd = [ appearance.texture.haralick_nmzd; curGlcmTextureStruct_nmzd ];
            
        end
        
    featureStruct.appearance = appearance;
    
    % geometric features
    
        % properties of Z-MIP
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
        geometry.Convexity = cellProps.Area / cellProps.ConvexArea;
    
    featureStruct.geometry = geometry;
    
    if nargout > 1
        varargout{1} = PARAMETERS;
    end
    
end

function [ featureValFixed ] = FixInvalidFeatureVal( featureVal, infVal )

    if isnan( featureVal )
        featureValFixed = 0;
    elseif isinf( featureVal )
        featureValFixed = infVal;
    else
        featureValFixed = featureVal;
    end
    
end