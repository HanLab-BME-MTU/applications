function experiment_gradient_weighted_watershed( imInput, metadata, flagParallelize )
%% segment cell using gradient weighted watershed 

    if ~exist( 'flagParallelize', 'var' )
        flagParallelize = false;
    end
        
    if flagParallelize
        flagPoolOpenedAlready = matlabpool( 'size' ) > 0;        
        if ~flagPoolOpenedAlready 
            matlabpool open;
        end            
    end

    minObjectDiameter = (6 * ones(1,3)) ./ metadata.voxelSpacing;
    ImageIntensityRange = ComputeImageDynamicRange( imInput, 99.0 );
    imAdjusted = AdjustImageIntensityRange( imInput, ImageIntensityRange );    

    % apply median filter to remove any spotty noise
    imAdjusted = matitk( 'FMEDIAN', [1,1,1], imAdjusted );    
    
    % threshold the image using some algorithm
    imThresh = zeros( size( imAdjusted ) );
    numHistogramBins = 256; 
    thresholdingAlgorithm = 'FPikazTopologicalStableStateThreshold';
    if flagParallelize
        parfor sliceId = 1:size( imInput, 3 )       
            imSlice = ComputeImageLogTransform( imAdjusted(:,:,sliceId) ); 
            imSlice = mat2gray( imSlice );
            imThresh(:,:,sliceId) = callMatitkFilter( thresholdingAlgorithm, { numHistogramBins } , imSlice );                        
        end 
    else        
        for sliceId = 1:size( imInput, 3 )       
            imSlice = ComputeImageLogTransform( imAdjusted(:,:,sliceId) ); 
            imSlice = mat2gray( imSlice );
            imThresh(:,:,sliceId) = callMatitkFilter( thresholdingAlgorithm, { numHistogramBins } , imSlice );                        
        end        
    end
    
    imThresh = imfill( imThresh );
    %imAdjusted( ~imThresh ) = min( imAdjusted(:) );    
    imseriesmaskshow( imInput, {imThresh} );
    set( gcf, 'Name', 'Thresholding Result' );
    
    % compute distance map
    imDistMap = callMatitkFilter( 'FSignedMaurerDistancemap',{} , imThresh, [], [], metadata.voxelSpacing );
    imDistMap( imDistMap > 0 ) = 0;
    imDistMap = -1 * imDistMap;
    
    imseriesmaskshow( imDistMap, imThresh );
    set( gcf, 'Name', 'Distance Map' );        

    % compute gradient magnitude image
    imGMag = callMatitkFilter( 'FGradientMagnitudeSmoothing', { 6 / 3.5 }, imAdjusted, [], [], metadata.voxelSpacing );
    
    imseriesshow( imGMag );
    set( gcf, 'Name', 'Gradient Magnitude' );
    
    % compute gradient weighted distance map
    minGrad = min( imGMag(:) );
    maxGrad = max( imGMag(:) );
    imDistGrad = imDistMap .* exp( 1 - ((imGMag - minGrad) / (maxGrad - minGrad)) );     
    imDistGrad( ~imThresh ) = max( imDistGrad(:) ) + 1;
    imDistGrad = max(imDistGrad(:)) - imDistGrad;
    imDistGrad = matitk( 'FMEDIAN', [2,2,1], imDistGrad );
    
    imseriesshow( imDistGrad );
    set( gcf, 'Name', 'Gradient-weighted Distance Transform' );
    
    % detect seed points
    [imLocalMax] = detect_cell_seeds_IntensityMaxima( imAdjusted, minObjectDiameter );
    imCellSeeds = imLocalMax;
    imCellSeeds( ~imThresh ) = 0;
    
    imseriesmaskshow( imInput, imCellSeeds );
    set( gcf, 'Name', 'Cell Seedpoints' );    
    
    % Overlay/impose the maxima on the inverted distance map image
    imMaximaImposed = imimposemin( imDistGrad, imCellSeeds );
    
    imMaximaImposedDisplay = imMaximaImposed;
    imMaximaImposedDisplay( imMaximaImposed == -Inf ) = min( imMaximaImposed( imMaximaImposed ~= -Inf ) ) - 3 * eps;
    imseriesmaskshow( imMaximaImposedDisplay, {imThresh} );    
    set( gcf, 'Name', 'Local Maxima imposed as minima on the Distance Map' );    

    % Apply Watershed
    L = watershed( imMaximaImposed );
    L( ~imThresh ) = 0;
    imCellMask = double( L > 0 );
    
    if flagParallelize
        parfor sliceId = 1:size( imInput, 3 )       
            imCellMask(:,:,sliceId) = imopen( imCellMask(:,:,sliceId), strel( 'disk', round(0.25 * min(minObjectDiameter(1:2))) ) );
        end                
    else
        for sliceId = 1:size( imInput, 3 )       
            imCellMask(:,:,sliceId) = imopen( imCellMask(:,:,sliceId), strel( 'disk', round(0.25 * min(minObjectDiameter(1:2))) ) );
        end            
    end
    
    L = bwlabeln( imCellMask );
    regstats = regionprops(L, 'Area');
    objAreaList = [regstats.Area];
    smallObjInd = find( objAreaList <= 0.15 * prod(minObjectDiameter) );
    L( ismember( L, smallObjInd ) ) = 0;
    
    segMaskRGB = label2rgbND(L);
    imseriesmaskshowrgb( imInput, segMaskRGB, 'maskAlphas', 0.2 );
    set( gcf, 'Name', 'Result of Marker-controlled Watershed Segmentation' );
    
    if flagParallelize && ~flagPoolOpenedAlready
        matlabpool close;
    end    

end

function old_code()
    % compute distance transform of the threshold result
    imDist = double( bwdist( 1 - imThresh ) );    
   
    % apply waterhshed algorithm on distance transform
    imDistWatershed = imDist;    
    imDistWatershed( ~imThresh ) = max( imDistWatershed(:) ) + 1;
    imDistWatershed = max( imDistWatershed(:) ) - imDistWatershed;    
    imDistWatershed = matitk( 'FMEDIAN', [2,2,1], imDistWatershed );
    
    imseriesmaskshow( imDistWatershed, {imThresh} );
    set( gcf, 'Name', 'Distance Transform' );
    
    imLabelDistanceMap = watershed( imDistWatershed );
    
    imseriesmaskshow( imInput, {imLabelDistanceMap > 1} );
    set( gcf, 'Name', 'Watershed Distance Map' );
    
    % compute gradient magnitude image
    imGMag = callMatitkFilter( 'FGradientMagnitude', {}, imAdjusted, [], [], voxelSpacingCorrected );
    
    imseriesshow( imGMag );
    set( gcf, 'Name', 'Gradient Magnitude' );
    
    % compute gradient weighted distance transform
    minGrad = min( imGMag(:) );
    maxGrad = max( imGMag(:) );
    imDistGrad = imDist .* exp( 1 - ((imGMag - minGrad) / (maxGrad - minGrad)) );     
    imDistGrad( ~imThresh ) = max( imDistGrad(:) ) + 1;
    imDistGrad = max(imDistGrad(:)) - imDistGrad;
    imDistGrad = matitk( 'FMEDIAN', [2,2,1], imDistGrad );
    
    imseriesshow( imDistGrad );
    set( gcf, 'Name', 'Gradient-weighted Distance Transform' );
    
    % apply waterhshed algorithm on gradient weighted distance transform
    imLabelGradWatershed = watershed( imDistGrad );
    
    imseriesmaskshow( imInput, {imLabelGradWatershed > 1} );
    set( gcf, 'Name', 'Watershed Gradient-Weighted Distance Map' );
end