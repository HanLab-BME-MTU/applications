function experiment_gvf_based_seed_detection( imInput, metadata, flagParallelize )

    if ~exist( 'flagParallelize', 'var' )
        flagParallelize = false;
    end
        
    if flagParallelize
        flagPoolOpenedAlready = matlabpool( 'size' ) > 0;        
        if ~flagPoolOpenedAlready 
            matlabpool open;
        end            
    end

    minObjectDiameter = [ 10, 10, 6 ] ./ metadata.voxelSpacing;
    ImageIntensityRange = ComputeImageDynamicRange( imInput, 99.0 );
    imAdjusted = AdjustImageIntensityRange( imInput, ImageIntensityRange );    

    % apply median filter to remove any spotty noise
    imAdjusted = matitk( 'FMEDIAN', [1,1,1], imAdjusted );    
    
    % threshold the image using some algorithm
    imThresh = zeros( size( imAdjusted ) );
    numHistogramBins = 256; 
    thresholdingAlgorithm = 'FOtsuThreshold';
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
    imAdjusted( ~imThresh ) = min( imAdjusted(:) );    
    imseriesmaskshow( imInput, {imThresh} );
    set( gcf, 'Name', 'Thresholding Result' );
    
    % smooth image before looking for locate maxima
    sigma = minObjectDiameter / 3.5;
    imBlurred = matitk_custom( 'FRecursiveGaussianSmoothing', { sigma }, imAdjusted );

    imseriesmaskshow( imBlurred, {imThresh} );
    set( gcf, 'Name', 'Blurred Image' );
    
    % detect cell seed points
    imCellSeedPoints = detect_cell_seed_points_gvf( imAdjusted, minObjectDiameter, 'debugMode', true );
    imCellSeedPoints( ~imThresh ) = 0;
    
        imseriesmaskshow( max(imAdjusted,[],3), max(imCellSeedPoints,[],3) );
        set( gcf, 'Name', 'MIP of Local Intensity Maxima in Blurred Image' );
    
    % Overlay/impose the maxima on the inverted input image
    imMaximaImposed = imimposemin( 1 - mat2gray( imBlurred ), imLocalMax );
    
    imMaximaImposedDisplay = imMaximaImposed;
    imMaximaImposedDisplay( imMaximaImposed == -Inf ) = min( imMaximaImposed( imMaximaImposed ~= -Inf ) ) - 3 * eps;
    imseriesmaskshow( imMaximaImposedDisplay, {imThresh} );    set( gcf, 'Name', 'Local Maxima imposed as minima in Input image' );    
    
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
    
%     % generate segmentation movie
%     ImageIntensityRange = ComputeImageDynamicRange( imInput, 99.0 );
%     imAdjustedLog = ComputeImageLogTransform( AdjustImageIntensityRange( imInput, ImageIntensityRange ) );        
%     displayrange = [min(imAdjustedLog(:)), max(imAdjustedLog(:))];
%     alpha = 0.2;
%     for i = 1:size(imInput,3)
%        imSlice = im2uint8( mat2gray(imAdjustedLog(:,:,i), displayrange) );
%        rgbSlice = cat(3, imSlice, imSlice, imSlice );
%        rgbMaskSlice = squeeze(segMaskRGB(:,:,i,:));
%        rgbSlice = uint8( double((1 - alpha) * rgbSlice) + 255 * alpha * rgbMaskSlice );
%        segMovie(i) = im2frame(rgbSlice); 
%     end
    
    if flagParallelize && ~flagPoolOpenedAlready
        matlabpool close;
    end
    
end