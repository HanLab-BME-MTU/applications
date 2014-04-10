function experiment_segment_adhesion_fibers( imInput, metadata )

    ImageIntensityRange = ComputeImageDynamicRange( imInput, 99.0 );
    imAdjusted = AdjustImageIntensityRange( imInput, ImageIntensityRange );    

    % apply median filter to remove any spotty noise
    imAdjusted = matitk( 'FMEDIAN', [1,1,1], imAdjusted );    
    figure, histogram( imAdjusted );
    
    % use thresholding to suppress background structures
    imThresh = zeros( size( imAdjusted ) );
    
        % apply global thresholding algorithm to compute a rough range for the
        % local threshold to waver in
        numHistogramBins = round( range(ImageIntensityRange)/5 );    
        imThreshOtsu = matitk_custom( 'FOtsuThreshold', { numHistogramBins } , imAdjusted );    
        OtsuThreshold = min( imAdjusted( imThreshOtsu > 0 ) );
        
        % apply local thesholding algorithm
        windowRadius = [50,50,1];
        allowedThresholdRange = [0.75, 1.5] * OtsuThreshold;
        imThresh = threshLocalMean( imAdjusted, windowRadius, 0.0, 'allowedThresholdRange', allowedThresholdRange );
    
    imseriesmaskshow( imInput, {imThresh, imThreshOtsu}, 'spacing', metadata.voxelSpacing );
    set( gcf, 'Name', 'Thresholding Result: Local (Mean) -- Global (Otsu)' );    
    %imAdjusted( ~imThresh ) = min( imAdjusted(:) );
    
    % apply steerable filter 
    zAnisotropyFactor = metadata.voxelSpacing(3)/metadata.voxelSpacing(1);
    [res, theta, nms, pixelScaleMap] = multiscaleSteerableDetector3D( imAdjusted, 1, 2.^(0.5:0.1:1.0), zAnisotropyFactor );
    
    % threshold the response of the steerable filter
    numHistogramBins = round( range(ImageIntensityRange)/5 );    
    imLineSegMask = matitk_custom( 'FOtsuThreshold', { numHistogramBins } , res );    
    LineThreshold = min( res( imLineSegMask > 0 ) );
    
    windowRadius = [50,50,1];
    allowedThresholdRange = [0.75, 1.5] * LineThreshold;
    imLineSegMask = threshLocalMean( res, windowRadius, 0.0, 'allowedThresholdRange', allowedThresholdRange );
    
%     % map line mask points to feature space
%     flagLineMask = imLineSegMask > 0;
%     [yind, xind, zind] = ind2sub( size( imLineSegMask ), find( flagLineMask ) );
%     ptLineMask = [ xind, yind, zind ];
%     lineOrientationVec = [ theta.x1( flagLineMask ), theta.x2( flagLineMask ), theta.x3( flagLineMask ) ];
%     perpVec = ptLineMask - repmat( sum( ptLineMask .* lineOrientationVec, 2 ), 1, 3 ) .* lineOrientationVec;
%     ptLineFeature = perpVec;
%     
%     figure, plot3( ptLineFeature(:,1), ptLineFeature(:,2), ptLineFeature(:,3) );
%     set( gcf, 'Name', 'Curvy Points in Feature Space' );
%     title( 'Curvy Points in Feature Space' );
%     xlabel( 'v1' );
%     ylabel( 'v3' );
%     zlabel( 'Perpendicular Distance to Origin' );
    
    % display result
    imLineSegRGBMask = zeros( [size(imLineSegMask), 3] );
    imLineSegRGBMask(:,:,:,1) = imLineSegMask;

    pixelScaleMap( ~imLineSegMask ) = 0;
    imScaleRGB = reshape( label2rgb( pixelScaleMap(:), 'jet', 'k' ) / 255.0, [size( imInput ), 3] );
   
    imseriesmaskshowrgb( res, imScaleRGB, 'spacing', metadata.voxelSpacing );
    set( gcf, 'Name', 'Response of the Multiscale Steerable Detector' );

    imseriesmaskshowrgb( nms, imScaleRGB, 'spacing', metadata.voxelSpacing );
    set( gcf, 'Name', 'Result of Non-maxima suppression' );

    imseriesmaskshow( imInput, {imLineSegMask, imThresh}, 'spacing', metadata.voxelSpacing );
    set( gcf, 'Name', 'Fiber Segmentation Result' );
    
end