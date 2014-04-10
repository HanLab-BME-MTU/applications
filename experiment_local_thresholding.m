function experiment_local_thresholding(imInput, metadata)
%% compare various local thresholding methods
    ImageIntensityRange = ComputeImageDynamicRange( imInput, 98.0 );
    imAdjusted = mat2gray( imInput, ImageIntensityRange ) * 4096;
    imAdjusted = matitk( 'FMEDIAN', [1,1,1], imAdjusted );
    
    % apply global thresholding algorithm to compute a rough range for the
    % local threshold to waver in
    numHistogramBins = round( range(ImageIntensityRange)/5 );    
    imThreshOtsu = matitk_custom( 'FOtsuThreshold', { numHistogramBins } , imAdjusted );    
    OtsuThreshold = min( imAdjusted( imThreshOtsu > 0 ) );
    figure, histogram( imAdjusted );
    
    % apply local thesholding algorithms
    windowRadius = [60,60,2];
    allowedThresholdRange = [0.50, 1.5] * OtsuThreshold;
    
    imThreshWellner = threshLocalWellner( imAdjusted, windowRadius, 0.85, 'allowedThresholdRange', allowedThresholdRange );
    imThreshNiblack = threshLocalNiblack( imAdjusted, windowRadius, -0.2, 'allowedThresholdRange', allowedThresholdRange );
    imThreshSauvola = threshLocalSauvola( imAdjusted, windowRadius, 0.34, 0.5, 'allowedThresholdRange', allowedThresholdRange );
    imThreshMean = threshLocalMean( imAdjusted, windowRadius, 0.0, 'allowedThresholdRange', allowedThresholdRange );
                                      
    imseriesmaskshow( imAdjusted, {imThreshOtsu, imThreshWellner, imThreshNiblack, imThreshSauvola, imThreshMean} );  
    set( gcf, 'Name', 'Local Threshold: Global (Otsu) - Wellner - Niblack - Sauvola - Mean' );

end