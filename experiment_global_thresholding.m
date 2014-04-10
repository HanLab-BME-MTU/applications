function experiment_global_thresholding( imInput, metadata )
%% compare various global thresholding methods

    ImageIntensityRange = ComputeImageDynamicRange( imInput, 99.0 );
    imAdjusted = AdjustImageIntensityRange( imInput, ImageIntensityRange );
    imAdjusted = mat2gray( log( 1 + imAdjusted ) ) * range( ImageIntensityRange );
    imAdjusted = matitk( 'FMEDIAN', [2,2,2], imAdjusted );
    numHistogramBins = round( range(ImageIntensityRange)/5 );
    figure, histogram( imAdjusted );
    
    imThreshKittlerMinerr = callMatitkFilter( 'FKittlerMinimumErrorThreshold', { numHistogramBins } , imAdjusted );    
    imThreshOtsu = callMatitkFilter( 'FOtsuThreshold', { numHistogramBins } , imAdjusted );    
    imThreshIsoData = callMatitkFilter( 'FIsoDataThreshold', { numHistogramBins } , imAdjusted );    
    imThreshPrewittIntermodes = callMatitkFilter( 'FPrewittIntermodesThreshold', { numHistogramBins } , imAdjusted );        
    
    imseriesmaskshow( imInput, {imThreshKittlerMinerr, imThreshOtsu, imThreshIsoData, imThreshPrewittIntermodes } );
    set( gcf, 'Name', 'Clustering: MinError - Otsu - IsoData - Intermodes' );
    
    imThreshKapurEntropy = callMatitkFilter( 'FKapurMaximumEntropyThreshold', { numHistogramBins } , imAdjusted );    
    imThreshRenyiEntropy = callMatitkFilter( 'FRenyiEntropyThreshold', { numHistogramBins } , imAdjusted );    
    imThreshLiEntropy = callMatitkFilter( 'FLiMinimumCrossEntropyThreshold', { numHistogramBins } , imAdjusted );    
    imThreshShanbhagEntropy = callMatitkFilter( 'FShanbhagFuzzyEntropicThreshold', { numHistogramBins } , imAdjusted );    
    imThreshYenEntropy = callMatitkFilter( 'FYenEntropyThreshold', { numHistogramBins } , imAdjusted );        
    
    imseriesmaskshow( imInput, {imThreshKapurEntropy, imThreshRenyiEntropy, imThreshLiEntropy, imThreshShanbhagEntropy, imThreshYenEntropy } );    
    set( gcf, 'Name', 'Entropy: Kapur - Renyi - Li - Shanbhag - Yen' );
    
    imThreshTsaiMoment = callMatitkFilter( 'FTsaiMomentsPreservingThreshold', { numHistogramBins } , imAdjusted );    
    imThreshHuangFuzzySimilarity = callMatitkFilter( 'FHuangFuzzySimilarityThreshold', { numHistogramBins } , imAdjusted );        
    imThreshTriangle = callMatitkFilter( 'FTriangleThreshold', { numHistogramBins } , imAdjusted );        
    imThreshPikazStableState = callMatitkFilter( 'FPikazTopologicalStableStateThreshold', { 500 } , imAdjusted );    
    
    imseriesmaskshow( imInput, {imThreshTsaiMoment, imThreshHuangFuzzySimilarity, imThreshTriangle, imThreshPikazStableState } );    
    set( gcf, 'Name', 'Other: Moment - Huang - Triangle - StableState' );
    
    % best/distinct four  
    imseriesmaskshow( imInput, {imThreshKittlerMinerr, imThreshOtsu, imThreshKapurEntropy, imThreshTsaiMoment} );  
    set( gcf, 'Name', 'Distinct: Minerr - Otsu - KapurEntropy - TsaiMoment' );
      
    % generate montage for 4 different slices
    dispSlices = round( linspace(1,metadata.volSize(3),4) );
    maskAlpha = 0.3;
    maskColor = [1,0,0];
    
    for i = 1:numel( dispSlices )
        
        sliceId = dispSlices(i);
        imSlice = log( imInput(:,:,sliceId) );       
        
        figure;
        set( gcf, 'Name', sprintf( 'Slice - %d - Threshold Results', sliceId ) );
        
        subplot(3,4,1)        
            imRGB = genImageMaskOverlay( imSlice, imThreshKittlerMinerr(:,:,sliceId), maskColor, maskAlpha );
            imshow( imRGB );
            title( 'MinimumError' );
            
        subplot(3,4,2)        
            imRGB = genImageMaskOverlay( imSlice, imThreshOtsu(:,:,sliceId), maskColor, maskAlpha );
            imshow( imRGB );
            title( 'Otsu' );

        subplot(3,4,3)        
            imRGB = genImageMaskOverlay( imSlice, imThreshPrewittIntermodes(:,:,sliceId), maskColor, maskAlpha );
            imshow( imRGB );
            title( 'Intermeans' );

        subplot(3,4,4)        
            imRGB = genImageMaskOverlay( imSlice, imThreshKapurEntropy(:,:,sliceId), maskColor, maskAlpha );
            imshow( imRGB );
            title( 'KapurEntropy' );

        subplot(3,4,5)        
            imRGB = genImageMaskOverlay( imSlice, imThreshRenyiEntropy(:,:,sliceId), maskColor, maskAlpha );
            imshow( imRGB );
            title( 'RenyiEntropy' );

        subplot(3,4,6)        
            imRGB = genImageMaskOverlay( imSlice, imThreshLiEntropy(:,:,sliceId), maskColor, maskAlpha );
            imshow( imRGB );
            title( 'LiEntropy' );

        subplot(3,4,7)        
            imRGB = genImageMaskOverlay( imSlice, imThreshShanbhagEntropy(:,:,sliceId), maskColor, maskAlpha );
            imshow( imRGB );
            title( 'ShanbhagEntropy' );

        subplot(3,4,8)        
            imRGB = genImageMaskOverlay( imSlice, imThreshYenEntropy(:,:,sliceId), maskColor, maskAlpha );
            imshow( imRGB );
            title( 'YenEntropy' );

        subplot(3,4,9)        
            imRGB = genImageMaskOverlay( imSlice, imThreshTsaiMoment(:,:,sliceId), maskColor, maskAlpha );
            imshow( imRGB );
            title( 'TsaiMoment' );

        subplot(3,4,10)        
            imRGB = genImageMaskOverlay( imSlice, imThreshHuangFuzzySimilarity(:,:,sliceId), maskColor, maskAlpha );
            imshow( imRGB );
            title( 'FuzzySimilarity' );

        subplot(3,4,11)        
            imRGB = genImageMaskOverlay( imSlice, imThreshTriangle(:,:,sliceId), maskColor, maskAlpha );
            imshow( imRGB );
            title( 'Triangle' );

        subplot(3,4,12)        
            imRGB = genImageMaskOverlay( imSlice, imThreshPikazStableState(:,:,sliceId), maskColor, maskAlpha );
            imshow( imRGB );
            title( 'StableState' );
            
    end
end