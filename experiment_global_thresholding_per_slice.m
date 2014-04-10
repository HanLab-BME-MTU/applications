function experiment_global_thresholding_per_slice( imInput, metadata )
%% apply global thresholding algorithms in each slice separately

    ImageIntensityRange = ComputeImageDynamicRange( imInput, 99.0 );
    imAdjusted = AdjustImageIntensityRange( imInput, ImageIntensityRange );
    imAdjusted = mat2gray( log( 1 + imAdjusted ) ) * range(ImageIntensityRange);
    imAdjusted = matitk( 'FMEDIAN', [2,2,2], imAdjusted );
    numHistogramBins = round( range(ImageIntensityRange)/5 );
    
    imThreshKittlerMinerr = zeros( size( imAdjusted ) );
    imThreshOtsu = zeros( size( imAdjusted ) );
    imThreshIsoData = zeros( size( imAdjusted ) );
    imThreshPrewittIntermodes = zeros( size( imAdjusted ) );

    imThreshKapurEntropy = zeros( size( imAdjusted ) );
    imThreshRenyiEntropy = zeros( size( imAdjusted ) );
    imThreshLiEntropy = zeros( size( imAdjusted ) );
    imThreshShanbhagEntropy = zeros( size( imAdjusted ) );
    imThreshYenEntropy = zeros( size( imAdjusted ) );

    imThreshTsaiMoment = zeros( size( imAdjusted ) );
    imThreshHuangFuzzySimilarity = zeros( size( imAdjusted ) );
    imThreshTriangle = zeros( size( imAdjusted ) );
    imThreshPikazStableState = zeros( size( imAdjusted ) );
    
    hHistPlot = figure;
    
    for sliceId = 1:size( imInput, 3 )
       
        imSlice = imAdjusted(:,:,sliceId); 

        [p,x] = hist( imSlice, numHistogramBins );
        p = p / sum(p);
        figure(hHistPlot);
        hold all;
            csp = csaps(x,p);
            fnplt(csp);
        
        imThreshKittlerMinerr(:,:,sliceId) = callMatitkFilter( 'FKittlerMinimumErrorThreshold', { numHistogramBins } , imSlice );    
        imThreshOtsu(:,:,sliceId) = callMatitkFilter( 'FOtsuThreshold', { numHistogramBins } , imSlice );    
        imThreshIsoData(:,:,sliceId) = callMatitkFilter( 'FIsoDataThreshold', { numHistogramBins } , imSlice );    
        imThreshPrewittIntermodes(:,:,sliceId) = callMatitkFilter( 'FPrewittIntermodesThreshold', { numHistogramBins } , imSlice );        

        imThreshKapurEntropy(:,:,sliceId) = callMatitkFilter( 'FKapurMaximumEntropyThreshold', { numHistogramBins } , imSlice );    
        imThreshRenyiEntropy(:,:,sliceId) = callMatitkFilter( 'FRenyiEntropyThreshold', { numHistogramBins } , imSlice );    
        imThreshLiEntropy(:,:,sliceId) = callMatitkFilter( 'FLiMinimumCrossEntropyThreshold', { numHistogramBins } , imSlice );    
        imThreshShanbhagEntropy(:,:,sliceId) = callMatitkFilter( 'FShanbhagFuzzyEntropicThreshold', { numHistogramBins } , imSlice );    
        imThreshYenEntropy(:,:,sliceId) = callMatitkFilter( 'FYenEntropyThreshold', { numHistogramBins } , imSlice );        

        imThreshTsaiMoment(:,:,sliceId) = callMatitkFilter( 'FTsaiMomentsPreservingThreshold', { numHistogramBins } , imSlice );    
        imThreshHuangFuzzySimilarity(:,:,sliceId) = callMatitkFilter( 'FHuangFuzzySimilarityThreshold', { numHistogramBins } , imSlice );        
        imThreshTriangle(:,:,sliceId) = callMatitkFilter( 'FTriangleThreshold', { numHistogramBins } , imSlice );        
        imThreshPikazStableState(:,:,sliceId) = callMatitkFilter( 'FPikazTopologicalStableStateThreshold', { 500 } , imSlice );    
        
    end
    
    imseriesmaskshow( imInput, {imThreshKittlerMinerr, imThreshOtsu, imThreshIsoData, imThreshPrewittIntermodes } );
    set( gcf, 'Name', 'Clustering: MinError - Otsu - IsoData - Intermodes' );
    
    imseriesmaskshow( imInput, {imThreshKapurEntropy, imThreshRenyiEntropy, imThreshLiEntropy, imThreshShanbhagEntropy, imThreshYenEntropy } );    
    set( gcf, 'Name', 'Entropy: Kapur - Renyi - Li - Shanbhag - Yen' );

    imseriesmaskshow( imInput, {imThreshTsaiMoment, imThreshHuangFuzzySimilarity, imThreshTriangle, imThreshPikazStableState } );    
    set( gcf, 'Name', 'Other: Moment - Huang - Triangle - StableState' );

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