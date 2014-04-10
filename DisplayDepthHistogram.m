function DisplayDepthHistogram( imInput, numIntensityLevels, alpha )

    if ~exist( 'numIntensityLevels', 'var' )
        numIntensityLevels = 256;
    end

    if ~exist( 'alpha', 'var' )
        alpha = 0.001;
    end
    
    imageIntensityRange = ComputeImageIntensityRange( imInput );
    binLocs = linspace( imageIntensityRange(1), imageIntensityRange(2), numIntensityLevels );
    binWidth = binLocs(2) - binLocs(1);
    
    imDepthHistogram = zeros( numel(binLocs), size(imInput,3) );

    qhigh = zeros( size(imInput,3) );
    imHigh = zeros( size(imInput) );
    
    for i = 1:size(imInput,3)
        
        imSlice = imInput(:,:,i);
        curSliceHistogram = histc( imSlice(:), binLocs );
        curSliceHistogram = curSliceHistogram / (eps + sum(curSliceHistogram));
        
        imDepthHistogram(:,i) = curSliceHistogram(:);
        
        imFgndMask = imSlice > thresholdOtsu(imSlice);
        
        qhigh(i) = quantile( imSlice(:), 1-alpha );
        
        imHighMask = imSlice > qhigh(i);
        imHigh(:,:,i) = imHighMask;
        
    end
    
    imseriesmaskshow( imInput, double(imHigh) );
    
    imseriesshow(imDepthHistogram * 100);
    hold on;
        plot(1 + (qhigh - imageIntensityRange(1))/binWidth, 'r-', 'LineWidth', 2.0);
    hold off;
   
    figure, plot( qhigh );
    
end