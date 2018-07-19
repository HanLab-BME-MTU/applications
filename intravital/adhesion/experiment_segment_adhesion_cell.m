function experiment_segment_adhesion_cell( imInput, metadata, flagParallelize )

    if ~exist( 'flagParallelize', 'var' )
        flagParallelize = false;
    end
        
    if flagParallelize
        flagPoolOpenedAlready = matlabpool( 'size' ) > 0;        
        if ~flagPoolOpenedAlready 
            matlabpool open;
        end            
    end
    
    ImageIntensityRange = ComputeImageDynamicRange( imInput, 99.0 );
    imAdjusted = AdjustImageIntensityRange( imInput, ImageIntensityRange );    

    % apply median filter to remove any spotty noise
    imAdjusted = matitk( 'FMEDIAN', [2,2,2], imAdjusted );    
    
    % threshold the image using some algorithm
    imThresh = zeros( size( imAdjusted ) );
    numHistogramBins = 256;    
    if flagParallelize
        parfor sliceId = 1:size( imInput, 3 )       
            imSlice = ComputeImageLogTransform( imAdjusted(:,:,sliceId) ); 
            imSlice = mat2gray( imSlice );
            imSliceThresh = callMatitkFilter( 'FRenyiEntropyThreshold', { numHistogramBins } , imSlice );    
            imSliceThresh = imopen( imSliceThresh, strel( 'disk', 5 ) );
            imSliceThresh = imclose( imSliceThresh, strel( 'disk', 5 ) );
            imSliceThresh = imdilate( imSliceThresh, strel( 'disk', 3 ) );
            imThresh(:,:,sliceId) = imSliceThresh;
        end        
    else
        for sliceId = 1:size( imInput, 3 )       
            imSlice = ComputeImageLogTransform( imAdjusted(:,:,sliceId) ); 
            imSlice = mat2gray( imSlice );
            imSliceThresh = callMatitkFilter( 'FRenyiEntropyThreshold', { numHistogramBins } , imSlice );    
            imSliceThresh = imopen( imSliceThresh, strel( 'disk', 5 ) );
            imSliceThresh = imclose( imSliceThresh, strel( 'disk', 5 ) );
            imSliceThresh = imdilate( imSliceThresh, strel( 'disk', 3 ) );
            imThresh(:,:,sliceId) = imSliceThresh;
        end                
    end
        
    
    imThresh = imfill( imThresh );    
    
    % pick the largest connected component   
    [L, numComponents] = bwlabeln( imThresh );    
    if numComponents > 1 
        regstats = regionprops( L, 'Area' );
        [ maxarea, maxobjind ] = max( [ regstats.Area ] );
        imThresh = double(L == maxobjind);
    end
    
    imseriesmaskshow( imInput, imThresh );
    set( gcf, 'Name', 'Thresholding Result' );    
    
    % compute distance map
    imDistMap = callMatitkFilter( 'FSignedMaurerDistancemap',{} , imThresh, [], [], metadata.voxelSpacing  );
    imDistMap( imDistMap > 0 ) = 0;
    imDistMap = -1 * imDistMap;
    imseriesmaskshow( imDistMap, imThresh );
    set( gcf, 'Name', 'Distance Map' );    
    
    % compute cell skeleton from distance map
    % display frontiers
    numQuantizationLevels = 5;
    imDistMapQuantized = ceil( mat2gray( imDistMap ) * numQuantizationLevels );    
    imDistMapFrontiersRGBMask = label2rgbND( imDistMapQuantized );
    imseriesmaskshowrgb( imInput, imDistMapFrontiersRGBMask, 'maskAlphas', 0.5 );
    set( gcf, 'Name', 'Distance Map Frontiers' );    

    % detect spots in each slice
    sigma = 4;
    METHOD = 3;    
    
    switch METHOD
        
        case 1
            
            imLoG = zeros( size( imInput ) );
            imLocalMaxima = zeros( size( imInput ) );
            imLOGPrunedMask = zeros( size( imInput ) );

            if flagParallelize
                parfor sliceId = 1:size( imInput, 3 )        
                    imSlice = imInput(:,:,sliceId);
                    [pstruct, imLOGPrunedMask(:,:,sliceId), imLocalMaxima(:,:,sliceId), imLoG(:,:,sliceId)] = pointSourceDetection(imSlice,sigma);            
                end
            else
                for sliceId = 1:size( imInput, 3 )        
                    imSlice = imInput(:,:,sliceId);
                    [pstruct, imLOGPrunedMask(:,:,sliceId), imLocalMaxima(:,:,sliceId), imLoG(:,:,sliceId)] = pointSourceDetection(imSlice,sigma);            
                end                
            end
            
            imseriesmaskshow( imLoG, {imLocalMaxima, imLOGPrunedMask, imThresh} );
            set( gcf, 'Name', sprintf( 'LOG Result, sigma - %f', sigma ) );    

            imseriesmaskshow( imInput, {imLocalMaxima, imLOGPrunedMask, imThresh} );
            set( gcf, 'Name', sprintf( 'Overlay of Detected spots, sigma - %f', sigma ) );        

            imseriesmaskshow( max(imInput,[],3), max(imLocalMaxima,[],3) );
            set( gcf, 'Name', sprintf( 'MIP Overlay of Detected spots, sigma - %f', sigma ) );            
            
        case 2

            MaxDetectionWindow = round( (3*sigma) * [1,1,metadata.voxelSpacing(1)/metadata.voxelSpacing(3)] );
            evenind = mod( MaxDetectionWindow, 2 ) == 0;
            MaxDetectionWindow(evenind) = MaxDetectionWindow(evenind) + 1;
            
            % run LOG filter
            imLOG = callMatitkFilter( 'FLaplacianOfGaussian', {sigma}, imInput );
            imLOG = max( imLOG(:) ) - imLOG;
            imseriesmaskshowrgb( imLOG, {imDistMapFrontiersRGBMask} );    
            set( gcf, 'Name', sprintf( 'LOG Result, sigma - %f', sigma ) );
            
            % detect local maxima in LOG
            imLOGMasked = imLOG;
            imLocalMaxima = locmax3d(imLOGMasked, MaxDetectionWindow);    
            imLocalMaxima( ~imThresh ) = 0;
            imLocalMaxima( imDistMapQuantized > 2 ) = 0;
            imLocalMaxima = double(imLocalMaxima > 0);

            imseriesmaskshow( max(imInput,[],3), max(imLocalMaxima,[],3) );
            set( gcf, 'Name', sprintf( 'MIP Overlay of Detected spots, sigma - %f', sigma ) );
            
            imseriesmaskshow( imInput, {imLocalMaxima, imThresh} );
            set( gcf, 'Name', sprintf( 'Overlay of Detected spots, sigma - %f', sigma ) );
            
        case 3
            
            MaxDetectionWindow = round( (3*sigma) * [1,1,metadata.voxelSpacing(1)/metadata.voxelSpacing(3)] );
            evenind = mod( MaxDetectionWindow, 2 ) == 0;
            MaxDetectionWindow(evenind) = MaxDetectionWindow(evenind) + 1;
            
            % run LOG filter
            sigmaVec = sigma-1:sigma+1;
            imLOG = ones(size(imInput));
            for i = 1:numel(sigmaVec)
                imLOG = imLOG .* callMatitkFilter( 'FLaplacianOfGaussian', {sigma}, imInput );
            end
            imLOG = max( imLOG(:) ) - imLOG;
            imseriesmaskshowrgb( imLOG, {imDistMapFrontiersRGBMask} );    
            set( gcf, 'Name', sprintf( 'LOG Result, sigma - %f', sigma ) );
            
            % detect local maxima in LOG
            imLOGMasked = imLOG;
            imLocalMaxima = locmax3d(imLOGMasked, MaxDetectionWindow);    
            imLocalMaxima( ~imThresh ) = 0;
            imLocalMaxima( imDistMapQuantized > 2 ) = 0;
            imLocalMaxima = double(imLocalMaxima > 0);

            imseriesmaskshow( max(imInput,[],3), max(imLocalMaxima,[],3) );
            set( gcf, 'Name', sprintf( 'MIP Overlay of Detected spots, sigma - %f', sigma ) );
            
            imseriesmaskshow( imInput, {imLocalMaxima, imThresh} );
            set( gcf, 'Name', sprintf( 'Overlay of Detected spots, sigma - %f', sigma ) );
            
    end   
   
    if flagParallelize && ~flagPoolOpenedAlready
        matlabpool close;
    end
    
end