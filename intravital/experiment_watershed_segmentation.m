function [L, segMaskRGB] = experiment_watershed_segmentation( imInput, metadata, flagParallelize )

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
    
    % locate local intensity maxima in gaussian blurred image
    MaximaSuppressionSize = round( minObjectDiameter/1.5 );
    evenind = (mod( MaximaSuppressionSize, 2 ) == 0);
    MaximaSuppressionSize( evenind ) = MaximaSuppressionSize( evenind ) + 1;    
    imLocalMax = locmax3d(imBlurred, MaximaSuppressionSize);    
    imLocalMax = double(imLocalMax > 0);
    imLocalMax( ~imThresh ) = 0;
    
    imseriesmaskshow( max(imBlurred,[],3), max(imLocalMax,[],3) );
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
    
    imCellSeedPoints(~L) = 0;
    [ seedY, seedX, seedZ ] = ind2sub( size(imCellSeedPoints), find( imCellSeedPoints > 0 ) );    
    cellSeedPointLocations = [ seedX, seedY, seedZ ];

    Display3DDataAndResultsInImaris( imInput, spacing, ...
                                     'rgbMask', segMaskRGB, ...
                                     'spotLocations', imLocalMax, ...
                                     'spotRadius', 3);
    
    if flagParallelize && ~flagPoolOpenedAlready
        matlabpool close;
    end

end

