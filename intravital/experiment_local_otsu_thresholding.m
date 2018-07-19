function experiment_local_otsu_thresholding(imInput, metadata, flagParallelize )

    if ~exist( 'flagParallelize', 'var' )
        flagParallelize = false;
    end

    % open pool of cpu cores if not open already
    if flagParallelize
        flagPoolOpenedAlready = matlabpool('size')> 0;        
        if ~flagPoolOpenedAlready 
            matlabpool open;
        end            
    end
    
    % pre-processing and standardization
    ImageIntensityRange = ComputeImageDynamicRange( imInput, 98.0 );
    %imAdjusted = mat2gray( imInput, ImageIntensityRange ) * 4096;
    imAdjusted = mat2gray( imInput ) * 4096;
    imAdjusted = matitk( 'FMEDIAN', [1,1,1], imAdjusted );

    % local thresholding in each slice
    localWindowRadius = round(30 ./ metadata.voxelSpacing(1));
    localWindowPace = round(localWindowRadius / 3);
    minLocalGlobalThresholdRatio = 0.6;
        
    [imCellForegroundMask_LocalOtsu] = segmentCellForegroundUsingLocalOtsu( imAdjusted, localWindowRadius, ...
                                                                            'localWindowPace', localWindowPace, ...
                                                                            'minLocalGlobalThresholdRatio', minLocalGlobalThresholdRatio, ...
                                                                            'flagParallelize', flagParallelize, ...
                                                                            'debugMode', true );
                                                                        
    [imCellForegroundMask_LocalMinerrGauss] = segmentCellForegroundUsingLocalMinError( imAdjusted, localWindowRadius, ...
                                                                                       'model', 'gaussian', ...  
                                                                                       'localWindowPace', localWindowPace, ...
                                                                                       'minLocalGlobalThresholdRatio', minLocalGlobalThresholdRatio, ...
                                                                                       'flagParallelize', flagParallelize, ...
                                                                                       'debugMode', true );

    [imCellForegroundMask_LocalMinerrPoisson] = segmentCellForegroundUsingLocalMinError( imAdjusted, localWindowRadius, ...
                                                                                         'model', 'poisson', ...  
                                                                                         'localWindowPace', localWindowPace, ...
                                                                                         'minLocalGlobalThresholdRatio', minLocalGlobalThresholdRatio, ...
                                                                                         'flagParallelize', flagParallelize, ...
                                                                                         'debugMode', true );
                                                                              
                                                                        
    % slice-by-slice otsu thresholding for comparison
    %imAdjustedLog = mat2gray( ComputeImageLogTransform(imAdjusted) );
    imAdjustedLog = imAdjusted;
    
    [imCellForegroundMask_SliceOtsu] = segmentCellForegroundUsingSliceBySliceOtsu( imAdjustedLog, ...
                                                                                   'minLocalGlobalThresholdRatio', minLocalGlobalThresholdRatio, ...
                                                                                   'flagParallelize', flagParallelize, ...
                                                                                   'debugMode', true );

    % Display the result
    imseriesmaskshow( imInput, ...
                      {imCellForegroundMask_LocalOtsu, imCellForegroundMask_LocalMinerrGauss, imCellForegroundMask_LocalMinerrPoisson, imCellForegroundMask_SliceOtsu}, ...
                      'maskAlphas', 0.2 );
                  
    set( gcf, 'Name', 'cellForegroundMask - Local Otsu, Local Minerr gauss, Local Minerr poisson, Otsu slice-by-slice' );
   
    % close matlab pool if not open before calling this function
    if flagParallelize && ~flagPoolOpenedAlready
        matlabpool close;
    end    
    
end