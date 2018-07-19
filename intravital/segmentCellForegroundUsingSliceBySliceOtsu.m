function [imCellForegroundMask] = segmentCellForegroundUsingSliceBySliceOtsu(imInput, varargin)

    p = inputParser;    
    p.addRequired('imInput', @(x) (isnumeric(x) && ismember(ndims(x), [2,3])));
    p.parse(imInput);
    
    p.addParamValue('minLocalGlobalThresholdRatio', 0.6, @(x) isscalar(x) && isnumeric(x));    
    p.addParamValue('debugMode', false, @(x) (isscalar(x) && islogical(x)));    
    p.addParamValue('flagParallelize', false, @(x) (isscalar(x) && islogical(x)));
    p.parse(imInput, varargin{:});
    
    minLocalGlobalThresholdRatio = p.Results.minLocalGlobalThresholdRatio;
    
    flagDebugMode = p.Results.debugMode;
    flagParallelize = p.Results.flagParallelize;

    % open pool of cpu cores if not open already
    if flagParallelize
        flagPoolOpenedAlready = matlabpool('size')> 0;        
        if ~flagPoolOpenedAlready 
            matlabpool open;
        end            
    end
    
    % local otsu thresholding in each slice
    imCellForegroundMask = zeros(size(imInput));
    
    if flagDebugMode
        sliceThreshVals = zeros(size(imInput,3), 1);
    end

    otsuGlobalThresh = thresholdOtsu( imInput );
    
    if flagParallelize

        if flagDebugMode
            totalTimer = tic;
            fprintf('\nPerforming Otsu Thresholding in each of the %d slices ... ', size(imInput, 3));
        end                                            
        
        parfor sliceId = 1:size(imInput, 3)
            
            imSlice = imInput(:,:,sliceId);
            
            curSliceOtsuThresh = thresholdOtsu( imSlice );

            if flagDebugMode
                sliceThreshVals(sliceId) = curSliceOtsuThresh;                
            end
            
            if curSliceOtsuThresh < (minLocalGlobalThresholdRatio * otsuGlobalThresh)
                curSliceOtsuThresh = otsuGlobalThresh;
            end
            
            imCellForegroundMask(:,:,sliceId) = imSlice > curSliceOtsuThresh;
            
        end

        if flagDebugMode
            timeElapsed = toc(totalTimer);
            fprintf('took %f seconds\n', timeElapsed);
        end            
        
    else
        
        if flagDebugMode
            totalTimer = tic;
            fprintf('\nPerforming Otsu Thresholding in each of the %d slices ... \n', size(imInput, 3));
        end
        
        otsuGlobalThresh = thresholdOtsu( imInput );
        
        for sliceId = 1:size(imInput, 3)

            if flagDebugMode
                tic
                fprintf('\n\tThresholding slice %d/%d ... ', sliceId, size(imInput, 3));
            end
            
            imSlice = imInput(:,:,sliceId);
            
            curSliceOtsuThresh = thresholdOtsu( imSlice );

            if flagDebugMode
                sliceThreshVals(sliceId) = curSliceOtsuThresh;                
            end
            
            if curSliceOtsuThresh < (minLocalGlobalThresholdRatio * otsuGlobalThresh)
                curSliceOtsuThresh = otsuGlobalThresh;
            end
            
            imCellForegroundMask(:,:,sliceId) = imSlice > curSliceOtsuThresh;
            
        end
    
        if flagDebugMode
            timeElapsed = toc(totalTimer);
            fprintf('\nTotal time taken - %f seconds\n', timeElapsed);
        end            
        
    end
    
    % post-processing
    diskRad = ones(1,ndims(imInput));
    diskRad(1:2) = 3;
    imCellForegroundMask = imopen(imCellForegroundMask, streldisknd(diskRad));

    % display stuff in debug mode
    if flagDebugMode        
        
        if ndims(imInput) > 2

            figure, plot(1:size(imInput,3), sliceThreshVals, 'b-', 'LineWidth', 2.0 );
            hold on;
                plot([1, size(imInput,3)], otsuGlobalThresh + zeros(1,2), 'r-', 'LineWidth', 2.0 );
                plot([1, size(imInput,3)], minLocalGlobalThresholdRatio * otsuGlobalThresh + zeros(1,2), ...
                     'g-', 'LineWidth', 2.0 );
            hold off;
            
            xlabel( 'Z-slice' );
            ylabel( 'Otsu threshold' );            
            title( 'Variation of otsu threshold with depth' );            
            legend( { 'slice-thresold', 'global-threshold', 'slice-threshold lower-bound' } );
            
        end
        
    end
        
    % close matlab pool if not open before calling this function
    if flagParallelize && ~flagPoolOpenedAlready
        matlabpool close;
    end    
    
end