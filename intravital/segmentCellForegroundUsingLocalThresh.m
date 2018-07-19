function [imCellForegroundMask] = segmentCellForegroundUsingLocalThresh(imInput, localWindowRadius, varargin)

    p = inputParser;    
    p.addRequired('imInput', @(x) (isnumeric(x) && ismember(ndims(x), [2,3])));
    p.addRequired('localWindowRadius', @(x) isscalar(x) && isnumeric(x));
    p.parse(imInput, localWindowRadius);
    
    p.addParamValue('localWindowPace', round(localWindowRadius / 3), @(x) isscalar(x) && isnumeric(x));
    p.addParamValue('minLocalGlobalThresholdRatio', 0.5, @(x) isscalar(x) && isnumeric(x));    
    p.addParamValue('debugMode', false, @(x) (isscalar(x) && islogical(x)));    
    p.addParamValue('flagParallelize', false, @(x) (isscalar(x) && islogical(x)));
    p.parse(imInput, localWindowRadius, varargin{:});
    
    localWindowPace = round(p.Results.localWindowPace);
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
        imLocalThresholdVals = zeros(size(imInput));
        imGlobalSliceThresholdVals = zeros(size(imInput));
    end
    
    strPdfModel = 'poisson';
    threshFunc = @(x) (thresholdMinimumError(x, 'model', strPdfModel));
    
    if flagParallelize

        if flagDebugMode
            tic
            fprintf('\nPerforming Local Thresholding in each of %d slices ... ', size(imInput, 3));
        end
        
        parfor sliceId = 1:size(imInput, 3)
            
            imSlice = imInput(:,:,sliceId);
            
            otsuGlobalThresh = thresholdOtsu( imSlice );
            
            [globalThreshVal, ...
             imMask, imLocalThreshVal] = thresholdLocalSeg(imSlice, threshFunc, ...
                                                            localWindowRadius, localWindowPace, ...
                                                            minLocalGlobalThresholdRatio * 100);
                                            
            imMask = imMask .* (imSlice >= minLocalGlobalThresholdRatio .* otsuGlobalThresh);            
            imCellForegroundMask(:,:,sliceId) = imMask;
            
            if flagDebugMode
                imLocalThresholdVals(:,:,sliceId) = imLocalThreshVal;
                imGlobalSliceThresholdVals(:,:,sliceId) = globalThreshVal;
            end
            
        end

        if flagDebugMode
            timeElapsed = toc;
            fprintf('took %f seconds\n', timeElapsed);
        end            
    
    else
        
        for sliceId = 1:size(imInput, 3)

            if flagDebugMode
                tic
                fprintf('\nThresholding slice %d/%d ... ', sliceId, size(imInput, 3));
            end
            
            imSlice = imInput(:,:,sliceId);
            
            otsuGlobalThresh = thresholdOtsu( imSlice );
            
            [globalThreshVal, ...
             imMask, imLocalThreshVal] = thresholdLocalSeg(imSlice, threshFunc, ...
                                                            localWindowRadius, localWindowPace, ...
                                                            minLocalGlobalThresholdRatio * 100);
                                            
            imMask = imMask .* (imSlice >= minLocalGlobalThresholdRatio .* otsuGlobalThresh);            
            imCellForegroundMask(:,:,sliceId) = imMask;
            
            if flagDebugMode
                imLocalThresholdVals(:,:,sliceId) = imLocalThreshVal;
                imGlobalSliceThresholdVals(:,:,sliceId) = globalThreshVal;
                timeElapsed = toc;
                fprintf('took %f seconds\n', timeElapsed);
            end            
            
        end
    
    end
    
    % post-processing
    diskRad = ones(1,ndims(imInput));
    diskRad(1:2) = 3;
    imCellForegroundMask = imopen(imCellForegroundMask, streldisknd(diskRad));

    % display stuff in debug mode
    if flagDebugMode        
        
        if ndims(imInput) > 2
            globalStackThresh = thresholdOtsu( imInput );
            sliceThreshVals = squeeze(imGlobalSliceThresholdVals(1,1,:));
            figure, plot(1:size(imInput,3), sliceThreshVals, 'b-', 'LineWidth', 2.0 );
            hold on;
                plot(1:numel(sliceThreshVals), globalStackThresh * ones(1,numel(sliceThreshVals)), ...
                             'g-', 'LineWidth', 2.0 );
            hold off;
            xlabel( 'Z-slice' );
            ylabel( 'Otsu threshold' );
            title( 'Variation of otsu threshold from slice to slice' );
            legend( { 'slice threshold', 'global threshold' } );
        end
        
        imseriesmaskshow(imGlobalSliceThresholdVals, imCellForegroundMask, 'maskAlphas', 0.2);
        set(gcf, 'Name', 'Slice Threshold Map');

        imseriesmaskshow(imLocalThresholdVals, imCellForegroundMask, 'maskAlphas', 0.2);
        set(gcf, 'Name', sprintf('Local Threshold Map: WindowRadius - %d, minLocalThreshRatio - %.3f', ...
                                    localWindowRadius, minLocalGlobalThresholdRatio));
    end
        
    % close matlab pool if not open before calling this function
    if flagParallelize && ~flagPoolOpenedAlready
        matlabpool close;
    end    
    
end