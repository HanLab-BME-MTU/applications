function [imCellForegroundMask] = segmentCellForegroundUsingLocalMinError(imInput, localWindowRadius, varargin)

    p = inputParser;    
    p.addRequired('imInput', @(x) (isnumeric(x) && ismember(ndims(x), [2,3])));
    p.addRequired('localWindowRadius', @(x) isscalar(x) && isnumeric(x));
    p.parse(imInput, localWindowRadius);

    p.addParamValue('model', 'poisson', @(x) (ismember(x, {'gaussian', 'poisson'})));
    p.addParamValue('localWindowPace', round(localWindowRadius / 3), @(x) isscalar(x) && isnumeric(x));
    p.addParamValue('minLocalGlobalThresholdRatio', 0.6, @(x) isscalar(x) && isnumeric(x));    
    p.addParamValue('minSliceToStackThresholdRatio', 0.4, @(x) isscalar(x) && isnumeric(x));
    p.addParamValue('numHistogramBins', 256, @(x) (isscalar(x) && (x-floor(x)) == 0));
    p.addParamValue('debugMode', false, @(x) (isscalar(x) && islogical(x)));    
    p.addParamValue('flagParallelize', false, @(x) (isscalar(x) && islogical(x)));
    p.parse(imInput, localWindowRadius, varargin{:});
    
    strPdfModel = p.Results.model;    
    localWindowPace = round(p.Results.localWindowPace);
    minSliceLocalGlobalThresholdRatio = p.Results.minLocalGlobalThresholdRatio;   
    numHistogramBins = p.Results.numHistogramBins;
    minSliceToStackThresholdRatio = p.Results.minSliceToStackThresholdRatio;
    
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
    
    threshFunc = @(x) (thresholdMinimumError(x, 'model', strPdfModel, ...
                                                'numHistogramBins', numHistogramBins));
    globalStackThresh = threshFunc( imInput );
    
    if flagDebugMode
        imLocalThresholdVals = zeros(size(imInput));
        imGlobalSliceThresholdVals = zeros(size(imInput));
    end
    
    % perform thresholding
    if flagParallelize

        if flagDebugMode
            totalTimer = tic;
            fprintf('\nPerforming Local Min-error Thresholding with a %s model in each of the %d slices ... ', strPdfModel, size(imInput, 3));
        end                                            
        
        parfor sliceId = 1:size(imInput, 3)
            
            imSlice = imInput(:,:,sliceId);

            threshFunc = @(x) (thresholdMinimumError(x, 'model', strPdfModel, ...
                                                        'numHistogramBins', numHistogramBins));

            
            [globalThreshVal, ...
             imMask, imLocalThreshVal] = thresholdLocalSeg(imSlice, threshFunc, ...
                                                            localWindowRadius, localWindowPace, ...
                                                            minSliceLocalGlobalThresholdRatio * 100);
                                            
            if globalThreshVal < minSliceToStackThresholdRatio * globalStackThresh
                imMask = imMask > globalStackThresh;
                fprintf('\nWARNING: Slice - %d: SliceThreshold/StackThreshold of %.2f is less than the allowed lower-bound of %.2f ... \n', ...
                         sliceId, globalThreshVal/globalStackThresh, minSliceToStackThresholdRatio);
            end
                                                        
            imCellForegroundMask(:,:,sliceId) = imMask;
            
            if flagDebugMode
                imLocalThresholdVals(:,:,sliceId) = imLocalThreshVal;
                imGlobalSliceThresholdVals(:,:,sliceId) = globalThreshVal;
            end
            
        end

        if flagDebugMode
            timeElapsed = toc(totalTimer);
            fprintf('took %f seconds\n', timeElapsed);
        end            
        
    else
        
        if flagDebugMode
            totalTimer = tic;
            fprintf('\nPerforming Local Min-error Thresholding with a %s model in each of the %d slices ... \n', strPdfModel, size(imInput, 3));
        end
        
        for sliceId = 1:size(imInput, 3)

            if flagDebugMode
                tic
                fprintf('\n\tThresholding slice %d/%d ... ', sliceId, size(imInput, 3));
            end
            
            imSlice = imInput(:,:,sliceId);

            threshFunc = @(x) (thresholdMinimumError(x, 'model', strPdfModel, ...
                                                        'numHistogramBins', numHistogramBins));

            [globalThreshVal, ...
             imMask, imLocalThreshVal] = thresholdLocalSeg(imSlice, threshFunc, ...
                                                            localWindowRadius, localWindowPace, ...
                                                            minSliceLocalGlobalThresholdRatio * 100);
                                            
            if globalThreshVal < minSliceToStackThresholdRatio * globalStackThresh
                imMask = imMask > globalStackThresh;
                fprintf('\nWARNING: Slice - %d: SliceThreshold/StackThreshold of %.2f is less than the allowed lower-bound of %.2f ... \n', ...
                         sliceId, globalThreshVal/globalStackThresh, minSliceToStackThresholdRatio);
            end
                                                        
            imCellForegroundMask(:,:,sliceId) = imMask;
            
            if flagDebugMode
                imLocalThresholdVals(:,:,sliceId) = imLocalThreshVal;
                imGlobalSliceThresholdVals(:,:,sliceId) = globalThreshVal;
            end
            
            if flagDebugMode
                timeElapsed = toc;
                fprintf('took %f seconds\n', timeElapsed);
            end            
            
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
            
            globalOtsuThresh = thresholdOtsu( imInput );
            
            threshFunc = @(x) (thresholdMinimumError(x, 'model', strPdfModel, ...
                                                        'numHistogramBins', numHistogramBins));
            sliceThreshVals = squeeze(imGlobalSliceThresholdVals(1,1,:));
            figure, plot(1:size(imInput,3), sliceThreshVals, 'b-', 'LineWidth', 2.0 );
            hold on;
                plot(1:numel(sliceThreshVals), globalStackThresh * ones(1,numel(sliceThreshVals)), ...
                             'g-', 'LineWidth', 2.0 );
                plot(1:numel(sliceThreshVals), minSliceToStackThresholdRatio * globalStackThresh * ones(1,numel(sliceThreshVals)), ...
                             'r-', 'LineWidth', 2.0 );                         
                plot(1:numel(sliceThreshVals), globalOtsuThresh * ones(1,numel(sliceThreshVals)), ...
                             'm-', 'LineWidth', 2.0 );
            hold off;
            xlabel( 'Z-slice' );
            ylabel( 'threshold value' );
            title( sprintf('Variation of %s minimum error threshold from slice to slice', strPdfModel) );
            legend( { 'slice threshold', 'global threshold', 'slice threshold lower-bnd', 'global otsu threshold' } );
            
        end
        
        imseriesmaskshow(imGlobalSliceThresholdVals, imCellForegroundMask, 'maskAlphas', 0.2);
        set(gcf, 'Name', 'Slice Threshold Map');
        
        imseriesmaskshow(imLocalThresholdVals, imCellForegroundMask, 'maskAlphas', 0.2);
        set(gcf, 'Name', sprintf('Local Threshold Map: WindowRadius - %d, minLocalThreshRatio - %.3f', ...
                                    localWindowRadius, minSliceLocalGlobalThresholdRatio));
    end
        
    % close matlab pool if not open before calling this function
    if flagParallelize && ~flagPoolOpenedAlready
        matlabpool close;
    end    
    
end