function [imCellForegroundMask] = segmentCellForegroundUsingLocalOtsu(imInput, localWindowRadius, varargin)

    p = inputParser;    
    p.addRequired('imInput', @(x) (isnumeric(x) && ismember(ndims(x), [2,3])));
    p.addRequired('localWindowRadius', @(x) isscalar(x) && isnumeric(x));
    p.parse(imInput, localWindowRadius);
    
    p.addParamValue('localWindowPace', round(localWindowRadius / 3), @(x) isscalar(x) && isnumeric(x));
    p.addParamValue('minLocalGlobalThresholdRatio', 0.6, @(x) isscalar(x) && isnumeric(x));    
    p.addParamValue('debugMode', false, @(x) (isscalar(x) && islogical(x)));    
    p.addParamValue('flagParallelize', false, @(x) (isscalar(x) && islogical(x)));
    p.parse(imInput, localWindowRadius, varargin{:});
    
    localWindowPace = round(p.Results.localWindowPace);
    minSliceLocalGlobalThresholdRatio = p.Results.minLocalGlobalThresholdRatio;
    minStackLocalGlobalThresholdRatio = 0.4;
    
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
    
    globalStackThresh = thresholdOtsu( imInput );
    
    if flagDebugMode
        imLocalThresholdVals = zeros(size(imInput));
        imGlobalSliceThresholdVals = zeros(size(imInput));
    end    
    
    if flagParallelize

        if flagDebugMode
            totalTimer = tic;
            fprintf('\nPerforming Local Otsu Thresholding in each of the %d slices ... ', size(imInput, 3));
        end
        
        parfor sliceId = 1:size(imInput, 3)
            
            imSlice = imInput(:,:,sliceId);
            
            [globalThreshVal, ...
             imMask, imLocalThreshVal] = thresholdLocalSeg(imSlice, 'Otsu', ...
                                                            localWindowRadius, localWindowPace, ...
                                                            minSliceLocalGlobalThresholdRatio * 100);
                                            
            if globalThreshVal < minStackLocalGlobalThresholdRatio * globalStackThresh
                imMask = imMask > globalStackThresh;
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
            fprintf('\nPerforming Local Otsu Thresholding in each of the %d slices ... \n', size(imInput, 3));
        end
        
        for sliceId = 1:size(imInput, 3)

            if flagDebugMode
                tic
                fprintf('\n\tThresholding slice %d/%d ... ', sliceId, size(imInput, 3));
            end
            
            imSlice = imInput(:,:,sliceId);
            
            [globalThreshVal, ...
             imMask, imLocalThreshVal] = thresholdLocalSeg(imSlice, 'Otsu', ...
                                                            localWindowRadius, localWindowPace, ...
                                                            minSliceLocalGlobalThresholdRatio * 100);
                                            
            if globalThreshVal < minStackLocalGlobalThresholdRatio * globalStackThresh
                imMask = imMask > globalStackThresh;
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
            sliceThreshVals = squeeze(imGlobalSliceThresholdVals(1,1,:));
            figure, plot(1:size(imInput,3), sliceThreshVals, 'b-', 'LineWidth', 2.0 );
            hold on;
                plot(1:numel(sliceThreshVals), globalStackThresh * ones(1,numel(sliceThreshVals)), ...
                             'g-', 'LineWidth', 2.0 );
                plot(1:numel(sliceThreshVals), minStackLocalGlobalThresholdRatio * globalStackThresh * ones(1,numel(sliceThreshVals)), ...
                             'r-', 'LineWidth', 2.0 );
                         
            hold off;
            xlabel( 'Z-slice' );
            ylabel( 'Otsu threshold' );
            title( 'Variation of otsu threshold from slice to slice' );
            legend( { 'slice threshold', 'global threshold', 'slice threshold lower-bnd' } );
        end
        
        imseriesmaskshow(imGlobalSliceThresholdVals, imCellForegroundMask, 'maskAlphas', 0.2);
        set(gcf, 'Name', sprintf('Slice Threshold Map: WindowRadius - %d, minLocalThreshRatio - %.3f', ...
                                  localWindowRadius, minSliceLocalGlobalThresholdRatio));
        
        imseriesmaskshow(imLocalThresholdVals, imCellForegroundMask, 'maskAlphas', 0.2);
        set(gcf, 'Name', sprintf('Local Threshold Map: WindowRadius - %d, minLocalThreshRatio - %.3f', ...
                                    localWindowRadius, minSliceLocalGlobalThresholdRatio));
    end
        
    % close matlab pool if not open before calling this function
    if flagParallelize && ~flagPoolOpenedAlready
        matlabpool close;
    end    
    
end