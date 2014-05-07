function [ imLabelCellSeg, varargout ] = segmentCellsInIntravitalData( imInput, spacing, varargin )

    p = inputParser;
    p.CaseSensitive = false;
    p.addRequired( 'imInput', @(x) ( ismember( ndims(x), [2,3] ) ) );
    p.parse( imInput );
    
    imdims = ndims(imInput);   
    
    p.addRequired( 'spacing', @(x) (numel(x) == imdims) );   
    
    p.addParamValue( 'thresholdingAlgorithm', ...
                     'MinErrorPoissonSliceBySliceLocal', ...
                     @(x) (ischar(x) && ismember(x, {'OtsuGlobalSliceBySliceHybrid', 'OtsuSliceBySliceLocal', 'MinErrorPoissonSliceBySliceLocal', 'BackgroudRemovalUsingMorphologicalOpening'}) ));
    p.addParamValue( 'localThresholdWindowRadiusPhysp', 30, @(x) (isnumeric(x) && isscalar(x)) );
    p.addParamValue( 'minLocalGlobalThresholdRatio', 0.6, @(x) (isnumeric(x) && isscalar(x)) );
    p.addParamValue( 'minSignalToBackgroundRatio', 2.0, @(x) (isnumeric(x) && isscalar(x)) );
    
    p.addParamValue( 'seedPointDetectionAlgorithm', ...
                     'AdaptiveMultiscaleLoG', ...
                     @(x) (ischar(x) && ismember(x,{'IntensityMaxima', 'MultiscaleLoG', 'MultiscaleLoBG', 'AdaptiveMultiscaleLoG'})) );
    p.addParamValue( 'seedDetectorResponseCutoff', [], @(x) (isnumeric(x) && isscalar(x)) );
                 
    p.addParamValue( 'cellDiameterRange', [12, 20], @(x) (isnumeric(x) && numel(x) == 2) );
    p.addParamValue( 'minCellVolume', [], @(x) (isnumeric(x) && isscalar(x)) );
    p.addParamValue( 'minCellBBoxSizePhysp', [3.5 3.5 8], @(x) (isnumeric(x) && numel(x) == 3));
    
    p.addParamValue( 'regionMergingModelFile', [], @(x) (isempty(x) || (ischar(x) && exist(x, 'file'))) );
    
    p.addParamValue( 'roiMask', [], @(x) ( (isnumeric(x) || islogical(x)) && ndims(x) == ndims(imInput) && all(size(x) == size(imInput)) ) )    
    p.addParamValue( 'minCellROIOverlap', 0.5, @(x) (isscalar(x) && x >= 0.0 && x <= 1.0) ) 
    p.addParamValue( 'flagIgnoreCellsOnXYBorder', true, @(x) (isscalar(x)) );
    
    p.addParamValue( 'flagParallelize', false, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue( 'flagDebugMode', false, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue( 'flagDisplayResultsInImaris', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse( imInput, spacing, varargin{:} );    

    % capture parameters
    PARAMETERS = p.Results;
    
    cellDiameterRange = PARAMETERS.cellDiameterRange;
    minCellVolume = PARAMETERS.minCellVolume;
    if isempty( PARAMETERS.minCellVolume )
        minCellVolume = 0.75 * (4/3) * pi * (0.5 * PARAMETERS.cellDiameterRange(1))^3;    
    end
    
    thresholdingAlgorithm = PARAMETERS.thresholdingAlgorithm;
    localThresholdWindowRadiusPhysp = PARAMETERS.localThresholdWindowRadiusPhysp;
    minLocalGlobalThresholdRatio = PARAMETERS.minLocalGlobalThresholdRatio;

    seedPointDetectionAlgorithm = PARAMETERS.seedPointDetectionAlgorithm;
    seedDetectorResponseCutoff = PARAMETERS.seedDetectorResponseCutoff;
    
    roiMask = PARAMETERS.roiMask;    
    minCellROIOverlap = PARAMETERS.minCellROIOverlap;
    flagIgnoreCellsOnXYBorder = PARAMETERS.flagIgnoreCellsOnXYBorder;
    
    regionMergingModelFile = PARAMETERS.regionMergingModelFile;
    
    PARAMETERS.minCellBBoxSizeImsp = max( [ PARAMETERS.minCellBBoxSizePhysp ./ spacing; 3 * ones(1,3) ] );
    
    flagDebugMode = PARAMETERS.flagDebugMode;
    flagDisplayResultsInImaris = PARAMETERS.flagDisplayResultsInImaris;
    flagParallelize = PARAMETERS.flagParallelize;
    if flagParallelize
        flagPoolOpenedAlready = matlabpool( 'size' ) > 0;        
        if ~flagPoolOpenedAlready 
            matlabpool open;
        end            
    end
    
    totalSegTimer = tic;

    % standardize the image
    imInputStandardized = mat2gray( imInput ) * 4096;    
    
    % apply median filter to remove any spotty noise
    fprintf( '\n\n>> Appling median filter to remove noise ...\n\n' );
    
    imAdjusted = matitk( 'FMEDIAN', [1, 1, 1], imInputStandardized );
    
    %% threshold the image using some algorithm
    fprintf( '\n\n>> Running a thresholding algorithm ...\n\n' );    
    
    localWindowRadius = round(localThresholdWindowRadiusPhysp ./ spacing(1));
    localWindowPace = round(localWindowRadius / 3);
    
    threshTimer = tic;
    
    switch thresholdingAlgorithm
        
        case 'OtsuGlobalSliceBySliceHybrid'

            imThresh = segmentCellForegroundUsingSliceBySliceOtsu( imAdjusted, ...
                                                                   'minLocalGlobalThresholdRatio', minLocalGlobalThresholdRatio, ...
                                                                   'flagParallelize', flagParallelize );          
                                                               
        case 'OtsuSliceBySliceLocal'
            
            imThresh = segmentCellForegroundUsingLocalOtsu( imAdjusted, localWindowRadius, ...
                                                            'localWindowPace', localWindowPace, ...
                                                            'minLocalGlobalThresholdRatio', minLocalGlobalThresholdRatio, ...
                                                            'flagParallelize', flagParallelize );
                                                               

        case 'MinErrorPoissonSliceBySliceLocal'
            
                                                        
            imThresh = segmentCellForegroundUsingLocalMinError( imAdjusted, localWindowRadius, ...
                                                                'model', 'poisson', ...  
                                                                'localWindowPace', localWindowPace, ...
                                                                'minLocalGlobalThresholdRatio', minLocalGlobalThresholdRatio, ...
                                                                'flagParallelize', flagParallelize );

                                                            
        case 'BackgroudRemovalUsingMorphologicalOpening'
            
            krnlMax = streldisknd( round(0.5 * max(cellDiameterRange) ./ spacing(1:2)) ); 
            imLocalBackground = imopen(imAdjusted, krnlMax);
            imSignalToBackgroundRatio = imAdjusted ./ (eps + imLocalBackground);
            imThresh = double(imSignalToBackgroundRatio >= PARAMETERS.minSignalToBackgroundRatio);
            
        otherwise
            
            error( 'ERROR: invalid thresholding method' );
        
    end

    imThresh = imfill( imThresh ); 
    imThresh = imclose(imThresh, streldisknd(2*ones(1,ndims(imInput))) );
    
    % remove regions with small and invalid/unusual sizes
    threshL = bwlabeln( imThresh );
    thRegStats = regionprops( threshL, {'Area', 'BoundingBox'} );
    
    flagIsBBoxSizeValid = false(1, numel(thRegStats));
    imsize = size(imInput);
    minBBox = max( [0.5 * min(cellDiameterRange) ./ spacing; PARAMETERS.minCellBBoxSizeImsp] );
    
    for i = 1:numel( thRegStats )
        
        bboxSideLength = thRegStats(i).BoundingBox((ndims(imInput)+1):end);
        flagIsCurBBoxBigEnough = all( bboxSideLength >= minBBox );
        if ~flagIsCurBBoxBigEnough
            continue;
        end
        
        % check if this is the region containing the grid using a heuristic
        bboxSideLengthXY = bboxSideLength(1:2);
        flagIsGridRegion = any(bboxSideLengthXY >= 0.5 * imsize([2,1])) && min(bboxSideLengthXY) / max(bboxSideLengthXY) < 0.25;
        flagIsBBoxSizeValid(i) = ~flagIsGridRegion;
        
    end
    
    regArea = [thRegStats.Area] * prod(spacing);
    
    smallRegInd = find( ~flagIsBBoxSizeValid | regArea < minCellVolume );
    threshL( ismember(threshL, smallRegInd) ) = 0;
    imThresh = double(threshL > 0);

    cleanDiskRad = max( [round(0.25 * min(cellDiameterRange) ./ spacing); [2,2,1]] );
    imThresh = imopen(imThresh, streldisknd(cleanDiskRad) ); % removes thin vessel like structures
    
    imAdjusted( ~imThresh ) = min( imAdjusted(:) );
    
    fprintf( '\nThresholding took a total of %f seconds\n', toc(threshTimer) );
    
    if flagDebugMode
        imseriesmaskshow( imInput, {imThresh}, 'spacing', spacing );
        set( gcf, 'Name', 'Thresholding Result' );
    end
    
    %% detect cell seed points using multiscale LoG   
    fprintf( '\n\n>> Detecting seed points ...\n\n' );

    seedTimer = tic;
    
    maxNumLoGScales = 10;
    
    if isempty(seedDetectorResponseCutoff)
        seedDetectionResponseCutoff = (65.0/4096) * max( imInputStandardized(:) );
    end
    
    seedToBackgroundDistanceCutoff = 0.5 * 0.5 * min(cellDiameterRange);
    
    imDistMap = bwdistsc( 1 - imThresh, spacing );
    
    switch seedPointDetectionAlgorithm

        case 'IntensityMaxima'
            
            [imCellSeedPoints, imResponse] = detect_cell_seeds_IntensityMaxima( imAdjusted, ...
                                                                               mean( cellDiameterRange ), ...
                                                                               'spacing', spacing, ...
                                                                               'debugMode', flagDebugMode );
                                                                           
        case 'MultiscaleLoG'            
            
            numLoGScales = min(maxNumLoGScales, cellDiameterRange(2)-cellDiameterRange(1)+1);
            [imCellSeedPoints, imResponse] = detectBlobsUsingMultiscaleLoG( imAdjusted, ...
                                                                            cellDiameterRange, ...
                                                                            'numLoGScales', numLoGScales, ...
                                                                            'spacing', spacing, ...
                                                                            'debugMode', flagDebugMode );

        case 'AdaptiveMultiscaleLoG'
            
            
            numLoGScales = min(maxNumLoGScales, cellDiameterRange(2)-cellDiameterRange(1)+1);
            [imCellSeedPoints, imResponse] = detectBlobsUsingAdaptiveMultiscaleLoG( imAdjusted, ...
                                                                                    imDistMap, ... 
                                                                                    'blobDiameterRange', cellDiameterRange, ...
                                                                                    'numLoGScales', numLoGScales, ...
                                                                                    'spacing', spacing, ...
                                                                                    'debugMode', flagDebugMode, ...
                                                                                    'showPlots', false);

        case 'MultiscaleLoBG'
            
            numLoGScales = min(maxNumLoGScales, cellDiameterRange(2)-cellDiameterRange(1)+1);
            [imCellSeedPoints, imResponse] = detectBlobsUsingMultiscaleLoBG( imAdjusted, ...
                                                                             cellDiameterRange, ...
                                                                             'numLoGScales', numLoGScales, ...
                                                                             'spacing', spacing, ...
                                                                             'debugMode', flagDebugMode );

    end

    imCellSeedPoints(~imThresh) = 0;    
    imCellSeedPoints = double( imCellSeedPoints > 0 );
    
    % remove weak seed points
    imCellSeedPoints( imCellSeedPoints > 0 & imResponse < seedDetectionResponseCutoff ) = 0; % makes things faster
    
    % remove seed points very close to the region borders
    imCellSeedPoints( imCellSeedPoints > 0 & imDistMap < seedToBackgroundDistanceCutoff) = 0;

    % remove regions without a seed point or with weak seed responses
    threshL = bwlabeln( imThresh );
    fgndRegProps = regionprops(threshL, 'PixelIdxList' );
    
    flagPruneRegion = false( numel(fgndRegProps), 1 );
    for i = 1:numel(fgndRegProps)
        
        curRgnSeedPixInd = fgndRegProps(i).PixelIdxList( imCellSeedPoints( fgndRegProps(i).PixelIdxList ) > 0 );
        
        % zero-out regions without seed-point
        if isempty(curRgnSeedPixInd)
            flagPruneRegion(i) = true;
            continue;
        end
        
        % zero-out regions with all seed points with response below a speficied cutoff
        if max( imResponse(curRgnSeedPixInd) ) < seedDetectionResponseCutoff
            flagPruneRegion(i) = true;             
            continue;
        end

        % zero-out regions with max distance map value below a specified cutoff
        if max(imDistMap(fgndRegProps(i).PixelIdxList)) < seedToBackgroundDistanceCutoff
            flagPruneRegion(i) = true;             
            continue;
        end
        
    end
    threshL( ismember(threshL, find(flagPruneRegion)) ) = 0;
    imThresh = double( threshL > 0 );
    
    fprintf( '\nSeed Detection took a total of %f seconds\n', toc(seedTimer) );

    seedMask = imdilate(imCellSeedPoints > 0, ones(3,3,3));
    if flagDebugMode
        imseriesmaskshow( imInput, {seedMask, imThresh} );
        set( gcf, 'Name', 'Seed Points Overlayed on the Input Image' );
        
        imseriesmaskshow( imResponse, {seedMask, imThresh} );
        set( gcf, 'Name', 'Seed Points Overlayed on Seed Detection Filter Response' );
    end
    
    %% Apply marker-based Watershed
    fprintf( '\n\n>> Running watershed algorithm ...\n\n' );

    imFeature = 1 - mat2gray( imResponse );
    imFeature(~imThresh) = 100; % some high value
    imMinimaImposedAtSeeds = imimposemin( imFeature, imCellSeedPoints );
    L = watershed( imMinimaImposedAtSeeds );
    L( ~imThresh ) = 0; % cookie-cut with threshold mask
    
%     if flagDebugMode
%         imMinimaImposedAtSeedsDisplay = imMinimaImposedAtSeeds;
%         imMinimaImposedAtSeedsDisplay( imMinimaImposedAtSeeds == -Inf ) = min( imMinimaImposedAtSeeds( imMinimaImposedAtSeeds ~= -Inf ) ) - 3 * eps;
%         imseriesmaskshow( imMinimaImposedAtSeedsDisplay, {imThresh}, 'spacing', spacing );    
%         set( gcf, 'Name', 'Minima Imposed at Seed Points' );
%     end
    
    minObjDiameterImsp = min(cellDiameterRange(1) ./ spacing(1:2));    
    cleanStrel = strel('disk', round(0.25 * minObjDiameterImsp));
    imCellMask = imopen(L > 0, cleanStrel);       
    L = bwlabeln( imCellMask );
    
    regLabelWithSeed = unique( L( imCellSeedPoints > 0 ) );
    L( ~ismember(L, regLabelWithSeed) ) = 0; 
    L = bwlabeln(L > 0);
    
    % perform region merging
    if ~isempty(regionMergingModelFile)
        
       fprintf('\n\n>> Correcting over-segmentation errors using a learning-based region-merging approach ... \n\n' );

       mergeTimer = tic;
       
       L = PerformRegionMerging(imAdjusted / max(imInputStandardized(:)), L, imCellSeedPoints, regionMergingModelFile, spacing, flagDebugMode, flagParallelize ); 

       fprintf( '\nRegion Merging took a total of %f seconds\n', toc(mergeTimer) );
       
    end

    % Clean up 
    fprintf( '\n\n>> Cleaning up ... \n\n' );
    
    objRegProps = regionprops(L, 'Area', 'PixelIdxList', 'BoundingBox' );

        % perform sanity checks
        flagPruneRegion = false( numel(objRegProps), 1 );
        for i = 1:numel(objRegProps)

            % prune regions whose bounding box is below a threshold
            bboxSideLength = objRegProps(i).BoundingBox((ndims(imInput)+1):end);
            flagIsCurBBoxBigEnough = all( bboxSideLength >= minBBox );
            if ~flagIsCurBBoxBigEnough
                flagPruneRegion(i) = true;
                continue;
            end
            
            % prune regions with volume below a certain threshold
            if objRegProps(i).Area * prod(spacing) < minCellVolume 
                flagPruneRegion(i) = true;
                continue;
            end
            
            % zero-out regions without seed-point
            curRgnSeedPixInd = objRegProps(i).PixelIdxList( imCellSeedPoints( objRegProps(i).PixelIdxList ) > 0 );
            if isempty(curRgnSeedPixInd)
                flagPruneRegion(i) = true;
                continue;
            end

            % zero-out regions with all seed points with response below a speficied cutoff
            if max( imResponse(curRgnSeedPixInd) ) < seedDetectionResponseCutoff
                flagPruneRegion(i) = true;             
                continue;
            end

            % zero-out regions with max distance map value below a specified cutoff
            if max(imDistMap(objRegProps(i).PixelIdxList)) < seedToBackgroundDistanceCutoff
                flagPruneRegion(i) = true;             
                continue;
            end
            
            % check if there is enough overlap with the specified roi
            if ~isempty(roiMask) && minCellROIOverlap > 0
                
                areaInsideROI = sum( roiMask( objRegProps(i).PixelIdxList ) );
                curRoiOverlap = areaInsideROI / objRegProps(i).Area;
                if curRoiOverlap <= minCellROIOverlap
                    flagPruneRegion(i) = true;      
                    continue;
                end
                
            end
            
        end
        L( ismember(L, find(flagPruneRegion)) ) = 0;
    
        % ignore cells touching the X or Y image border
        if flagIgnoreCellsOnXYBorder
            
            imXYBorder = zeros( size(imInput) );
            indBorder = cell(1, ndims(imInput));
            for i = 1:ndims(imInput)
                indBorder{i} = ':';
            end
            for i = 1:2
                curIndBorder =  indBorder;
                curIndBorder{i} = 1;
                imXYBorder( curIndBorder{:} ) = 1;
                curIndBorder{i} = size(imInput, i);
                imXYBorder( curIndBorder{:} ) = 1;
            end        
            borderObjLabels = unique( L( imXYBorder > 0 ) );
            L( ismember(L, borderObjLabels) ) = 0;
            
        end
        
        % ignore cells lying completely outside the roiMask
        if ~isempty(roiMask)           
            
            insideRoiLabels = unique( L( roiMask ) );
            L( ~ismember(L, insideRoiLabels) ) = 0;
            
        end
        
    % Generate labelled segmentation image
    fprintf( '\n\n>> Generating Labelled Segmentation Mask \n\n' );
    
    imLabelCellSeg = bwlabeln( L > 0 );
    
    [segMaskRGB, cellSegColorMap] = label2rgbND( imLabelCellSeg );
    seedMaskRGB = seedMask;
    seedMaskRGB(:,:,:,2:3) = 0;
    
    if flagDebugMode        
        imWatershedSurface = bwperim(imLabelCellSeg > 0);
        imWatershedSurface(:,:,:,2:3) = 0;
        
        imseriesmaskshowrgb( imInput, {segMaskRGB, seedMaskRGB, imWatershedSurface}, 'maskAlphas', [0.2, 0.2, 0.5], 'spacing', spacing );
        set( gcf, 'Name', 'Final Segmentation Result' );
    end
    
    % display segmentation result in imaris
    if flagDisplayResultsInImaris
        
        imCellSeedPoints(~imLabelCellSeg) = 0;
        DisplayCellSegmentationResultInImaris(imInput, imLabelCellSeg, cellSegColorMap, imCellSeedPoints, ...
                                              'spacing', spacing );
                                              
    end

    % display computation time
    timeElapsed = toc( totalSegTimer );
    fprintf( '\n\n>> Segmentation took %.2f seconds ... \n\n', timeElapsed );
    
    if nargout > 1
        imCellSeedPoints(~L) = 0;
        varargout{1} = imCellSeedPoints .* imResponse;
    end
    
    if nargout > 2
        varargout{2} = p.Results;
    end
    
    if flagParallelize && ~flagPoolOpenedAlready
        matlabpool close;
    end    
    
end

function [imLabelCellSegProcessed] = PerformRegionMerging(imInput, imLabelCellSeg, imCellSeedPoints, wekaModelFile, spacing, flagDebugMode, flagParallelize)

    % load weka model file
    import weka.*;
    import weka.core.*;    
    import weka.core.SerializationHelper.*;
    wekaModel = weka.core.SerializationHelper.readAll( wekaModelFile );
    
    % compute basic cell properties 
    numCells = max( imLabelCellSeg(:) );
    flagIsDegenerateRegion = false(1, numCells);
    
    if flagParallelize
    
        regPropsArr = cell(numCells, 1);
        parfor cellId = 1:numCells

            curRegProps = ComputeRegionProperties(imLabelCellSeg, cellId);

            % add seed point location
            indCurRegionSeedPoint = curRegProps.PixelIdxList( imCellSeedPoints(curRegProps.PixelIdxList) > 0 );
            if isempty( indCurRegionSeedPoint )            
                curRegProps.ptCellSeedPoint = curRegProps.ptCentroid;
            else
                [~, maxind] = max( imCellSeedPoints( indCurRegionSeedPoint ) );
                indCurRegionSeedPoint = indCurRegionSeedPoint( maxind );    

                curRegProps.ptCellSeedPoint = ind2submat( size(imCellSeedPoints), indCurRegionSeedPoint );
            end

            regPropsArr{cellId} = curRegProps;

        end

        cellProps = [];
        for cellId = 1:numCells 
           cellProps = [cellProps; regPropsArr{cellId}];   
           flagIsDegenerateRegion(cellId) = (regPropsArr{cellId}.ConvexArea == 0);
        end
        clear regPropsArr;
        
    else
        
        cellProps = [];
        for cellId = 1:numCells

            curRegProps = ComputeRegionProperties(imLabelCellSeg, cellId);

            % add seed point location
            indCurRegionSeedPoint = curRegProps.PixelIdxList( imCellSeedPoints(curRegProps.PixelIdxList) > 0 );
            if isempty( indCurRegionSeedPoint )            
                curRegProps.ptCellSeedPoint = curRegProps.ptCentroid;
            else
                [~, maxind] = max( imCellSeedPoints( indCurRegionSeedPoint ) );
                indCurRegionSeedPoint = indCurRegionSeedPoint( maxind );    

                curRegProps.ptCellSeedPoint = ind2submat( size(imCellSeedPoints), indCurRegionSeedPoint );
            end

            cellProps = [cellProps; curRegProps];
            flagIsDegenerateRegion(cellId) = (curRegProps.ConvexArea == 0);

        end

    end

    indDegenerateRegions = find(flagIsDegenerateRegion); % these will be dropped
    
    fprintf( '\nTotal number of regions before region merging: %d\n', numel(cellProps) );
    
    % build region adjacency graph
    adjMatrixInitial = zeros(numCells, numCells);

    neighStrel = streldisknd(2*ones(1, ndims(imInput)));

    for i = 1:numCells

       if flagIsDegenerateRegion(i)
           continue;
       end
       
       % get cropped label mask        
       imCurLabelCellSegCropped = imLabelCellSeg( cellProps(i).indCropBBox{:} );
       imCurCellMask = (imCurLabelCellSegCropped == i);

       % check if the cell has significant boundary overlap with any other cells
       imDilatedCellMask = imdilate( imCurCellMask, neighStrel );
       neighPixLabels = imCurLabelCellSegCropped( imDilatedCellMask - imCurCellMask > 0 );
       neighCellIdList = setdiff( unique(neighPixLabels), [0, i] );            
       neighCellIdList( neighCellIdList < i) = [];
       neighCellIdList( flagIsDegenerateRegion(neighCellIdList) ) = [];
       
       % add edges to adjacency matrix
       adjMatrixInitial(i, neighCellIdList) = 1;
       adjMatrixInitial(neighCellIdList, i) = 1;

    end        

    fprintf( '\nTotal number of edges in the region adjacency graph: %d\n', 0.5 * sum(adjMatrixInitial(:)) );

    fprintf( '\nAverage degree of a vertex in the region adjacency graph: %.2f +/- %.2f \n', mean(sum(adjMatrixInitial,2)), std(sum(adjMatrixInitial,2)) );
    
    % initialize the similarity matrix
    simMatrixInitial = adjMatrixInitial;
    [v1, v2, ~] = find( adjMatrixInitial );

    if flagParallelize

       edgeInd = (find(v1 < v2))';
       simFeatures = cell( numel(edgeInd), 1 );
       
       parfor i = 1:numel(edgeInd)

           c1 = v1( edgeInd(i) );
           c2 = v2( edgeInd(i) );
           
           % compute region merging features
           simFeatures{i} = ComputeRegionMergingFeatures(imInput, imLabelCellSeg, cellProps, c1, c2, 'spacing', spacing);          
           
       end
       
       for i = 1:numel(edgeInd)

            c1 = v1( edgeInd(i) );
            c2 = v2( edgeInd(i) );
           
            % compute similarity based on probabilistic prediction of a trained classifier
            simVal = ComputeRegionSimilarity(wekaModel, simFeatures{i});
            
            simMatrixInitial(c1, c2) = simVal;
            simMatrixInitial(c2, c1) = simVal;
            
       end
       
    else

       for k = (find(v1 < v2))'

           featureStruct = ComputeRegionMergingFeatures(imInput, imLabelCellSeg, cellProps, v1(k), v2(k), 'spacing', spacing);          
                      
           simVal = ComputeRegionSimilarity(wekaModel, featureStruct);

           simMatrixInitial(v1(k), v2(k)) = simVal; 
           simMatrixInitial(v2(k), v1(k)) = simVal;

       end

    end

    % perform iterative region merging
    pq = PriorityQueue();       
    imLabelCellSegProcessed = imLabelCellSeg;
    
        % add all max-similarity edges into a priority queue
        [v1, v2, simMatrixVal] = find( simMatrixInitial );        

        for k = 1:numel(v1)

            i = v1(k);
            j = v2(k);

            if i >= j 
                continue;
            end

            if IsEdgeLocallyMaximal(simMatrixInitial, i, j)                               
                pq.insert( simMatrixVal(k), [i, j] );
            end                                
        end

        % merge until queue becomes empty
        simMatrix = simMatrixInitial;
        adjMatrix = adjMatrixInitial;
        
        flagWasVertexMerged = zeros(numCells, 1);
        vertexMergeTreeHeight = zeros(numCells, 1);
        
        numMerges = 0;

        [segMaskRGB, segColorMap] = label2rgbND( imLabelCellSeg );
 
        if flagDebugMode
        
            objSeedMask = imdilate(imCellSeedPoints, ones(5 * ones(1, ndims(imInput))) );
            objSeedMaskRGB = cat( ndims(imInput)+1, objSeedMask, zeros( [size(objSeedMask), 2] ) );

            imseriesmaskshowrgb(imInput, {segMaskRGB, objSeedMaskRGB});
            set( gcf, 'Name', 'Result Before Merging' );

            figMergeProgress = figure;        
            set(gcf, 'Name', 'MergeProgress');
        
        end
        
        mergeProgVec = [];
        
        regProps = cellProps;
        maxRegionInd = numCells;
        
        segColorMapProcessed = segColorMap;
        segColorMapProcessed(end, :) = [];

        while ~pq.isEmpty()

            % find the most similar edge and merge it
            [curEdgeSimilarity, nodeIndices] = pq.pop();

            % check if any of nodes have been merged into another earlier
            if any( flagWasVertexMerged(nodeIndices) )
                continue;
            end

            % create a new region
            maxRegionInd = maxRegionInd + 1;
            curRegionId = maxRegionInd;
            
            flagWasVertexMerged( nodeIndices ) = curRegionId;
            vertexMergeTreeHeight( curRegionId ) = max( vertexMergeTreeHeight(nodeIndices) ) + 1; 
            
            % display regions being merged
            if flagDebugMode
                
                figure(figMergeProgress);                

                mergeProgVec(end+1, :) = [nodeIndices, curRegionId, curEdgeSimilarity];
                plot(mergeProgVec(:,4), 'b-', 'LineWidth', 2.0);
                ylim([0, 1]);
                
            end
            
            % compute merged region 
            mergeStrel = ones(5 * ones(1, ndims(imInput)));
            
                % crop portion of the image enclosing the two cells
                ptBBoxCorners = [ regProps(nodeIndices(1)).BoundingBox, regProps(nodeIndices(2)).BoundingBox ];
                indCropBox = cell(1, ndims(imInput));
                for i = 1:ndims(imInput)
                    indCropBox{i} = min(ptBBoxCorners(i,:)):max(ptBBoxCorners(i,:));        
                end
                imCurLabelCellSegCropped = imLabelCellSegProcessed( indCropBox{:} );                
                imCurComponentCellMaskUnion = ismember(imCurLabelCellSegCropped, nodeIndices);
                imCurWatershed = imclose(imCurComponentCellMaskUnion, mergeStrel) & ~imCurComponentCellMaskUnion & ...
                                 imdilate(imCurLabelCellSegCropped == nodeIndices(1), mergeStrel) & ...
                                 imdilate(imCurLabelCellSegCropped == nodeIndices(2), mergeStrel);                
                imCurMergedRegionMask = imCurComponentCellMaskUnion | imCurWatershed;
                imCurLabelCellSegCropped( imCurMergedRegionMask > 0 ) = curRegionId;                
                imLabelCellSegProcessed( indCropBox{:} ) = imCurLabelCellSegCropped;
            
            segColorMapProcessed(curRegionId, :) = segColorMapProcessed(nodeIndices(2), :);
            
            % compute region properties
            curRegProps = ComputeRegionProperties(imLabelCellSegProcessed, curRegionId, spacing);

            indCurRegionSeedPoint = curRegProps.PixelIdxList( imCellSeedPoints(curRegProps.PixelIdxList) > 0 );
            if isempty( indCurRegionSeedPoint )            
                ptCurRegionSeedPoint = curRegProps.Centroid;
            else
                [~, maxind] = max( imCellSeedPoints( indCurRegionSeedPoint ) );
                indCurRegionSeedPoint = indCurRegionSeedPoint( maxind );    
                ptCurRegionSeedPoint = ind2submat( size(imCellSeedPoints), indCurRegionSeedPoint );
            end
            curRegProps.ptCellSeedPoint = ptCurRegionSeedPoint;
            regProps(curRegionId) = curRegProps;

            % add a node in the graph for the merged region and compute similarity to adjacent regions
            neighRegionIndices = setdiff( find(adjMatrix(nodeIndices(1) , :) > 0 | adjMatrix(nodeIndices(2) , :) > 0), nodeIndices );

            simMatrix(curRegionId, :) = 0;
            simMatrix(:, curRegionId) = 0;
            
            adjMatrix(curRegionId, :) = 0;
            adjMatrix(:, curRegionId) = 0;
            
            flagWasVertexMerged(curRegionId) = 0;
            
            for nid = neighRegionIndices

                if flagWasVertexMerged(nid)
                    error('error: code shouldnt come here');
                end

                adjMatrix(curRegionId, nid) = 1;
                adjMatrix(nid, curRegionId) = 1;
                
                featureStruct = ComputeRegionMergingFeatures(imInput, imLabelCellSegProcessed, regProps, curRegionId, nid, 'spacing', spacing);          
                curNeighSimVal = ComputeRegionSimilarity(wekaModel, featureStruct);
                
                if curNeighSimVal > 0
                    simMatrix(curRegionId, nid) = curNeighSimVal;
                    simMatrix(nid, curRegionId) = curNeighSimVal;
                end

            end

            % remove the nodes and all associated edges from the graph
            simMatrix(nodeIndices, :) = 0; 
            simMatrix(:, nodeIndices) = 0;        
            
            adjMatrix(nodeIndices, :) = 0;
            adjMatrix(:, nodeIndices) = 0;
            
            % add locally maximal edges to priority queue    
            for nid = neighRegionIndices

                if IsEdgeLocallyMaximal(simMatrix, curRegionId, nid)                               
                    pq.insert( full(simMatrix(curRegionId, nid)), sort([curRegionId, nid]) );
                end                                

            end

            numMerges = numMerges + 1;

        end           
        
        numMerges         
        fprintf( '\nTotal number of regions after region merging: %d\n', numel( unique(imLabelCellSegProcessed(:)) ) );

        if flagDebugMode
            
            segColorMapProcessed(end+1,:) = [0, 0, 0];
            imseriesmaskshowrgb(imInput, {label2rgbND(imLabelCellSegProcessed, segColorMapProcessed), objSeedMaskRGB});
            set(gcf, 'Name', 'Result After Merging');
            
        end

        imLabelCellSegProcessed = bwlabeln(imopen(imLabelCellSegProcessed > 0, ones(3*ones(1,ndims(imInput))) ));
        
end


function [ simVal ] = ComputeRegionSimilarity(wekaModel, featureStruct)

    % convert feature structure to a linear feature vector
    [featureVec , featureNameList] = ConvertFeatureStructToFeatureVec( featureStruct );        

    % apply region merging classifier
    [predictedClassLabel, predictionProbability] = ApplyWekaRegionMergingClassifier(wekaModel, featureVec, featureNameList);

    if predictedClassLabel > 0
        simVal = predictionProbability;
    else
        simVal = 0; % this edge will be skipped by the merging process
    end
    
end

function [flagMax] = IsEdgeLocallyMaximal( similarityMatrix, v1, v2 )

    matDim = size(similarityMatrix, 1);
    
    if similarityMatrix(v1, v2) <= 0 || ...
       any( similarityMatrix(v1, [1:v2-1, v2+1:matDim]) > similarityMatrix(v1, v2) ) || ...
       any( similarityMatrix(v2, [1:v1-1, v1+1:matDim]) > similarityMatrix(v1, v2) )            

        flagMax = false;

    else
        
        flagMax = true;
        
    end

end
