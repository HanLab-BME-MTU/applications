function test_overseg_merge_2d

clc
clear 
close all

%***********************************************************************
% 
%                           PARAMETERS 
%
%***********************************************************************

    dataFileDir = 'C:\deepak\data\Stefan_July2012';
    
    %SEG_ALGO = 'RANDOM_WALKS';
    SEG_ALGO = 'WATERSHED';

    % critical parameters (needs tuning)
    cellDiameterRange = [8, 20];
    mergeSimilarityThreshold = 0.8;
    flagUseConvexity = false;
    logResponseCutoff = 2.0;
    flagLearnThreshold = true;
    
    % less-critical parameters (leave defaults)
    threshNeighOverlap = 5;
    neighSearchRad = 2;    
    numHistogramBins = 16.0;
    
%***********************************************************************

% load data
if ~exist( 'dataFilePath', 'var' )
    [fileName,pathName] = uigetfile( fullfile( dataFileDir, '*.mat' ), 'Select the data file' );   
    dataFilePath = fullfile( pathName, fileName )
end

dataFileContents = load( dataFilePath );
imInput = dataFileContents.im;
spacing = [0.5, 0.5];

% for parallelization
if matlabpool( 'size' ) == 0
    matlabpool open;
end            

% adjust the intensity range of the image 
imageDynamicRange = ComputeImageDynamicRange( imInput, 99.9 );
imAdjusted = mat2gray(imInput, imageDynamicRange) * 2000;    

% apply median filter to remove any spotty noise
imAdjusted = medfilt2( imAdjusted, [3,3] );    

% compute threshold
fprintf( '\nThresholding the Image using a Global Thresholding Algorithm ...\n' );

localWindowRadius = round(30 ./ spacing(1));
localWindowPace = round(localWindowRadius / 3);
minLocalGlobalThresholdRatio = 0.6;

imThresh = segmentCellForegroundUsingLocalMinError( imAdjusted, localWindowRadius, ...
                                                    'model', 'poisson', ...  
                                                    'localWindowPace', localWindowPace, ...
                                                    'minLocalGlobalThresholdRatio', minLocalGlobalThresholdRatio, ...
                                                    'flagParallelize', true );

                                                
imseriesmaskshow( imAdjusted, imThresh );
set( gcf, 'Name', 'Result of Thresholding' );

if false % strcmp(SEG_ALGO, 'WATERSHED')
    imAdjusted( ~imThresh ) = 0;
end

% compute cell seed points
imDistMap = callMatitkFilter( 'FSignedMaurerDistancemap', {} , imThresh, [], [], spacing  );
imDistMap( imDistMap > 0 ) = 0;
imDistMap = -1 * imDistMap;

imseriesmaskshow( imDistMap, imThresh );
set( gcf, 'Name', 'Foreground Distance Map' );    

numLoGScales = min(20, cellDiameterRange(2)-cellDiameterRange(1)+1);
[imCellSeedPoints, imResponse] = detectBlobsUsingAdaptiveMultiscaleLoG( imAdjusted, ...
                                                                        imDistMap, ... 
                                                                        'blobDiameterRange', cellDiameterRange, ...
                                                                        'numLoGScales', numLoGScales, ...
                                                                        'spacing', spacing, ...
                                                                        'debugMode', true );
                                                                    
imCellSeedPoints(~imThresh) = 0; % remove seed points in background
meanResponseBgnd = mean( imResponse(~imThresh) );
stdResponseBgnd = std( imResponse(~imThresh) );
imCellSeedPoints( imCellSeedPoints < (meanResponseBgnd + logResponseCutoff * stdResponseBgnd) ) = 0; % remove seed points with weak LoG response

objSeedMask = double(imdilate(imCellSeedPoints > 0, streldisknd([3,3])));
objSeedMaskRGB = cat(3, objSeedMask, zeros([size(objSeedMask), 2]));

imseriesmaskshow( imAdjusted, objSeedMask );
set( gcf, 'Name', 'Cell Seed Point Detection Result' );

% Compute image gradient
imGradMag = callMatitkFilter( 'FGradientMagnitude', {}, imAdjusted );

% Apply marker-based Watershed
fprintf( '\n\n>> Running segmentation algorithm ...\n\n' );

switch SEG_ALGO

    case 'WATERSHED'

        imBackgroundSeeds = bwmorph( ~imThresh, 'shrink', 1 );

        imFeature = 1 - mat2gray( imResponse );
        %imFeature = mat2gray( imGradMag );

        imMinimaImposedAtSeeds = imimposemin( imFeature, imCellSeedPoints | imBackgroundSeeds );
        imMinimaImposedAtSeedsDisplay = imMinimaImposedAtSeeds;
        imMinimaImposedAtSeedsDisplay( imMinimaImposedAtSeeds == -Inf ) = min( imMinimaImposedAtSeeds( imMinimaImposedAtSeeds ~= -Inf ) ) - 3 * eps;

        imseriesmaskshow( imMinimaImposedAtSeedsDisplay, {objSeedMask, imBackgroundSeeds}, 'spacing', spacing );    
        set( gcf, 'Name', 'Minima Imposed at Seed Points' );

        L = watershed( imMinimaImposedAtSeeds );
        objLabelInd = unique( L( imCellSeedPoints > 0 ) );
        L( ~ismember(L, objLabelInd) ) = 0;

        % Post-processing - remove very small objects
        minCellVolume = pi * (0.5 * 0.5 * min(cellDiameterRange))^2;
        regstats = regionprops(L, 'Area');
        objAreaList = [regstats.Area] * prod(spacing);
        smallObjInd = find( objAreaList < minCellVolume );
        L( ismember( L, smallObjInd ) ) = 0;
        imLabelCellSeg = bwlabeln( L > 0 );
        
    case 'RANDOM_WALKS'
        
        imBackgroundSeeds = bwmorph( ~imThresh, 'shrink', 2 );
        bgndSeedIndices = find(imBackgroundSeeds);
        
        imForegroundSeeds = bwmorph(imCellSeedPoints > 0, 'thicken', 5);
        L = bwlabeln( imForegroundSeeds );
        fgndSeedIndices = find(L > 0);

        seedIndices = [fgndSeedIndices; bgndSeedIndices];
        seedLabels = [L(fgndSeedIndices); zeros(size(bgndSeedIndices))];
        numLabels = max(L(:)) + 1;
        
        imseriesmaskshow( imAdjusted, {imBackgroundSeeds, imForegroundSeeds} );
        set( gcf, 'Name', 'Foreground-background seeds' );
        
%         [res, theta, nmsEdge, scaleMap] = multiscaleSteerableDetector(mat2gray(imAdjusted), 3, [1, 2, 4]);
%         [res, theta, nmsRidge, scaleMap] = multiscaleSteerableDetector(1 - mat2gray(imAdjusted), 2, [1, 2, 4]);        
%         imFeature = mat2gray( (1 - mat2gray(imAdjusted)) .* (mat2gray(nmsEdge) + mat2gray(nmsRidge)) );

        [res, theta, nmsEdge, scaleMap] = multiscaleSteerableDetector(mat2gray(imAdjusted), 3, [1, 2, 4]);
        %imFeature = mat2gray( mat2gray( (1 - mat2gray(imAdjusted)) ) .* mat2gray(nmsEdge) );
        imFeature = mat2gray(nmsEdge);
        
        [L, ...
         pixToLabelProbabilityMap] = random_walker_segmentation( imFeature, seedIndices, seedLabels, ...
                                                                 'sigma', 0.05, ...
                                                                 'flagNormalize', false);
        
        imLabelCellSeg = L;
        for i = 1:max(L(:))      
        
            imCurCellMask = ( imLabelCellSeg == i );
            imCurCellPerim = imdilate(imCurCellMask, streldisknd(1 * ones(1,ndims(imInput)))) & ~imCurCellMask;
            imLabelCellSeg( imCurCellPerim > 0 ) = 0;
            
        end
        
end

[segMaskRGB, segColorMap] = label2rgbND( imLabelCellSeg );

[gx, gy] = gradient(imLabelCellSeg); 
imSegBoundaryMask = (gx ~= 0 | gy ~= 0);

imseriesmaskshow(imInput, {imSegBoundaryMask, objSeedMask});
set(gcf, 'Name', 'cell segmentation boundary mask overlay');

% Post-processing - correct oversegmentation 
numCells = max(imLabelCellSeg(:));    
adjMatrix = sparse(numCells, numCells);

    % get some stats for each cell region
    cellStats = [];    
    for i = 1:numCells                       
        cellStats = [ cellStats; ComputeRegionProperties( imAdjusted, imLabelCellSeg, imCellSeedPoints, i, spacing ) ];        
    end
    
    % build cell neighborhood graph    
    displayAdjacencyGraphX = [];
    displayAdjacencyGraphY = [];    
    neighStrel = streldisknd(neighSearchRad*ones(1, ndims(imInput)));
    
    for i = 1:numCells               
        
        % get cropped label mask        
        imCurLabelCellSegCropped = imLabelCellSeg( cellStats(i).indCropBBox{:} );
        imCurCellMask = (imCurLabelCellSegCropped == i);
        
        % check if it has any neighboring cells
        imDilatedCellMask = imdilate( imCurCellMask, neighStrel );
        neighPixLabels = imCurLabelCellSegCropped( imDilatedCellMask - imCurCellMask > 0 );
        neighCellIdList = setdiff( unique(neighPixLabels), [0, i] );
        
        if ~isempty(neighCellIdList)
            
            for j = neighCellIdList

                if i >= j
                    continue;
                end
                
                % first check overlap is significant
                numOverlapPixels = numel( find(neighPixLabels == j) );                
                percentOverlap = 100.0 * numOverlapPixels / min( cellStats(i).Perimeter, cellStats(j).Perimeter );
                
                if percentOverlap < threshNeighOverlap
                   continue;
                end
                    
                adjMatrix(i, j) = 1;
                adjMatrix(j, i) = 1;                
                
                displayAdjacencyGraphX(:, end+1 ) = [cellStats(i).ptCellSeedPoint(2); cellStats(j).ptCellSeedPoint(2)];
                displayAdjacencyGraphY(:, end+1 ) = [cellStats(i).ptCellSeedPoint(1); cellStats(j).ptCellSeedPoint(1)];
                
            end
            
        end
        
    end  

    imseriesmaskshowrgb(imInput, {segMaskRGB, objSeedMaskRGB});
    set(gcf, 'Name', 'cell neighborhood graph overlayed on gradient magnitude');
    hold on;
    plot( displayAdjacencyGraphX, displayAdjacencyGraphY, 'go-', 'LineWidth', 2.0 );
    hold off;
    
    imseriesmaskshowrgb(imGradMag, {segMaskRGB, objSeedMaskRGB})
    hold on;
    plot( displayAdjacencyGraphX, displayAdjacencyGraphY, 'go-', 'LineWidth', 2.0 );
    hold off;

    imseriesmaskshow(imLabelCellSeg, {imSegBoundaryMask, objSeedMask});
    set(gcf, 'Name', 'segmentation label map');
    hold on;
    plot( displayAdjacencyGraphX, displayAdjacencyGraphY, 'go-', 'LineWidth', 2.0 );
    hold off;
    
    % perform merging
    simMatrixConvexity = adjMatrix;    
    simMatrixWatershedSaliency = adjMatrix;
    simMatrixAppearanceHistogram = adjMatrix; 
    
        % initialize similarity matrices
        [v1, v2, ~] = find( adjMatrix );
        
        for k = 1:numel(v1)
            
            i = v1(k);
            j = v2(k);
            
            if i >= j 
                continue;
            end

            [curSimApp, curSimWatershed, curSimShape] = ComputeRegionSimilarity( imAdjusted, imLabelCellSeg, cellStats, i, j, neighStrel, numHistogramBins );

            simMatrixAppearanceHistogram(i, j) = curSimApp;
            simMatrixAppearanceHistogram(j, i) = curSimApp;
            
            simMatrixWatershedSaliency(i, j) = curSimWatershed; % higher the better
            simMatrixWatershedSaliency(j, i) = curSimWatershed;
            
            simMatrixConvexity(i, j) = curSimShape; % higher the better
            simMatrixConvexity(j, i) = curSimShape;
            
        end

        [v1, v2, ~] = find( adjMatrix );        
        [~, ~, simApp] = find( simMatrixAppearanceHistogram );
        [~, ~, simWatershed] = find( simMatrixWatershedSaliency );
        [~, ~, simConvexity] = find( simMatrixConvexity );

        if flagUseConvexity
            simCombined = min( [simApp, simWatershed, simConvexity], [], 2 );        
        else
            simCombined = min( [simApp, simWatershed], [], 2 );        
        end
        
        simCombinedMax = max( simCombined );
        simCombinedMin = min( simCombined );
        simMatrixCombined = sparse(v1, v2, simCombined, numCells, numCells);

        % display similarity matrices for all criterion
        imseriesmaskshow(imAdjusted, {imSegBoundaryMask, objSeedMask});        
        DisplayWeightedRegionAdjacencyGraph( gca, simMatrixAppearanceHistogram, cellStats );
        set( gcf, 'Name', 'Appearance Similarity' );

        imseriesmaskshow(imAdjusted, {imSegBoundaryMask, objSeedMask});        
        DisplayWeightedRegionAdjacencyGraph( gca, simMatrixWatershedSaliency, cellStats );
        set( gcf, 'Name', 'Watershed Saliency' );

        imseriesmaskshow(imAdjusted, {imSegBoundaryMask, objSeedMask});        
        DisplayWeightedRegionAdjacencyGraph( gca, simMatrixConvexity, cellStats );
        set( gcf, 'Name', 'Convexity' );

        imseriesmaskshow(imAdjusted, {imSegBoundaryMask, objSeedMask});        
        DisplayWeightedRegionAdjacencyGraph( gca, simMatrixCombined, cellStats );
        set( gcf, 'Name', 'Combined' );
        
        % iteratively merge edges that are locally maximal
        pq = PriorityQueue();       
        imLabelCellSegProcessed = imLabelCellSeg;
        regStats = cellStats;
        
            % add all max-similarity edges into a priority queue
            [v1, v2, simCombined] = find( simMatrixCombined );        
            
            for k = 1:numel(v1)
                
                i = v1(k);
                j = v2(k);

                if i >= j 
                    continue;
                end
                
                if IsEdgeLocallyMaximal(simMatrixCombined, i, j)                               
                    pq.insert( simCombined(k), [i, j] );
                end                                
            end

            % merge until queue becomes empty
            flagWasVertexMerged = zeros(numCells, 1);
            simMatrixCombinedProcessed = simMatrixCombined;
            numMerges = 0;
            mergeProgSimVec = [];
            
            imseriesmaskshowrgb(imInput, {segMaskRGB, objSeedMaskRGB});
            set( gcf, 'Name', 'Result Before Merging' );
            
            figMergeProgress = figure;
            set(gcf, 'Name', 'MergeProgress');
            subplot(1,2,1)
            [ imMaskOverlay ] = genImageRGBMaskOverlay( imAdjusted, {segMaskRGB, objSeedMaskRGB}, [0.5, 0.5] );
            imshow( imMaskOverlay );
            
            jetColorMap = colormap( 'jet' );
            numColors = size(jetColorMap, 1);
            
            displayeMergedEdgesX = [];
            displayeMergedEdgesY = [];
            
            while ~pq.isEmpty()

                % find the most similar edge and merge it
                [curEdgeSimilarity, nodeIndices] = pq.pop();

                i = nodeIndices(1);
                j = nodeIndices(2);
                
                % check if any of nodes have been merged into another earlier
                if any( flagWasVertexMerged(nodeIndices) ) || simMatrixCombinedProcessed(i,j) ~= curEdgeSimilarity
                    continue;
                end
                
                % check if it is ok to merge these vertices
                if curEdgeSimilarity < mergeSimilarityThreshold
                    continue;                    
                end
                
                % display regions being merged
                figure(figMergeProgress);                
                curEdgeSimilarityNmzd = mat2gray(curEdgeSimilarity, [simCombinedMin, simCombinedMax]);
                subplot(1, 2, 1)
                hold on;                
                lineColor = jetColorMap( 1 + round(curEdgeSimilarityNmzd * (numColors-1)), : );
                plot( [regStats(i).ptCellSeedPoint(2), regStats(j).ptCellSeedPoint(2)], ...
                      [regStats(i).ptCellSeedPoint(1), regStats(j).ptCellSeedPoint(1)], ...
                      '-', 'Color', lineColor, 'LineWidth', 1.0 + 8.0 * curEdgeSimilarityNmzd );
                hold off;               
                
                mergeProgSimVec(end+1) = curEdgeSimilarityNmzd;
                subplot(1, 2, 2), plot(mergeProgSimVec, 'b-', 'LineWidth', 2.0);
                                
                mergeMovie(numMerges+1) = getframe(figMergeProgress);
                
                displayeMergedEdgesX(:, end+1) = [regStats(i).ptCellSeedPoint(2), regStats(j).ptCellSeedPoint(2)];
                displayeMergedEdgesY(:, end+1) = [regStats(i).ptCellSeedPoint(1), regStats(j).ptCellSeedPoint(1)];
                
                % do merging
                if regStats(i).Area > regStats(j).Area
                    pid = i;
                    chid = j;
                else
                    pid = j;
                    chid = i;
                end
                
                flagWasVertexMerged( chid ) = pid;
                simMatrixCombinedProcessed(chid, :) = 0;
                simMatrixCombinedProcessed(:, chid) = 0;
                
                imRegionUnionMask = ismember( imLabelCellSegProcessed, [i, j] );
                imRegionUnionMask = imclose(imRegionUnionMask, neighStrel);            
                imLabelCellSegProcessed( imRegionUnionMask > 0 ) = pid;
                
                regStats(pid) = ComputeRegionProperties( imAdjusted, imLabelCellSegProcessed, imCellSeedPoints, pid, spacing );
                
                neighRegionIndices = setdiff( find(simMatrixCombinedProcessed(i , :) > 0 | simMatrixCombinedProcessed(j , :) > 0), [i,j] );
                
                for nid = neighRegionIndices
                
                    if flagWasVertexMerged(nid)
                        error('error: code shouldnt come here');
                    end
                    
                    [curSimApp, curSimWatershed, curSimShape] = ComputeRegionSimilarity( imAdjusted, imLabelCellSegProcessed, regStats, pid, nid, neighStrel, numHistogramBins );
                    
                    if flagUseConvexity
                        curSimCombined = min( [curSimApp, curSimWatershed, curSimShape ] );
                    else
                        curSimCombined = min( [curSimApp, curSimWatershed] );
                    end
                    curSimCombined = curSimCombined;
                    
                    simMatrixCombinedProcessed( pid, nid ) = curSimCombined;
                    simMatrixCombinedProcessed( nid, pid ) = curSimCombined;
                    
                end
                
                for nid = neighRegionIndices
                    
                    if IsEdgeLocallyMaximal(simMatrixCombinedProcessed, pid, nid)                               
                        pq.insert( full(simMatrixCombinedProcessed(pid, nid)), sort([pid, nid]) );
                    end                                
                    
                end
                
                numMerges = numMerges + 1
                
            end
            
            imseriesmaskshowrgb(imInput, {label2rgbND(imLabelCellSegProcessed, segColorMap), objSeedMaskRGB});
            set(gcf, 'Name', 'Result After Merging');
            hold on;
            plot( displayeMergedEdgesX, displayeMergedEdgesY, 'g-', 'LineWidth', 3.0 );
            hold off;
            
            movie2avi( mergeMovie, fullfile( 'C:\deepak\results', 'RegionMerging.avi' ), 'FPS', 5 );
end

function [regProps] = ComputeRegionProperties( imInput, imLabelCellSeg, imCellSeedPoints, cellId, spacing )

    regProps.PixelIdxList = find(imLabelCellSeg == cellId); 

    regProps.Area = numel( regProps.PixelIdxList );
    
    % get list of coordinates of cell pixels
    ptPixelLocations = cell(1, ndims(imInput));
    [ptPixelLocations{:}] = ind2sub( size(imInput), regProps.PixelIdxList );
    ptPixelLocations = cell2mat( ptPixelLocations );    
    regProps.ptPixelLocations = ptPixelLocations;
    
    regProps.Centroid = round( mean( ptPixelLocations ) );
    
    % compute region bounding box and construct indices to crop region bounding box from any image  
    bboxPadding = 4;
    regProps.BoundingBox = ([ min(ptPixelLocations) ; max(ptPixelLocations) ])';
    
    indCropBBox = cell(1,ndims(imInput));
    for j = 1:ndims(imInput)
        indStart = max(1, regProps.BoundingBox(j,1) - bboxPadding);
        indEnd = min( size(imInput,j), regProps.BoundingBox(j,2) + bboxPadding );
        indCropBBox{j} = indStart:indEnd;
        regProps.BoundingBox(j, :) = [indStart, indEnd];        
    end
    regProps.indCropBBox = indCropBBox;
    
    % get list of coordinates of pixels on region boundary
    regProps.BoundaryPixelIdxList = find( bwperim(imLabelCellSeg == cellId) ); % (can be faster)
    
    ptBoundaryPixelLocations = cell(1, ndims(imInput));
    [ptBoundaryPixelLocations{:}] = ind2sub( size(imInput), regProps.BoundaryPixelIdxList );
    ptBoundaryPixelLocations = cell2mat( ptBoundaryPixelLocations );    

    regProps.ptBoundaryPixelLocations = ptBoundaryPixelLocations;
    regProps.Perimeter = numel( regProps.BoundaryPixelIdxList );
    
    % get list of coordinates of pixels on cell convex hull
    [K, v] = convhulln( ptBoundaryPixelLocations(:, [2, 1, 3:ndims(imInput)]) );

    regProps.ConvexArea = v;          
    regProps.ptConvexHull = ptBoundaryPixelLocations(K, :);
    
    % contrast
    meanIntensityInterior = mean( imInput( regProps.PixelIdxList ) );
    regProps.Contrast = std( imInput(regProps.BoundaryPixelIdxList)  / meanIntensityInterior );

    % ellipticity
    regProps.Ellipticity = ComputeEllipticVariance( regProps.ptPixelLocations, spacing );

    % seed point location
    indCurRegionSeedPoint = find( imCellSeedPoints(regProps.PixelIdxList) > 0 );
    indCurRegionSeedPoint = regProps.PixelIdxList( indCurRegionSeedPoint );
    if isempty( indCurRegionSeedPoint )            
        ptCurRegionSeedPoint = regProps.Centroid;
    else
        
        [maxval, maxind] = max( imCellSeedPoints( indCurRegionSeedPoint ) );
        indCurRegionSeedPoint = indCurRegionSeedPoint( maxind );    
        
        ptCurRegionSeedPoint = cell(1, ndims(imCellSeedPoints));
        [ptCurRegionSeedPoint{:}] = ind2sub( size(imCellSeedPoints), indCurRegionSeedPoint );
        ptCurRegionSeedPoint = cell2mat( ptCurRegionSeedPoint );
        
    end
    regProps.ptCellSeedPoint = ptCurRegionSeedPoint;
    
end

function DisplayWeightedRegionAdjacencyGraph( hAxis, simMatrix, regProps, flagShowMaximallySimilarEdgesOnly, simThresh )

    if ~exist( 'flagShowMaximallySimilarEdgesOnly', 'var' )
        flagShowMaximallySimilarEdgesOnly = false;
    end
    
    axes(hAxis);

    maxLineWidth = 6.0;
    
    [v1, v2, simVal] = find( simMatrix );        
    simValNmzd = mat2gray(simVal);

    jetColorMap = colormap( 'jet' );
    numColors = size(jetColorMap, 1);
    
    for k = 1:numel(v1)

        i = v1(k);
        j = v2(k);

        if i >= j 
            continue;
        end

        if flagShowMaximallySimilarEdgesOnly && ~IsEdgeLocallyMaximal(simMatrix, i, j)           

            hold on;
            plot( [regProps(i).ptCellSeedPoint(2), regProps(j).ptCellSeedPoint(2)], ...
                  [regProps(i).ptCellSeedPoint(1), regProps(j).ptCellSeedPoint(1)], ...
                  '-', 'Color', 'w', 'LineWidth', 1.0 + maxLineWidth * simValNmzd(k) );
            hold off;

        else
            
            if exist( 'simThresh', 'var' ) && simValNmzd(k) < simThresh
                lineColor = 1 - simValNmzd(k) * ones(1,3);                
            else
                lineColor = jetColorMap( 1 + round(simValNmzd(k) * (numColors-1)), : );
            end

            hold on;
            plot( [regProps(i).ptCellSeedPoint(2), regProps(j).ptCellSeedPoint(2)], ...
                  [regProps(i).ptCellSeedPoint(1), regProps(j).ptCellSeedPoint(1)], ...
                  '-', 'Color', lineColor, 'LineWidth', 1.0 + maxLineWidth * simValNmzd(k) );
            hold off;
            
        end

    end

end 

function [simApp, simWatershed, simShape] = ComputeRegionSimilarity( imAdjusted, imLabelCellSeg, regProps, i, j, neighStrel, numHistogramBins )

    % appearance -- histogram      
    binSize = range(imAdjusted(:)) / numHistogramBins;        
    imQuantized = floor( imAdjusted / binSize ) + 1;
    imQuantized( imQuantized > numHistogramBins ) = numHistogramBins;
    
    hist1 = tabulate( imQuantized( regProps(i).PixelIdxList ) );
    hist1 = hist1(:,3) / 100;
    hist1( numel(hist1)+1 : numHistogramBins ) = 0;            

    hist2 = tabulate( imQuantized( regProps(j).PixelIdxList ) );
    hist2 = hist2(:,3) / 100;
    hist2( numel(hist2)+1 : numHistogramBins ) = 0;

    simApp = eps + dot( sqrt( hist1 ), sqrt( hist2 ) );
    
    % contrast between cell interior and boundary
    ptBBoxCorners = [ regProps(i).BoundingBox, regProps(j).BoundingBox ];
    cropind = { min(ptBBoxCorners(1,:)):max(ptBBoxCorners(1,:)), min(ptBBoxCorners(2,:)):max(ptBBoxCorners(2,:)) };
    imLabelCellSegCropped = imLabelCellSeg( cropind{:} );
    imInputCropped = imAdjusted( cropind{:} );
    
    imComponentCellMaskUnion = ismember( imLabelCellSegCropped, [i, j] );
    imMergedCellMask = imclose(imComponentCellMaskUnion, neighStrel);            
    L = bwlabeln(imMergedCellMask & ~imComponentCellMaskUnion); 
    stats = regionprops( L, 'Area' );
    [maxval, maxind] = max( [stats.Area] );
    imWatershed = imdilate(L == maxind, neighStrel);
    
    medianCellIntensity1 = median( imInputCropped( imLabelCellSegCropped == i ) );
    medianCellIntensity2 = median( imInputCropped( imLabelCellSegCropped == j ) );
    medianCellIntensity = median( imInputCropped( imComponentCellMaskUnion > 0 ) );
    medianWatershedIntensity = median( imInputCropped(imWatershed > 0) );

    simWatershed = eps + min(medianWatershedIntensity / medianCellIntensity1, medianWatershedIntensity / medianCellIntensity2);
    %simWatershed = eps + medianWatershedIntensity / medianCellIntensity;
    
    % shape/convexity 
    c1 = regProps(i).Area / regProps(i).ConvexArea;
    c2 = regProps(j).Area / regProps(j).ConvexArea;     % ispinesib       

    ptMergedCellHull = [ regProps(i).ptConvexHull; regProps(j).ptConvexHull ];
    ptMergedCellHull = ptMergedCellHull(:, [2, 1, 3:ndims(imAdjusted)]);
    [~, v] = convhulln( ptMergedCellHull );
    c12 = (regProps(i).Area + regProps(j).Area) / v;            

    simShape = eps + c12 / max([c1, c2]);
   
end

function [flagMax] = IsEdgeLocallyMaximal( similarityMatrix, v1, v2 )

    matDim = size(similarityMatrix, 1);
    
    if any( similarityMatrix(v1, [1:v2-1, v2+1:matDim]) > similarityMatrix(v1, v2) ) || ...
       any( similarityMatrix(v2, [1:v1-1, v1+1:matDim]) > similarityMatrix(v1, v2) )            

        flagMax = false;

    else
        
        flagMax = true;
        
    end

end