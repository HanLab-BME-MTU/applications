function test_overseg_merge_2d_criterion

clc
clear 
close all

%***********************************************************************
% 
%                           PARAMETERS 
%
%***********************************************************************

    dataFileDir = 'C:\deepak\data\Stefan_July2012';

    SEG_ALGO = 'RANDOM_WALKS';
    %SEG_ALGO = 'WATERSHED';
    
    threshNeighOverlap = 5;
    neighSearchRad = 2;
    cellDiameterRange = [5, 20];

    
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

if strcmp(SEG_ALGO, 'WATERSHED')
    imAdjusted( ~imThresh ) = 0;
end

% compute cell seed points
imDistMap = callMatitkFilter( 'FSignedMaurerDistancemap',{} , imThresh, [], [], spacing  );
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
imCellSeedPoints( imCellSeedPoints < (meanResponseBgnd + 2.0 * stdResponseBgnd) ) = 0; % remove seed points with weak LoG response

% meanResponseBgnd = mean( imResponse(~imThresh & imCellSeedPoints > 0) );
% stdResponseBgnd = std( imResponse(~imThresh & imCellSeedPoints > 0) );
%imCellSeedPoints( imCellSeedPoints < (meanResponseBgnd + 2.0 * stdResponseBgnd) ) = 0; % remove seed points with weak LoG response
% imCellSeedPoints(~imThresh) = 0; % remove seed points in background

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
        imFeature = mat2gray( mat2gray( (1 - mat2gray(imAdjusted)) ) .* mat2gray(nmsEdge) );

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

segMaskRGB = label2rgbND( imLabelCellSeg );

% Compute cell stats
cellStats = regionprops(imLabelCellSeg, {'Area', 'BoundingBox', 'Centroid', 'PixelIdxList'} );
numCells = numel(cellStats);

% Post-processing - correct oversegmentation using graph matching
adjMatrix = sparse(numCells, numCells);

    % get additional stats for each cell
    for i = 1:numCells               
        
        % seed point location
        indCurCellSeedPoint = find( imCellSeedPoints(cellStats(i).PixelIdxList) > 0 );
        indCurCellSeedPoint = cellStats(i).PixelIdxList( indCurCellSeedPoint );
        if isempty( indCurCellSeedPoint )            
            ptCurCellSeedPoint = round(cellStats(i).Centroid([2,1,3:ndims(imInput)]));
        else
            ptCurCellSeedPoint = cell(1, ndims(imCellSeedPoints));
            [ptCurCellSeedPoint{:}] = ind2sub( size(imCellSeedPoints), indCurCellSeedPoint(1) );
            ptCurCellSeedPoint = cell2mat( ptCurCellSeedPoint );
        end
        cellStats(i).ptCellSeedPoint = ptCurCellSeedPoint;
        
        % reorder elements in the bounding box from [x, y, z, ...] to [y, x, z, ...]
        curCellBoundingBox = round(cellStats(i).BoundingBox([[2,1,3:ndims(imInput)], ndims(imInput)+[2,1,3:ndims(imInput)]]));
        
        cellStats(i).BoundingBox = curCellBoundingBox;
        
        % construct indices to crop cell bounding box from any image
        bboxPadding = 4;
        
        indCropBBox = cell(1,ndims(imInput));
        for j = 1:ndims(imInput)
            indStart = max(1, curCellBoundingBox(j) - bboxPadding);
            indEnd = min( size(imInput,j), curCellBoundingBox(j) + curCellBoundingBox(ndims(imInput) + j) - 1 + bboxPadding );
            indCropBBox{j} = indStart:indEnd;
        end
        cellStats(i).indCropBBox = indCropBBox;
        
        % get list of coordinates of cell pixels
        ptCellPixels = cell(1, ndims(imInput));
        [ptCellPixels{:}] = ind2sub( size(imInput), cellStats(i).PixelIdxList );
        ptCellPixels = cell2mat( ptCellPixels );    
        cellStats(i).ptCellPixels = ptCellPixels;
        
        % get list of coordinates of pixels on cell boundary
%         imCurCellMask = (imLabelCellSeg == i);
%         imCurCellBoundaryMask = imdilate(imCurCellMask, streldisknd(3 * ones(1,ndims(imInput)))) - imCurCellMask;
%         cellStats(i).BoundaryPixelIdxList = find( imCurCellBoundaryMask ); % (can be faster)

        cellStats(i).BoundaryPixelIdxList = find( bwperim(imLabelCellSeg == i) ); % (can be faster)

        ptCellBoundaryPixels = cell(1, ndims(imInput));
        [ptCellBoundaryPixels{:}] = ind2sub( size(imInput), cellStats(i).BoundaryPixelIdxList );
        ptCellBoundaryPixels = cell2mat( ptCellBoundaryPixels );    
        
        cellStats(i).ptCellBoundaryPixels = ptCellBoundaryPixels;
        
        cellStats(i).Perimeter = numel( cellStats(i).BoundaryPixelIdxList );
        
        % get list of coordinates of pixels on cell convex hull
        [K, v] = convhulln( ptCellBoundaryPixels(:, [2, 1, 3:ndims(imInput)]) );
        
        cellStats(i).ConvexArea = v;          
        cellStats(i).ptCellConvexHull = ptCellBoundaryPixels(K, :);
        
        % contrast
        meanIntensityInterior = mean( imAdjusted( cellStats(i).PixelIdxList ) );
        meanIntensityExterior = mean( imAdjusted( cellStats(i).BoundaryPixelIdxList ) );
        cellStats(i).Contrast = std( imAdjusted(cellStats(i).BoundaryPixelIdxList)  / meanIntensityInterior );
        
        % ellipticity
        cellStats(i).Ellipticity = ComputeEllipticVariance( cellStats(i).ptCellPixels, spacing );
        
    end
    
    cellContrast = [ cellStats.Contrast ];
    imCellContrast = imLabelCellSeg;
    imCellContrast( imLabelCellSeg > 0 ) = cellContrast( imLabelCellSeg( imLabelCellSeg > 0 ) );
    imseriesmaskshowrgb( imCellContrast, segMaskRGB );
    set(gcf, 'Name', 'Cell Contrast' );

    cellEllipticity = [ cellStats.Ellipticity ];
    imCellEllipticity = imLabelCellSeg;
    imCellEllipticity( imLabelCellSeg > 0 ) = cellEllipticity( imLabelCellSeg( imLabelCellSeg > 0 ) );
    imseriesmaskshowrgb( imCellEllipticity, segMaskRGB );
    set(gcf, 'Name', 'Cell Ellipticity' );
    
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
    
    imseriesmaskshowrgb( imGradMag, {segMaskRGB, objSeedMaskRGB} )
    hold on;
    plot( displayAdjacencyGraphX, displayAdjacencyGraphY, 'go-', 'LineWidth', 2.0 );
    hold off;
    
    % perform merging
    mergeCostMatrixConvexity = adjMatrix;    
    mergeCostMatrixWatershedSaliency = adjMatrix;
    mergeCostMatrixAppearance = adjMatrix;
    mergeCostMatrixAppearanceHistogram = adjMatrix; 
    
        % initialize merge cost matrix
        [v1, v2, val] = find( adjMatrix );
        
        imLocalRange = rangefilt( imAdjusted );
        imLocalStd = stdfilt( imAdjusted );
        imLocalEntropy = entropyfilt( mat2gray(imAdjusted) );
        
        numHistogramBins = 16;
        binSize = range(imAdjusted(:)) / numHistogramBins;        
        imQuantized = floor( imAdjusted / binSize ) + 1;
        imQuantized( imQuantized > numHistogramBins ) = numHistogramBins;
        
        for k = 1:numel(v1)
            
            i = v1(k);
            j = v2(k);
            
            if i >= j 
                continue;
            end

            % convexity 
            c1 = cellStats(i).Area / cellStats(i).ConvexArea;            

            c2 = cellStats(j).Area / cellStats(j).ConvexArea;            

            ptMergedCellHull = [ cellStats(i).ptCellConvexHull; cellStats(j).ptCellConvexHull ];
            ptMergedCellHull = ptMergedCellHull(:, [2, 1, 3:ndims(imInput)]);
            [K, v] = convhulln( ptMergedCellHull );
            c12 = (cellStats(i).Area + cellStats(j).Area) / v;            
            
            curConvexityFactor = eps + c12 / min([c1, c2]);
            
            mergeCostMatrixConvexity(i, j) = curConvexityFactor; % higher the better
            mergeCostMatrixConvexity(j, i) = curConvexityFactor;

%             % contrast between cell interior and boundary
%             meanIntensityInterior1 = mean( imAdjusted( cellStats(i).PixelIdxList ) );
%             meanIntensityExterior1 = mean( imAdjusted( cellStats(i).BoundaryPixelIdxList ) );
%             %contrast1 = meanIntensityInterior1 / meanIntensityExterior1;
%             contrast1 = 1 / std( imAdjusted(cellStats(i).BoundaryPixelIdxList) / meanIntensityInterior1 );
%             
%             meanIntensityInterior2 = mean( imAdjusted( cellStats(j).PixelIdxList ) );
%             meanIntensityExterior2 = mean( imAdjusted( cellStats(j).BoundaryPixelIdxList ) );
%             %contrast2 = meanIntensityInterior2 / meanIntensityExterior2;                
%             contrast2 = 1 / std( imAdjusted(cellStats(j).BoundaryPixelIdxList)  / meanIntensityInterior2 );
%             
%             imMergedCellMask = imclose(ismember( imLabelCellSeg, [i, j] ), neighStrel);            
%             imMergedCellBoundary = bwperim( imMergedCellMask );
%             meanIntensityInterior12 = mean( imAdjusted( imMergedCellMask ) );
%             meanIntensityExterior12 = mean( imAdjusted( imMergedCellBoundary ) );
%             %contrast12 = meanIntensityInterior12 / meanIntensityExterior12;
%             contrast12 = 1 / std( imAdjusted( imMergedCellBoundary ) / meanIntensityInterior12 );
%             
%             mergeCostMatrixContrast(i, j) = contrast12 / min([contrast1, contrast2]);
%             mergeCostMatrixContrast(j, i) = contrast12 / min([contrast1, contrast2]);

            % contrast between cell interior and boundary
            imComponentCellMaskUnion = ismember( imLabelCellSeg, [i, j] );
            imMergedCellMask = imclose(imComponentCellMaskUnion, neighStrel);            
            imWatershed = imdilate(imMergedCellMask & ~imComponentCellMaskUnion, neighStrel);

            medianCellIntensity1 = median( imAdjusted( cellStats(i).PixelIdxList ) );
            medianCellIntensity2 = median( imAdjusted( cellStats(j).PixelIdxList ) );
            medianWatershedIntensity = median( imAdjusted(imWatershed) );
            
            curWatershedSaliency = eps + medianWatershedIntensity / max(medianCellIntensity1, medianCellIntensity2);
            
            mergeCostMatrixWatershedSaliency(i, j) = curWatershedSaliency; % higher the better
            mergeCostMatrixWatershedSaliency(j, i) = curWatershedSaliency;

            % appearance - texture
            featureMomentFunc = @(pixInd) ( [ mean(imAdjusted(pixInd)), ...
                                              std(imAdjusted(pixInd)), ...
                                              skewness(imAdjusted(pixInd)), ... 
                                              kurtosis(imAdjusted(pixInd)), ...
                                              mean(imLocalRange(pixInd)), ...
                                              std(imLocalRange(pixInd)), ...
                                              mean(imLocalStd(pixInd)), ... 
                                              std(imLocalStd(pixInd)) ] );
                                          
            p1 = featureMomentFunc( cellStats(i).PixelIdxList );
            p2 = featureMomentFunc( cellStats(j).PixelIdxList );
            
            p1 = p1 / norm(p1);
            p2 = p1 / norm(p2);
            
            appSimilarity = eps + dot(p1, p2);
            
            mergeCostMatrixAppearance(i, j) = appSimilarity;
            mergeCostMatrixAppearance(j, i) = appSimilarity;
            
            % appearance -- histogram            
            hist1 = tabulate( imQuantized( cellStats(i).PixelIdxList ) );
            hist1 = hist1(:,3) / 100;
            hist1( numel(hist1)+1 : numHistogramBins ) = 0;            
            
            hist2 = tabulate( imQuantized( cellStats(j).PixelIdxList ) );
            hist2 = hist2(:,3) / 100;
            hist2( numel(hist2)+1 : numHistogramBins ) = 0;
            
            bhatcoef = eps + dot( sqrt( hist1 ), sqrt( hist2 ) );
            
            mergeCostMatrixAppearanceHistogram(i, j) = bhatcoef;
            mergeCostMatrixAppearanceHistogram(j, i) = bhatcoef;
            
        end

        % visualize convexity metric
        imseriesmaskshowrgb( imLabelCellSeg, {segMaskRGB, objSeedMaskRGB} );
        set(gcf, 'Name', 'Convexity Factor' );
        colorbar;
        
        jetColorMap = colormap( 'jet' );
        numColors = size(jetColorMap, 1 );

        [v1, v2, mergeCostConvexity] = find( mergeCostMatrixConvexity );        
        mergeCostConvexityNmzd = mat2gray(mergeCostConvexity);

        for k = 1:numel(v1)
            
            i = v1(k);
            j = v2(k);
            
            if i >= j 
                continue;
            end

            lineColor = jetColorMap( 1 + round(mergeCostConvexityNmzd(k) * (numColors-1)), : );
            
            hold on;
            plot( [cellStats(i).ptCellSeedPoint(2), cellStats(j).ptCellSeedPoint(2)], ...
                  [cellStats(i).ptCellSeedPoint(1), cellStats(j).ptCellSeedPoint(1)], ...
                  '-', 'Color', lineColor, 'LineWidth', 1.0 + 10.0 * mergeCostConvexityNmzd(k) );
            hold off;
            
        end

        % visualize appearance similarity metric
        imseriesmaskshowrgb( imLabelCellSeg, {segMaskRGB, objSeedMaskRGB} );
        set(gcf, 'Name', 'Appearance Similarity' );
        colorbar;
        
        [v1, v2, mergeCostAppearance] = find( mergeCostMatrixAppearance );        
        mergeCostAppearanceNmzd = mat2gray(mergeCostAppearance);

        for k = 1:numel(v1)
            
            i = v1(k);
            j = v2(k);
            
            if i >= j 
                continue;
            end

            lineColor = jetColorMap( 1 + round(mergeCostAppearanceNmzd(k) * (numColors-1)), : );
            
            hold on;
            plot( [cellStats(i).ptCellSeedPoint(2), cellStats(j).ptCellSeedPoint(2)], ...
                  [cellStats(i).ptCellSeedPoint(1), cellStats(j).ptCellSeedPoint(1)], ...
                  '-', 'Color', lineColor, 'LineWidth', 1.0 + 10.0 * mergeCostAppearanceNmzd(k) );
            hold off;
            
        end
        
        % visualize appearance similarity -- histogram
        imseriesmaskshowrgb( imLabelCellSeg, {segMaskRGB, objSeedMaskRGB} );
        set(gcf, 'Name', 'Appearance Similarity Histogram' );
        colorbar;
        
        [v1, v2, mergeCostAppearanceHistogram] = find( mergeCostMatrixAppearanceHistogram );        
        mergeCostAppearanceHistogramNmzd = mat2gray(mergeCostAppearanceHistogram);

        for k = 1:numel(v1)
            
            i = v1(k);
            j = v2(k);
            
            if i >= j 
                continue;
            end

            lineColor = jetColorMap( 1 + round(mergeCostAppearanceHistogram(k) * (numColors-1)), : );
            
            hold on;
            plot( [cellStats(i).ptCellSeedPoint(2), cellStats(j).ptCellSeedPoint(2)], ...
                  [cellStats(i).ptCellSeedPoint(1), cellStats(j).ptCellSeedPoint(1)], ...
                  '-', 'Color', lineColor, 'LineWidth', 1.0 + 10.0 * mergeCostAppearanceHistogram(k) );
            hold off;
            
        end
        
        % visualize watershed saliency metric
        imseriesmaskshowrgb( imLabelCellSeg, {segMaskRGB, objSeedMaskRGB} );
        set(gcf, 'Name', 'Watershed Saliency' );
        colorbar;

        [v1, v2, mergeCostWatershedSaliency] = find( mergeCostMatrixWatershedSaliency );        
        mergeCostWatershedSaliencyNmzd = mat2gray(mergeCostWatershedSaliency);

        for k = 1:numel(v1)
            
            i = v1(k);
            j = v2(k);
            
            if i >= j 
                continue;
            end

            lineColor = jetColorMap( 1 + round(mergeCostWatershedSaliencyNmzd(k) * (numColors-1)), : );
            
            hold on;
            plot( [cellStats(i).ptCellSeedPoint(2), cellStats(j).ptCellSeedPoint(2)], ...
                  [cellStats(i).ptCellSeedPoint(1), cellStats(j).ptCellSeedPoint(1)], ...
                  '-', 'Color', lineColor, 'LineWidth', 1.0 + 10.0 * mergeCostWatershedSaliencyNmzd(k) );
            hold off;
            
        end
        
        % visualize combined merge likelihood
        imseriesmaskshowrgb( imLabelCellSeg, {segMaskRGB, objSeedMaskRGB} );
        set(gcf, 'Name', 'Combined Merge Likelihood' );
        colorbar;

        [v1, v2, ~] = find( adjMatrix );        
        %mergeCostCombinedNmzd = mat2gray( mergeCostWatershedSaliencyNmzd + mergeCostConvexityNmzd + 0.25 * mergeCostAppearanceHistogramNmzd );
        %mergeCostCombinedNmzd = mat2gray( mergeCostWatershedSaliencyNmzd + 0.25 * mergeCostAppearanceHistogramNmzd );
        %mergeCostCombinedNmzd = mat2gray( mergeCostWatershedSaliencyNmzd .* mergeCostConvexityNmzd .* mergeCostAppearanceHistogramNmzd );
        mergeCostCombinedNmzd = mat2gray( min( [mergeCostWatershedSaliencyNmzd, mergeCostAppearanceHistogramNmzd], [], 2) );
        
        mergeCostMatrixCombined = sparse(v1, v2, mergeCostCombinedNmzd, numCells, numCells);
        
        for k = 1:numel(v1)
            
            i = v1(k);
            j = v2(k);
            
            if i >= j 
                continue;
            end

            lineColor = jetColorMap( 1 + round(mergeCostCombinedNmzd(k) * (numColors-1)), : );
            
            hold on;
            plot( [cellStats(i).ptCellSeedPoint(2), cellStats(j).ptCellSeedPoint(2)], ...
                  [cellStats(i).ptCellSeedPoint(1), cellStats(j).ptCellSeedPoint(1)], ...
                  '-', 'Color', lineColor, 'LineWidth', 1.0 + 10.0 * mergeCostCombinedNmzd(k) );
            hold off;
            
        end        
        
        
        % visualize maximal similarity edges
        imseriesmaskshowrgb( imLabelCellSeg, {segMaskRGB, objSeedMaskRGB} );
        set(gcf, 'Name', 'Maximal Similarity Edges' );
        colorbar;
        
        [v1, v2, mergeCostCombinedNmzd] = find( mergeCostMatrixCombined );        
        
        for k = 1:numel(v1)
            
            i = v1(k);
            j = v2(k);
            
            if i >= j 
                continue;
            end
            
            % check if edge is locally maximal
            if ~IsEdgeLocallyMaximal(mergeCostMatrixCombined, i, j)           
                
                hold on;
                plot( [cellStats(i).ptCellSeedPoint(2), cellStats(j).ptCellSeedPoint(2)], ...
                      [cellStats(i).ptCellSeedPoint(1), cellStats(j).ptCellSeedPoint(1)], ...
                      '-', 'Color', 'w', 'LineWidth', 1.0 + 10.0 * mergeCostCombinedNmzd(k) );
                hold off;
           
            else

                if mergeCostCombinedNmzd(k) > 0.5 
                    lineColor = jetColorMap( 1 + round(mergeCostCombinedNmzd(k) * (numColors-1)), : );
                else
                    lineColor = 1 - mergeCostCombinedNmzd(k) * ones(1,3);
                end
                
                
                hold on;
                plot( [cellStats(i).ptCellSeedPoint(2), cellStats(j).ptCellSeedPoint(2)], ...
                      [cellStats(i).ptCellSeedPoint(1), cellStats(j).ptCellSeedPoint(1)], ...
                      '-', 'Color', lineColor, 'LineWidth', 1.0 + 10.0 * mergeCostCombinedNmzd(k) );
                hold off;
                
                
            end
            
        end

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
    
    % region bounding box    
    regProps.BoundingBox = ([ min(ptPixelLocations) ; max(ptPixelLocations) ])';
    
    % construct indices to crop region bounding box from any image
    bboxPadding = 4;

    indCropBBox = cell(1,ndims(imInput));
    for j = 1:ndims(imInput)
        indStart = max(1, regProps.BoundingBox(j,1) - bboxPadding);
        indEnd = min( size(imInput,j), regProps.BoundingBox(j,2) + bboxPadding );
        indCropBBox{j} = indStart:indEnd;
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
    imComponentCellMaskUnion = ismember( imLabelCellSeg, [i, j] );
    imMergedCellMask = imclose(imComponentCellMaskUnion, neighStrel);            
    imWatershed = imdilate(imMergedCellMask & ~imComponentCellMaskUnion, neighStrel);

    medianCellIntensity1 = median( imAdjusted( regProps(i).PixelIdxList ) );
    medianCellIntensity2 = median( imAdjusted( regProps(j).PixelIdxList ) );
    medianWatershedIntensity = median( imAdjusted(imWatershed) );

    simWatershed = eps + medianWatershedIntensity / max(medianCellIntensity1, medianCellIntensity2);

    % shape/convexity 
    c1 = regProps(i).Area / regProps(i).ConvexArea;
    c2 = regProps(j).Area / regProps(j).ConvexArea;            

    ptMergedCellHull = [ regProps(i).ptConvexHull; regProps(j).ptConvexHull ];
    ptMergedCellHull = ptMergedCellHull(:, [2, 1, 3:ndims(imAdjusted)]);
    [~, v] = convhulln( ptMergedCellHull );
    c12 = (regProps(i).Area + regProps(j).Area) / v;            

    simShape = eps + c12 / min([c1, c2]);
   
end

function DisplayWeightedRegionAdjacencyGraph( hAxis, simMatrix, regProps, flagShowMaximallySimilarEdgesOnly, simThresh )

    if ~exist( 'flagShowMaximallySimilarEdgesOnly', 'var' )
        flagShowMaximallySimilarEdgesOnly = false;
    end
    
    gca(hAxis);
    colorbar;

    [v1, v2, simVal] = find( simMatrix );        
    simValNmzd = mat2gray(simVal);

    jetColorMap = colormap( 'jet' );
    
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
                  '-', 'Color', 'w', 'LineWidth', 1.0 + 10.0 * simValNmzd(k) );
            hold off;

        else
        
            
            if exist( 'simThresh', 'var' ) && simValNmzd > simThresh
                lineColor = jetColorMap( 1 + round(simValNmzd(k) * (numColors-1)), : );
            else
                lineColor = 1 - mergeCostCombinedNmzd(k) * ones(1,3);
            end

            hold on;
            plot( [regProps(i).ptCellSeedPoint(2), regProps(j).ptCellSeedPoint(2)], ...
                  [regProps(i).ptCellSeedPoint(1), regProps(j).ptCellSeedPoint(1)], ...
                  '-', 'Color', lineColor, 'LineWidth', 1.0 + 10.0 * simValNmzd(k) );
            hold off;
            
        end

    end

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

    