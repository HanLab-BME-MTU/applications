function testLearnOversegMergeMetric

clc
clear 
close all

%***********************************************************************
% 
%                           PARAMETERS 
%
%***********************************************************************

    dataFileDir = 'C:\deepak\data\regionMergingTraining2D';
    
    % critical parameters (needs tuning)
    cellDiameterRange = [6, 20];
    numPathQuantLevels = 64;
    numHistQuantLevels = 64;
    
    % less-critical parameters (leave defaults)
    threshNeighOverlap = 5;
    neighSearchRad = 2;    
    bndContStrelRad = 5;
    
%***********************************************************************

% load data
if ~exist( 'dataFilePath', 'var' )
    [fileName,pathName] = uigetfile(fullfile(dataFileDir, '*.mat; *.rann'), 'Select the data file' );   
    dataFilePath = fullfile( pathName, fileName )
    [~,~,fileExt] = fileparts(fileName);
end

% for parallelization
if matlabpool( 'size' ) == 0
    matlabpool open;
end            

switch fileExt
    
    case '.mat'
        
        dataFileContents = load( dataFilePath );
        imInput = dataFileContents.im;
        spacing = [0.5, 0.5];

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


        imAdjusted( ~imThresh ) = 0;

        % compute cell seed points
        imDistMap = callMatitkFilter( 'FSignedMaurerDistancemap', {} , imThresh, [], [], spacing  );
        imDistMap( imDistMap > 0 ) = 0;
        imDistMap = -1 * imDistMap;

        numLoGScales = min(20, cellDiameterRange(2)-cellDiameterRange(1)+1);
        [imCellSeedPoints, imResponse] = detectBlobsUsingAdaptiveMultiscaleLoG( imAdjusted, ...
                                                                                imDistMap, ... 
                                                                                'blobDiameterRange', cellDiameterRange, ...
                                                                                'numLoGScales', numLoGScales, ...
                                                                                'spacing', spacing, ...
                                                                                'debugMode', false );

        imCellSeedPoints(~imThresh) = 0; % remove seed points in background
        meanResponseBgnd = mean( imResponse(~imThresh) );
        stdResponseBgnd = std( imResponse(~imThresh) );
        imCellSeedPoints( imCellSeedPoints < (meanResponseBgnd + 2.0 * stdResponseBgnd) ) = 0; % remove seed points with weak LoG response

        objSeedMask = double(imdilate(imCellSeedPoints > 0, streldisknd([3,3])));
        objSeedMaskRGB = cat(3, objSeedMask, zeros([size(objSeedMask), 2]));

        % segment nuclei
        imBackgroundSeeds = bwmorph( ~imThresh, 'shrink', 1 );

        imFeature = 1 - mat2gray( imResponse );
        imMinimaImposedAtSeeds = imimposemin( imFeature, imCellSeedPoints | imBackgroundSeeds );
        imMinimaImposedAtSeedsDisplay = imMinimaImposedAtSeeds;
        imMinimaImposedAtSeedsDisplay( imMinimaImposedAtSeeds == -Inf ) = min( imMinimaImposedAtSeeds( imMinimaImposedAtSeeds ~= -Inf ) ) - 3 * eps;

        imseriesmaskshow( imMinimaImposedAtSeedsDisplay, {objSeedMask, imBackgroundSeeds}, 'spacing', spacing );    
        set( gcf, 'Name', 'Minima Imposed at Seed Points' );

        L = watershed( imMinimaImposedAtSeeds );
        objLabelInd = unique( L( imCellSeedPoints > 0 ) );
        L( ~ismember(L, objLabelInd) ) = 0;
        L(~imThresh) = 0;

        % Post-processing - remove very small objects
        minCellVolume = pi * (0.5 * 0.5 * min(cellDiameterRange))^2;
        regstats = regionprops(L, 'Area');
        objAreaList = [regstats.Area] * prod(spacing);
        smallObjInd = find( objAreaList < minCellVolume );
        L( ismember( L, smallObjInd ) ) = 0;
        imLabelCellSeg = bwlabeln( L > 0 );

        % segmentation boundary mask
        [segMaskRGB, segColorMap] = label2rgbND( imLabelCellSeg );

        [gx, gy] = gradient(imLabelCellSeg); 
        imSegBoundaryMask = (gx ~= 0 | gy ~= 0);

        imseriesmaskshow(imInput, {imSegBoundaryMask, objSeedMask});
        set(gcf, 'Name', 'cell segmentation boundary mask overlay');

        % cell stats
        numCells = max(imLabelCellSeg(:));    
        cellStats = [];    
        for cellId = 1:numCells                       

            regProps = ComputeRegionProperties(imLabelCellSeg, cellId);

            % seed point location
            indCurRegionSeedPoint = regProps.PixelIdxList( imCellSeedPoints(regProps.PixelIdxList) > 0 );
            if isempty( indCurRegionSeedPoint )            
                regProps.ptCellSeedPoint = regProps.Centroid;
            else

                [maxval, maxind] = max( imCellSeedPoints( indCurRegionSeedPoint ) );
                indCurRegionSeedPoint = indCurRegionSeedPoint( maxind );    

                regProps.ptCellSeedPoint = ind2submat( size(imCellSeedPoints), indCurRegionSeedPoint );
                regProps.ptIndCellSeedPoint = indCurRegionSeedPoint;

            end

            cellStats = [cellStats; regProps];

        end

        % build region adjacency graph (RAG)
        adjMatrix = sparse(numCells, numCells);

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

        imseriesmaskshow(imInput, {imSegBoundaryMask, objSeedMask});
        set(gcf, 'Name', 'cell neighborhood graph overlayed on input image');
        hold on;
        plot( displayAdjacencyGraphX, displayAdjacencyGraphY, 'go-', 'LineWidth', 2.0 );
        hold off;

        imseriesmaskshow(imLabelCellSeg, {imSegBoundaryMask, objSeedMask});
        set(gcf, 'Name', 'cell neighborhood graph overlayed on label map');
        hold on;
        plot( displayAdjacencyGraphX, displayAdjacencyGraphY, 'go-', 'LineWidth', 2.0 );
        hold off;

        % compute path profiles for all edges
        costMatrix = (1 - mat2gray(imDistMap)) + (1 - mat2gray(imAdjusted));

        [v1, v2, ~] = find( adjMatrix );        

        edgePathList = cell(numel(v1),1);
        
        pathClock = tic;
        parfor k = 1:numel(v1)

            i = v1(k);
            j = v2(k);

            if i >= j 
                continue;
            end

            srcNode = submat2ind( size(imInput), cellStats(i).ptCentroid );
            destNode = submat2ind( size(imInput), cellStats(j).ptCentroid );

            imRegionUnionMask = ismember( imLabelCellSeg, [i, j] );
            imRegionUnionMask = imclose(imRegionUnionMask, neighStrel);            

            [path, cost] = findMinimalPathOnCostMatrix(costMatrix, srcNode, destNode, ...
                                                       'roiMask', imRegionUnionMask, 'debugMode', false, ...
                                                       'costWeight', 0.5);

            edgePathList{k} = path;    
        end
        timeElapsed = toc(pathClock)

        edgePath = cell(numCells, numCells);
        pathMask = zeros(size(imInput));
        for k = 1:numel(v1)

            i = v1(k);
            j = v2(k);

            if i >= j 
                continue;
            end

            path = edgePathList{k};
            edgePath{i,j} = path;
            edgePath{j,i} = path;
            pathMask(path) = 1;

        end

        diskstrel = streldisknd(ones(1,ndims(imAdjusted)));
        pathMaskDisp = imdilate(pathMask, diskstrel);
        imseriesmaskshow( imAdjusted, {imSegBoundaryMask, objSeedMask, pathMaskDisp} );                           

        % ask user to point to edges which have to be merged
        imseriesmaskshow( imAdjusted, {imSegBoundaryMask, objSeedMask} );
        set(gcf, 'Name', 'Annotate Regions to be Merged');
        hold on;
        plot( displayAdjacencyGraphX, displayAdjacencyGraphY, 'go-', 'LineWidth', 2.0 );
        hold off;

        [x, y] = getpts(gcf);
        ptSeed = [x, y];
        ptSeedInd = sub2ind(size(imAdjusted), round(ptSeed(:,2)), round(ptSeed(:,1))); 
        c1 = imLabelCellSeg( ptSeedInd(1:2:numel(ptSeedInd)) );
        c2 = imLabelCellSeg( ptSeedInd(2:2:numel(ptSeedInd)) );
        adjMatrixMerge = 2 * sparse([c1; c2], [c2; c1], ones(numel(c1) * 2, 1), numCells, numCells) - spones(adjMatrix);

        imseriesmaskshow(imInput, {imSegBoundaryMask, objSeedMask});
        set(gcf, 'Name', 'User Selected Merge Edges');
        hold on;
        plot( displayAdjacencyGraphX, displayAdjacencyGraphY, 'go-', 'LineWidth', 2.0 );
        hold off;

        for k = 1:numel(c1)

            i = c1(k);
            j = c2(k);

            hold on;
            plot( [cellStats(i).ptCellSeedPoint(2), cellStats(j).ptCellSeedPoint(2)], ...
                  [cellStats(i).ptCellSeedPoint(1), cellStats(j).ptCellSeedPoint(1)], ...
                  '-', 'Color', 'r', 'LineWidth', 4.0 );
            hold off;

        end

        pathMaskMerge = zeros(size(imAdjusted));
        pathMaskDontMerge = zeros(size(imAdjusted));

        for k = 1:numel(v1)

            i = v1(k);
            j = v2(k);

            if i >= j 
                continue;
            end

            path = edgePath{i,j};
            if adjMatrixMerge(i,j) > 0
                pathMaskMerge(path) = path;
            else
                pathMaskDontMerge(path) = path;
            end    

        end

        diskstrel = streldisknd(ones(1,ndims(imAdjusted)));
        imseriesmaskshow( imAdjusted, {imSegBoundaryMask, objSeedMask, ...
                                       imdilate(pathMaskMerge, diskstrel), ...
                                       imdilate(pathMaskDontMerge, diskstrel)} );                           
        
        
    case '.rann'        
        
        % load data from annotation file
        curAnnotationData = load(dataFilePath, '-mat', 'annotationData' );    
        curAnnotationData = curAnnotationData.annotationData;

        imInput = curAnnotationData.imInput;
        imThresh = curAnnotationData.imThresh;
        imCellSeedPoints = curAnnotationData.imCellSeedPoints;
        imLabelCellSeg = curAnnotationData.imLabelCellSeg;
        numCells = curAnnotationData.numCells;
        adjMatrix = curAnnotationData.adjMatrix;
        adjMatrixMerge = curAnnotationData.adjMatrixMerge;
        spacing = [0.5, 0.5];
        
        neighStrel = streldisknd(neighSearchRad*ones(1, ndims(imInput)));
        
        % standardize the intensity range of the image 
        imageDynamicRange = ComputeImageDynamicRange( imInput, 99.9 );
        imAdjusted = mat2gray(imInput, imageDynamicRange);    

        % apply median filter to remove any spotty noise
        imAdjusted = medfilt2( imAdjusted, [3,3] );    

        % compute cell seed points
        imDistMap = callMatitkFilter( 'FSignedMaurerDistancemap', {} , imThresh, [], [], spacing  );
        imDistMap( imDistMap > 0 ) = 0;
        imDistMap = -1 * imDistMap;
        
        % compute basic properties of cell region
        cellStats = [];    
        for cellId = 1:numCells                       

            regProps = ComputeRegionProperties(imLabelCellSeg, cellId);

            % add seed point location
            indCurRegionSeedPoint = regProps.PixelIdxList( imCellSeedPoints(regProps.PixelIdxList) > 0 );
            if isempty( indCurRegionSeedPoint )            
                regProps.ptCellSeedPoint = regProps.Centroid;
            else
                [maxval, maxind] = max( imCellSeedPoints( indCurRegionSeedPoint ) );
                indCurRegionSeedPoint = indCurRegionSeedPoint( maxind );    

                regProps.ptCellSeedPoint = ind2submat( size(imCellSeedPoints), indCurRegionSeedPoint );
            end

            cellStats = [cellStats; regProps];

        end        
        
        % display stuff
        [segMaskRGB, segColorMap] = label2rgbND( imLabelCellSeg );
        objSeedMask = double(imdilate(imCellSeedPoints > 0, streldisknd([3,3])));
        
        [gx, gy] = gradient(imLabelCellSeg); 
        imSegBoundaryMask = (gx ~= 0 | gy ~= 0);

        [v1, v2, ~] = find( adjMatrix );        
        displayAdjacencyGraphX = [];
        displayAdjacencyGraphY = [];    
        
        for k = 1:numel(v1)

            i = v1(k);
            j = v2(k);

            if i >= j 
                continue;
            end

            displayAdjacencyGraphX(:, end+1 ) = [cellStats(i).ptCellSeedPoint(2); cellStats(j).ptCellSeedPoint(2)];
            displayAdjacencyGraphY(:, end+1 ) = [cellStats(i).ptCellSeedPoint(1); cellStats(j).ptCellSeedPoint(1)];
            
        end        
        
        imseriesmaskshow(imInput, {imSegBoundaryMask, objSeedMask});
        set(gcf, 'Name', 'cell segmentation boundary mask overlay');        
        hold on;
        plot( displayAdjacencyGraphX, displayAdjacencyGraphY, 'go-', 'LineWidth', 2.0 );
        hold off;
        
        imseriesmaskshow(imLabelCellSeg, {imSegBoundaryMask, objSeedMask});
        set(gcf, 'Name', 'cell segmentation label map');
        hold on;
        plot( displayAdjacencyGraphX, displayAdjacencyGraphY, 'go-', 'LineWidth', 2.0 );
        hold off;

        % compute path profiles for all edges
        costMatrix = (1 - mat2gray(imDistMap)) + (1 - mat2gray(imAdjusted));

        [v1, v2, ~] = find( adjMatrix );        
        edgePathList = cell(numel(v1),1);
        
        pathClock = tic;
        parfor k = 1:numel(v1)

            i = v1(k);
            j = v2(k);

            if i >= j 
                continue;
            end

            srcNode = submat2ind( size(imInput), cellStats(i).ptCentroid );
            destNode = submat2ind( size(imInput), cellStats(j).ptCentroid );

            imRegionUnionMask = ismember( imLabelCellSeg, [i, j] );
            imRegionUnionMask = imclose(imRegionUnionMask, neighStrel);            

            [path, cost] = findMinimalPathOnCostMatrix(costMatrix, srcNode, destNode, ...
                                                       'roiMask', imRegionUnionMask, 'debugMode', false, ...
                                                       'costWeight', 0.5);

            edgePathList{k} = path;    
        end
        timeElapsed = toc(pathClock)

        edgePath = cell(numCells, numCells);
        pathMask = zeros(size(imInput));
        for k = 1:numel(v1)

            i = v1(k);
            j = v2(k);

            if i >= j 
                continue;
            end

            path = edgePathList{k};
            edgePath{i,j} = path;
            edgePath{j,i} = path;
            pathMask(path) = 1;

        end

        diskstrel = streldisknd(ones(1,ndims(imAdjusted)));
        pathMaskDisp = imdilate(pathMask, diskstrel);
        imseriesmaskshow( imAdjusted, {imSegBoundaryMask, objSeedMask, pathMaskDisp} );  
        
end

% compute features to help a classifier decide whether or not to merge the regions            
[v1, v2,~] = find( adjMatrix ); 

uniModalFitErrorMat = spones(adjMatrix);
watershedSaliencyMat = spones(adjMatrix);
brightnessSimilarityMat = spones(adjMatrix);
boundaryContinuityMat = spones(adjMatrix);
overallBrightnessSimilarityMat = spones(adjMatrix);

imAdjustedQuant = round( mat2gray(imAdjusted) * (numPathQuantLevels-1) );

bndContStrel = streldisknd( bndContStrelRad * ones(1,ndims(imInput)) ); 
bndContStrel = bndContStrel / sum(bndContStrel(:));

hProfile = figure;
hFeature = figure;

for k = 1:numel(v1)
    
    i = v1(k);
    j = v2(k);

    if i >= j 
        continue;
    end
    
    curPath = edgePath{i,j};
    curPathMask = zeros(size(imInput));
    curPathMask(curPath) = 1;
    curPathMask = imdilate(curPathMask, streldisknd(2 * ones(1,ndims(imAdjusted))));    

    % compute unimodality feature
    pathIntensities = imAdjustedQuant( curPathMask > 0 );
    pathPixelLabels = imLabelCellSeg( curPathMask > 0 );
    
    pixval_Ri = pathIntensities(pathPixelLabels == i);
    med_Ri = median( pixval_Ri ); % robust center estimate
    mad_Ri = median( abs(pixval_Ri - med_Ri) ); % robust scale estimate
    sigma_Ri = 1.4826 * mad_Ri;
    uniModalFitError_Ri = numel( find( abs(pathIntensities - med_Ri) > 2.5 * sigma_Ri ) ) / numel(pathIntensities);
    
    pixval_Rj = pathIntensities(pathPixelLabels == j);
    med_Rj = median( pixval_Rj );
    mad_Rj = median( abs(pixval_Rj - med_Rj) ) + eps;
    sigma_Rj = 1.4826 * mad_Rj;
    uniModalFitError_Rj = numel( find( abs(pathIntensities - med_Rj) > 2.5 * sigma_Rj ) ) / numel(pathIntensities);
    
    %combUnimodelFitError = eps + ((numel(pixval_Ri) * uniModalFitError_Ri) + (numel(pixval_Rj) * uniModalFitError_Rj)) / numel(pathPixelLabels);
    combUnimodelFitError = eps + max( [uniModalFitError_Ri, uniModalFitError_Rj] );
    uniModalFitErrorMat(i,j) = combUnimodelFitError;
    uniModalFitErrorMat(j,i) = combUnimodelFitError;
    
    % compute watershed saliency measure
    ptBBoxCorners = [ cellStats(i).BoundingBox, cellStats(j).BoundingBox ];
    cropind = { min(ptBBoxCorners(1,:)):max(ptBBoxCorners(1,:)), min(ptBBoxCorners(2,:)):max(ptBBoxCorners(2,:)) };
    imLabelCellSegCropped = imLabelCellSeg( cropind{:} );
    imInputCropped = imAdjusted( cropind{:} );
    
    mergeStrel = ones((2*neighSearchRad+1)*ones(1, ndims(imInput)));
    imComponentMask1 = imclose(imLabelCellSegCropped == i, mergeStrel);
    imComponentMask2 = imclose(imLabelCellSegCropped == j, mergeStrel);
    imComponentCellMaskUnion = imfill( imComponentMask1 | imComponentMask2, 'holes' );
    imMergedCellMask = imclose(imComponentCellMaskUnion, mergeStrel);            
    L = bwlabeln(imMergedCellMask & ~imComponentCellMaskUnion); 
    stats = regionprops( L, 'Area' );
    [maxval, maxind] = max( [stats.Area] );
    imWatershed = imdilate(L == maxind, neighStrel) & imMergedCellMask;
    
    medianCellIntensity1 = median( imInputCropped( imLabelCellSegCropped == i ) );
    medianCellIntensity2 = median( imInputCropped( imLabelCellSegCropped == j ) );
    medianCellIntensity = median( imInputCropped( imComponentCellMaskUnion > 0 ) );
    medianWatershedIntensity = median( imInputCropped(imWatershed > 0) );

    watershedSaliency = eps + max(medianCellIntensity1, medianCellIntensity2) / medianWatershedIntensity;
    
    watershedSaliencyMat(i,j) = watershedSaliency;
    watershedSaliencyMat(j,i) = watershedSaliency;

    % boundary continuation
    imLocalCurvature = imfilter( double(imMergedCellMask), bndContStrel);
    imLocalCurvature(imLocalCurvature > 0.99) = 0;    
    imJunction = bwperim(imMergedCellMask) & (L == maxind);
    bndContinuity = eps + max(abs(imLocalCurvature(imJunction > 0) - 0.5));
    
    boundaryContinuityMat(i,j) = bndContinuity;
    boundaryContinuityMat(j,i) = bndContinuity;
    
    % brightness similarity
    hist1 = tabulate( imAdjustedQuant( cellStats(i).PixelIdxList ) );
    hist1 = hist1(:,3) / 100;
    hist1( numel(hist1)+1 : numPathQuantLevels ) = 0;            

    hist2 = tabulate( imAdjustedQuant( cellStats(j).PixelIdxList ) );
    hist2 = hist2(:,3) / 100;
    hist2( numel(hist2)+1 : numPathQuantLevels ) = 0;

    brightnessSimilarity = eps + dot( sqrt( hist1 ), sqrt( hist2 ) );
    
    brightnessSimilarityMat(i,j) = brightnessSimilarity;
    brightnessSimilarityMat(i,j) = brightnessSimilarity;
    
    % overall brightness difference
    overallBrightnessSimilarity = min(medianCellIntensity1, medianCellIntensity2) / max(medianCellIntensity1, medianCellIntensity2);
    
    overallBrightnessSimilarityMat(i,j) = overallBrightnessSimilarity;
    overallBrightnessSimilarityMat(j,i) = overallBrightnessSimilarity;
    
    % plot stuff
    figure( hFeature );    
    hold on;
    if adjMatrixMerge(i,j) > 0        
        plot(watershedSaliency, bndContinuity, 'g+' );
    else        
        plot(watershedSaliency, bndContinuity, 'r+' );
    end    
    hold off;
    xlabel( 'Watershed Saliency' );
    ylabel( 'Junction Curvature' );    
    
    figure(hProfile);
    [imMaskOverlay] = genImageMaskOverlay(imInputCropped, {imMergedCellMask, L == maxind, imJunction, curPathMask(cropind{:})}, [1 0 0; 0 1 0; 0 0 1; 1 0 1], [0.2, 0.2, 0.5, 0.5]);
    subplot(1,2,1), imshow( imMaskOverlay );
    title( { sprintf('boundary continuity = %.2f, brightness similarity = %.2f', bndContinuity, brightnessSimilarity), ...
             sprintf('watershed saliency = %.2f, unimodel fit error = %.2f', watershedSaliency, combUnimodelFitError)} );
    subplot(1,2,2);
    plot( imAdjustedQuant(curPath), 'b-o' );
    ylim([0, numPathQuantLevels]);
    
    hold on;
    
    bndInd = find( imLabelCellSeg(curPath) == 0 );
    bndInd = bndInd(1);
    plot( bndInd, imAdjustedQuant(curPath(bndInd)), 'ro' );

    plot( [1; bndInd], [med_Ri; med_Ri], 'm-' );
    plot( [1; bndInd], [med_Ri - 2.5 * sigma_Ri; med_Ri - 2.5 * sigma_Ri], 'm-' );
    plot( [1; bndInd], [med_Ri + 2.5 * sigma_Ri; med_Ri + 2.5 * sigma_Ri], 'm-' );
    
    plot( [bndInd; numel(curPath)], [med_Rj; med_Rj], 'g-' );
    plot( [bndInd; numel(curPath)], [med_Rj - 2.5 * sigma_Rj; med_Rj - 2.5 * sigma_Rj], 'g-' );
    plot( [bndInd; numel(curPath)], [med_Rj + 2.5 * sigma_Rj; med_Rj + 2.5 * sigma_Rj], 'g-' );    
    
    if adjMatrixMerge(i,j) > 0 
        title( 'Merge' );
    else
        title( 'Dont Merge' );
    end
    
    hold off;
    
end

figure, nhist( { full(uniModalFitErrorMat(adjMatrixMerge > 0)), full(uniModalFitErrorMat(adjMatrixMerge < 0))}, 'legend', { 'Merge', 'DontMerge' } );
title( 'Unimodal Fit Error' );

figure, nhist( { full(watershedSaliencyMat(adjMatrixMerge > 0)), full(watershedSaliencyMat(adjMatrixMerge < 0))}, 'legend', { 'Merge', 'DontMerge' } );
title( 'Watershed Saliency' );

figure, nhist( { full(boundaryContinuityMat(adjMatrixMerge > 0)), full(boundaryContinuityMat(adjMatrixMerge < 0))}, 'legend', { 'Merge', 'DontMerge' } );
title( 'Junction Curvature' );

figure, nhist( { full(overallBrightnessSimilarityMat(adjMatrixMerge > 0)), full(overallBrightnessSimilarityMat(adjMatrixMerge < 0))}, 'legend', { 'Merge', 'DontMerge' } );
title( 'Difference in Overall Brightness' );

imseriesmaskshow( imAdjusted, {imSegBoundaryMask, objSeedMask, pathMaskDisp} );  
DisplayWeightedRegionAdjacencyGraph( gca, uniModalFitErrorMat, cellStats );
set(gcf, 'Name', 'Path profile Unimodal Fit Error' );

imseriesmaskshow(imInput, {imSegBoundaryMask, objSeedMask});
DisplayWeightedRegionAdjacencyGraph( gca, boundaryContinuityMat, cellStats );
set(gcf, 'Name', 'Junction Curvature' );

imseriesmaskshow(imInput, {imSegBoundaryMask, objSeedMask});
DisplayWeightedRegionAdjacencyGraph( gca, overallBrightnessSimilarityMat, cellStats );
set(gcf, 'Name', 'Difference in Overall Brightness' );

end


function DisplayWeightedRegionAdjacencyGraph( hAxis, simMatrix, regProps, flagShowMaximallySimilarEdgesOnly, simThresh )

    if ~exist( 'flagShowMaximallySimilarEdgesOnly', 'var' )
        flagShowMaximallySimilarEdgesOnly = false;
    end

    if ~exist( 'simThresh', 'var' )
        simThresh = 0.0;
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
