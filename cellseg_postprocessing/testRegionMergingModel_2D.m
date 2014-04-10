function testRegionMergingModel_2D

clc
clear 
close all

%***********************************************************************
% 
%                           PARAMETERS 
%
%***********************************************************************

    dataFileDir = 'C:\deepak\data\Stefan_July2012';
    modelFileDir = 'C:\deepak\data\regionMergingTraining2D\features';
    
    %PARAMETERS.SEG_ALGO = 'RANDOM_WALKS';
    PARAMETERS.SEG_ALGO = 'WATERSHED';

    % critical parameters (needs tuning)
    PARAMETERS.cellDiameterRange = [8, 20];
    
    % less-critical parameters (recommended to leave defaults)
    PARAMETERS.threshNeighOverlap = 5; % percentage of overlap between touching regions to consider it as an edge
    PARAMETERS.neighSearchRad = 2;    
    PARAMETERS.logResponseCutoff = 1.0;
    PARAMETERS.minLocalGlobalThresholdRatio = 0.6;
    PARAMETERS.localThresholdWindowRadiusPhysp = 30;
    
%***********************************************************************

% Make sure weka is in the path
if isempty( cell2mat( strfind( javaclasspath('-dynamic'), 'weka.jar' ) ) )
    javaaddpath( fullfile('C:\Program Files\Weka-3-7', 'weka.jar') );
end

% open matlab pool for parallelization
if matlabpool( 'size' ) == 0
    matlabpool open;
end           

% load weka model file
[fileName, pathName] = uigetfile( fullfile(modelFileDir, '*.model'), 'Select Model File' );
modelFilePath = fullfile(pathName, fileName);

import weka.core.SerializationHelper;
wekaModel = weka.core.SerializationHelper.readAll( modelFilePath );

% load data
if ~exist( 'dataFilePath', 'var' )
    [fileName,pathName] = uigetfile( fullfile( dataFileDir, '*.mat' ), 'Select the data file' );   
    dataFilePath = fullfile( pathName, fileName )
end

dataFileContents = load( dataFilePath );
imInput = dataFileContents.im;
spacing = [0.5, 0.5];

%% Preprocessing
fprintf( '\nPreprocessing the input image ...\n' );

    % adjust the intensity range of the image 
    imageDynamicRange = ComputeImageDynamicRange( imInput, 99.9 );
    imAdjusted = mat2gray(imInput, imageDynamicRange);    

    % apply median filter to remove any spotty noise
    imAdjusted = medfilt2( imAdjusted, [3,3] );    

%% Cell Foreground Extraction
fprintf( '\nExtracting cell foreground ...\n' );

localWindowRadius = round(PARAMETERS.localThresholdWindowRadiusPhysp ./ spacing(1));
localWindowPace = round(localWindowRadius / 3);

imThresh = segmentCellForegroundUsingLocalMinError( imAdjusted, localWindowRadius, ...
                                                    'model', 'poisson', ...  
                                                    'localWindowPace', localWindowPace, ...
                                                    'minLocalGlobalThresholdRatio', PARAMETERS.minLocalGlobalThresholdRatio, ...
                                                    'flagParallelize', true );
                                         
imThresh = imfill( imThresh );        

imseriesmaskshow( imAdjusted, imThresh );
set( gcf, 'Name', 'Result of Thresholding' );

if strcmp(PARAMETERS.SEG_ALGO, 'WATERSHED')
    imAdjusted( ~imThresh ) = 0;
end

%% Cell Seedpoint Detection
imDistMap = callMatitkFilter( 'FSignedMaurerDistancemap', {} , imThresh, [], [], spacing  );
imDistMap( imDistMap > 0 ) = 0;
imDistMap = -1 * imDistMap;

imseriesmaskshow( imDistMap, imThresh );
set( gcf, 'Name', 'Foreground Distance Map' );    

numLoGScales = min(20, PARAMETERS.cellDiameterRange(2)-PARAMETERS.cellDiameterRange(1)+1);
[imCellSeedPoints, imResponse] = detectBlobsUsingAdaptiveMultiscaleLoG( imAdjusted, ...
                                                                        imDistMap, ... 
                                                                        'blobDiameterRange', PARAMETERS.cellDiameterRange, ...
                                                                        'numLoGScales', numLoGScales, ...
                                                                        'spacing', spacing, ...
                                                                        'debugMode', false );
                                                                    
imCellSeedPoints(~imThresh) = 0; % remove seed points in background
meanResponseBgnd = mean( imResponse(~imThresh) );
stdResponseBgnd = std( imResponse(~imThresh) );
imCellSeedPoints( imCellSeedPoints < (meanResponseBgnd + PARAMETERS.logResponseCutoff * stdResponseBgnd) ) = 0; % remove seed points with weak LoG response

objSeedMask = double(imdilate(imCellSeedPoints > 0, streldisknd([3,3])));
objSeedMaskRGB = cat(3, objSeedMask, zeros([size(objSeedMask), 2]));

imseriesmaskshow( imAdjusted, objSeedMask );
set( gcf, 'Name', 'Cell Seed Point Detection Result' );

%% Cell Segmentation 
fprintf( '\n\n>> Segmenting cells ...\n\n' );

imGradMag = callMatitkFilter( 'FGradientMagnitude', {}, imAdjusted );

switch PARAMETERS.SEG_ALGO

    case 'WATERSHED'

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
            imCurCellPerim = bwperim( imCurCellMask );
            imLabelCellSeg( imCurCellPerim > 0 ) = 0;
            
        end
        
end

% remove very small objects
minCellVolume = 0.5 * pi * (0.5 * min(PARAMETERS.cellDiameterRange))^2;
L = bwlabeln(L > 0);
regstats = regionprops(L, 'Area');
objAreaList = [regstats.Area] * prod(spacing);
smallObjInd = find( objAreaList < minCellVolume );
L( ismember( L, smallObjInd ) ) = 0;
imLabelCellSeg = bwlabeln( L > 0 );

numCells = max( imLabelCellSeg(:) );

[segMaskRGB, segColorMap] = label2rgbND( imLabelCellSeg );

[gx, gy] = gradient(imLabelCellSeg); 
imSegBoundaryMask = (gx ~= 0 | gy ~= 0);

imseriesmaskshow(imInput, {imSegBoundaryMask, objSeedMask});
set(gcf, 'Name', 'cell segmentation boundary mask overlay');

%% Post-processing to correct over-segmentation

    % Compute basic properties of cell regions
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
    
    % build region adjacency graph (RAG)
    adjMatrix = sparse(numCells, numCells);

    displayAdjacencyGraphX = [];
    displayAdjacencyGraphY = [];    
    neighStrel = streldisknd(PARAMETERS.neighSearchRad*ones(1, ndims(imInput)));

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

                if percentOverlap < PARAMETERS.threshNeighOverlap
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
    set(gcf, 'Name', 'cell neighborhood graph overlayed on segmentation label map');
    hold on;
    plot( displayAdjacencyGraphX, displayAdjacencyGraphY, 'go-', 'LineWidth', 2.0 );
    hold off;
    
    % initialize similarity matrix
    simMatrixInitial = adjMatrix;
    
    [v1, v2, ~] = find( adjMatrix );

    for k = (find(v1 < v2))'
       
        c1 = v1(k);
        c2 = v2(k);
        
        [ simVal ] = ComputeRegionSimilarity(imAdjusted, imLabelCellSeg, cellStats, c1, c2, wekaModel );
        
        simMatrixInitial(c1,c2) = simVal; 
        simMatrixInitial(c2,c1) = simVal;

    end
    
    imseriesmaskshow(imInput, {imSegBoundaryMask, objSeedMask});
    set(gcf, 'Name', 'cell neighborhood graph after pruning edges using a classifier');
    DisplayWeightedRegionAdjacencyGraph(gca, simMatrixInitial, cellStats );
    
    % iteratively merge edges that are locally maximal
    pq = PriorityQueue();       
    imLabelCellSegProcessed = imLabelCellSeg;
    regStats = cellStats;

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
        flagWasVertexMerged = zeros(numCells, 1);
        simMatrix = simMatrixInitial;
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

        regStats = cellStats;
        
        while ~pq.isEmpty()

            % find the most similar edge and merge it
            [curEdgeSimilarity, nodeIndices] = pq.pop();

            i = nodeIndices(1);
            j = nodeIndices(2);

            % check if any of nodes have been merged into another earlier
            if any( flagWasVertexMerged(nodeIndices) ) || simMatrix(i,j) ~= curEdgeSimilarity
                continue;
            end

            % display regions being merged
            figure(figMergeProgress);                
            curEdgeSimilarityNmzd = mat2gray(curEdgeSimilarity, [0.5, 1.0]);
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

            % decide which region becomes parent and child
            if regStats(i).Area > regStats(j).Area
                parentId = i;
                childId = j;
            else
                parentId = j;
                childId = i;
            end

            flagWasVertexMerged( childId ) = parentId;
            simMatrix(childId, :) = 0;
            simMatrix(:, childId) = 0;

            % compute merged region
            imRegionUnionMask = ismember( imLabelCellSegProcessed, [i, j] );
            imRegionUnionMask = imclose(imRegionUnionMask, neighStrel);            
            imLabelCellSegProcessed( imRegionUnionMask > 0 ) = parentId;

            % compute region properties
            regProps = ComputeRegionProperties(imLabelCellSegProcessed, parentId, spacing);

            indCurRegionSeedPoint = regProps.PixelIdxList( find( imCellSeedPoints(regProps.PixelIdxList) > 0 ) );
            if isempty( indCurRegionSeedPoint )            
                ptCurRegionSeedPoint = regProps.Centroid;
            else
                [~, maxind] = max( imCellSeedPoints( indCurRegionSeedPoint ) );
                indCurRegionSeedPoint = indCurRegionSeedPoint( maxind );    
                ptCurRegionSeedPoint = ind2submat( size(imCellSeedPoints), indCurRegionSeedPoint );
            end
            regProps.ptCellSeedPoint = ptCurRegionSeedPoint;
            regStats(parentId) = regProps;
            
            % recompute similarity to adjacent regions
            neighRegionIndices = setdiff( find(simMatrix(i , :) > 0 | simMatrix(j , :) > 0), [i,j] );

            for nid = neighRegionIndices

                if flagWasVertexMerged(nid)
                    error('error: code shouldnt come here');
                end

                curNeighSimVal = ComputeRegionSimilarity(imAdjusted, imLabelCellSegProcessed, cellStats, c1, c2, wekaModel );

                simMatrix(parentId, nid) = curNeighSimVal;
                simMatrix(nid, parentId) = curNeighSimVal;
                
            end

            for nid = neighRegionIndices

                if IsEdgeLocallyMaximal(simMatrix, parentId, nid)                               
                    pq.insert( full(simMatrix(parentId, nid)), sort([parentId, nid]) );
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

function DisplayWeightedRegionAdjacencyGraph( hAxis, simMatrix, regProps, flagShowMaximallySimilarEdgesOnly)

    if ~exist( 'flagShowMaximallySimilarEdgesOnly', 'var' )
        flagShowMaximallySimilarEdgesOnly = false;
    end
    
    axes(hAxis);

    maxLineWidth = 6.0;
    
    [v1, v2, simVal] = find( simMatrix );        
    simValNmzd = mat2gray(simVal, [0.5, 1]);

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
            
            lineColor = jetColorMap( 1 + round(simValNmzd(k) * (numColors-1)), : );

            hold on;
            plot( [regProps(i).ptCellSeedPoint(2), regProps(j).ptCellSeedPoint(2)], ...
                  [regProps(i).ptCellSeedPoint(1), regProps(j).ptCellSeedPoint(1)], ...
                  '-', 'Color', lineColor, 'LineWidth', 1.0 + maxLineWidth * simValNmzd(k) );
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

function [ simVal ] = ComputeRegionSimilarity(imAdjusted, imLabelCellSeg, cellStats, c1, c2, wekaModel )

    % compute region merging features
    featureDataStruct = ComputeRegionMergingFeatures(imAdjusted, imLabelCellSeg, cellStats, c1, c2);          
    [featureVec , featureNameList] = ConvertFeatureStructToFeatureVec( featureDataStruct );        

    % apply region merging classifier
    [predictedClassLabel, predictionProbability] = ApplyWekaRegionMergingClassifier(wekaModel, featureVec, featureNameList);

    if predictedClassLabel > 0
        simVal = predictionProbability;
    else
        simVal = 0; % this edge will be skipped by the merging process
    end
    
end
