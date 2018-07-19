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
    PARAMETERS.cellDiameterRange = [6, 20];
    
    % less-critical parameters (leave defaults)
    PARAMETERS.threshNeighOverlap = 5;
    PARAMETERS.neighSearchRad = 2;    
    PARAMETERS.logResponseCutoff = 1.0;
    
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

imThresh = imfill(imThresh);                                                
imAdjusted( ~imThresh ) = 0;

% compute cell seed points
imDistMap = callMatitkFilter( 'FSignedMaurerDistancemap', {} , imThresh, [], [], spacing  );
imDistMap( imDistMap > 0 ) = 0;
imDistMap = -1 * imDistMap;

numLoGScales = min(20, PARAMETERS.cellDiameterRange(2)-PARAMETERS.cellDiameterRange(1)+1);
[imCellSeedPoints, imResponse] = detectBlobsUsingAdaptiveMultiscaleLoG( imAdjusted, ...
                                                                        imDistMap, ... 
                                                                        'blobDiameterRange', PARAMETERS.cellDiameterRange, ...
                                                                        'numLoGScales', numLoGScales, ...
                                                                        'spacing', spacing, ...
                                                                        'debugMode', false );
                                                                    
imCellSeedPoints(~imThresh) = 0; % remove seed points in background
meanResponseBgnd = median( imResponse(~imThresh) );
stdResponseBgnd = 1.4826 * mad(imResponse(~imThresh), 1);
imCellSeedPoints( imCellSeedPoints < (meanResponseBgnd + PARAMETERS.logResponseCutoff * stdResponseBgnd) ) = 0; % remove seed points with weak LoG response

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
minCellVolume = 0.5 * pi * (0.5 * min(PARAMETERS.cellDiameterRange))^2;
L = bwlabeln(L > 0);
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
set(gcf, 'Name', 'cell neighborhood graph overlayed on gradient magnitude');
hold on;
plot( displayAdjacencyGraphX, displayAdjacencyGraphY, 'go-', 'LineWidth', 2.0 );
hold off;

imseriesmaskshow(imLabelCellSeg, {imSegBoundaryMask, objSeedMask});
set(gcf, 'Name', 'cell neighborhood graph overlayed on gradient magnitude');
hold on;
plot( displayAdjacencyGraphX, displayAdjacencyGraphY, 'go-', 'LineWidth', 2.0 );
hold off;

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

hVisAnnotation = figure;
imMaskOverlay = genImageMaskOverlay(imInput, {imSegBoundaryMask, objSeedMask}, [1 0 0; 0 1 0], [0.5, 0.2] );
imshow(imMaskOverlay);
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

% save annotation
annotationData.PARAMETERS = PARAMETERS;
annotationData.dataFilePath = dataFilePath;
annotationData.imInput = imInput;
annotationData.imThresh = imThresh;
annotationData.imCellSeedPoints = imCellSeedPoints;
annotationData.imLabelCellSeg = imLabelCellSeg;
annotationData.numCells = numCells;
annotationData.adjMatrix = adjMatrix;
annotationData.adjMatrixMerge = adjMatrixMerge;

[pathName, fileName, ext] = fileparts(dataFilePath);
outFileName_Prefix = sprintf('%s_%d_%d', fileName, PARAMETERS.cellDiameterRange(1), PARAMETERS.cellDiameterRange(2));
[outFileName, outFilePath] = uiputfile( fullfile(pathName, [outFileName_Prefix, '.rann']) );
save( fullfile(outFilePath, outFileName), 'annotationData' );
SaveFigure(hVisAnnotation, fullfile(outFilePath, [outFileName_Prefix '_VisAnnotation.png']), 'png');
