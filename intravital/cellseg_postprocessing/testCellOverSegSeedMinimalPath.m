clc
clear 
%close all

%***********************************************************************
% 
%                           PARAMETERS 
%
%***********************************************************************

    dataFileDir = 'C:\deepak\data\Stefan_July2012';
    
    % critical parameters (needs tuning)
    cellDiameterRange = [6, 20];
    
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
                                                
if true % strcmp(SEG_ALGO, 'WATERSHED')
    imAdjusted( ~imThresh ) = 0;
end

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
                                                                        'debugMode', true );
                                                                    
imCellSeedPoints(~imThresh) = 0; % remove seed points in background
meanResponseBgnd = mean( imResponse(~imThresh) );
stdResponseBgnd = std( imResponse(~imThresh) );
imCellSeedPoints( imCellSeedPoints < (meanResponseBgnd + 2.0 * stdResponseBgnd) ) = 0; % remove seed points with weak LoG response

objSeedMask = double(imdilate(imCellSeedPoints > 0, streldisknd([3,3])));
objSeedMaskRGB = cat(3, objSeedMask, zeros([size(objSeedMask), 2]));

% check path profiles
[ imMaskOverlay ] = genImageRGBMaskOverlay( imAdjusted, objSeedMaskRGB, [0.5, 0.5] );
figure, imshow( imMaskOverlay, [] );
[x, y] = getpts(gcf);
ptSeed = [x, y];
ptSeedInd = sub2ind(size(imAdjusted), round(ptSeed(:,2)), round(ptSeed(:,1))); 

pathMask = zeros(size(imAdjusted));

numQuantLevels = 64;
imAdjustedQuant = round( mat2gray(imAdjusted) * (numQuantLevels-1) );
costMatrix = (1 - mat2gray(imDistMap)) + (1 - mat2gray(imAdjusted));

hProfile = figure;
ylim([0, numQuantLevels]);
set( gcf, 'Name', 'Interseed intensity profiles' );

hUniform = figure;
plot(0:0.1:1,0:0.1:1, 'r--', 'LineWidth', 2.0);
set( gcf, 'Name', 'Deviation from a uniform profile' );
ylim([0, 1]);

imseriesmaskshow( imAdjusted, objSeedMask );
hSeedNeighbors = gcf;
set( gcf, 'Name', 'Requested seed profiles' );

for i = 1:2:numel(ptSeedInd)

    figure( hSeedNeighbors );
    hold all;
    plot( [ptSeed(i,1), ptSeed(i+1,1)], ...
          [ptSeed(i,2), ptSeed(i+1,2)], ...  
          '-', 'LineWidth', 2.0 );
    hold off;

    srcNode = ptSeedInd(i);
    destNode = ptSeedInd(i+1);

    [path, cost] = findMinimalPathOnCostMatrix(costMatrix, srcNode, destNode, ...
                                               'roiMask', imThresh, 'debugMode', false, ...
                                               'costWeight', 0.5);
    
    pathMask(path) = 1;
    pathProfile = imAdjustedQuant(path);
    
    figure( hProfile );
    hold all;
    x = (0:numel(path)-1) / numel(path);
    pathMedianIntensity = median(pathProfile);
    plot(x, pathProfile, '-o'); 
    hold off;

    figure( hUniform );
    hold all;
    pathDist = pathProfile / sum(pathProfile);
    plot(x, cumsum(pathDist), '-'); 
    hold off;
    
end

imseriesmaskshow( imAdjusted, {imThresh, objSeedMask, imdilate(pathMask, streldisknd(ones(1,ndims(imAdjusted))))} );

