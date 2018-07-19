clc
clear 
close all

% load data
stefanOldDataDir = 'C:\deepak\data\Stefan_June2012';
stefanNewDataDir = 'C:\deepak\data\Stefan_July2012';
ralphDataDir = 'C:\deepak\data\Ralph_June2012';
adhesionDir = 'C:\deepak\data\adhesions';

dataFileDir = stefanNewDataDir;

if ~exist( 'dataFilePath', 'var' )
    [fileName,pathName] = uigetfile( fullfile( dataFileDir, '*.mat' ), 'Select the data file' );   
    dataFilePath = fullfile( pathName, fileName )
end

dataFileContents = load( dataFilePath );
imInput = dataFileContents.im;
spacing = [0.5, 0.5];

% adjust the intensity range of the image 
imAdjusted = mat2gray(imInput) * 4096;    

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

%imAdjusted( ~imThresh ) = 0;

% compute image gradient
[gx, gy] = gradient( imAdjusted );
imGradMag = sqrt( gx.^2 + gy.^2 );

% compare different methods for cells seed point detection
seed_detection_type = 'automatic';

switch seed_detection_type
    
    case 'interactive'
        
        [ fgnd_seed_points , bgnd_seed_points ] = get_fgnd_bgnd_seeds_2d_strokes( imGradMag );
        
        SEED_NEIGH = 3;
        objSeedMask = zeros(size(imInput));
        for i = 1:numel(fgnd_seed_points)
            curSeedIndices = sub2ind(size(imInput), fgnd_seed_points{i}(:,2), fgnd_seed_points{i}(:,1)); 
            objSeedMask(curSeedIndices) = i;
        end
        objSeedMask = imdilate( objSeedMask, streldisknd(SEED_NEIGH*ones(1,ndims(imInput))) );
        objSeedIndices = find(objSeedMask);

        seedIndices = objSeedIndices;
        seedLabels = objSeedMask(objSeedIndices) - 1;
        numLabels = numel(fgnd_seed_points);
        
    case 'automatic'
        
        cellDiameterRange = [8, 20];

        imDistMap = bwdist(imThresh);
        imDistMap( imDistMap > 0 ) = 0;
        imDistMap = -1 * imDistMap;

        imseriesmaskshow( imDistMap, imThresh );
        set( gcf, 'Name', 'Foreground Distance Map' );    

        numLoGScales = min(10, cellDiameterRange(2)-cellDiameterRange(1)+1);
        [imCellSeedPoints, imResponse] = detectBlobsUsingAdaptiveMultiscaleLoG( imAdjusted, ...
                                                                                imDistMap, ... 
                                                                                'blobDiameterRange', cellDiameterRange, ...
                                                                                'numLoGScales', numLoGScales, ...
                                                                                'spacing', spacing, ...
                                                                                'debugMode', true );
        imCellSeedPoints(~imThresh) = 0; % remove seed points in background
        meanResponseBgnd = mean( imResponse(~imThresh) );
        stdResponseBgnd = std( imResponse(~imThresh) );
        imCellSeedPoints( imCellSeedPoints < (meanResponseBgnd + 2.5 * stdResponseBgnd) ) = 0; % remove seed points with weak LoG response

        objSeedMask = double(imdilate(imCellSeedPoints > 0, streldisknd([3,3])));
        
        imseriesmaskshow( imAdjusted, objSeedMask );
        set( gcf, 'Name', 'Cell Seed Point Detection Result' );

        erodeAmount = 5;
        stel = streldisknd(erodeAmount * ones(1, ndims(imInput)));
        imBackgroundSeeds = imerode(~imThresh, stel) | bwmorph( ~imThresh, 'shrink', 5 );
        bgndSeedIndices = find(imBackgroundSeeds);
        
        imForegroundSeeds = bwmorph(imCellSeedPoints > 0, 'thicken', 5);
        L = bwlabeln( imForegroundSeeds );
        fgndSeedIndices = find(L > 0);

        seedIndices = [fgndSeedIndices; bgndSeedIndices];
        seedLabels = [L(fgndSeedIndices); zeros(size(bgndSeedIndices))];
        numLabels = max(L(:)) + 1;
        
        imseriesmaskshow( imAdjusted, {imBackgroundSeeds, imForegroundSeeds} );
        set( gcf, 'Name', 'Foreground-background seeds' );
        
%         seedIndices = fgndSeedIndices;
%         seedLabels = L(fgndSeedIndices);
%         numLabels = max(L(:));
        
end

% run random walker segmentation algorithm
tic

% [segMask, ...
%  pixToLabelProbabilityMap] = random_walker_segmentation( mat2gray(imResponse), seedIndices, seedLabels, ...
%                                                           'sigma', 0.05, ...
%                                                          'flagNormalize', true);
% segMask(~imThresh) = 0;

[segMask, ...
 pixToLabelProbabilityMap] = random_walker_segmentation( imAdjusted, seedIndices, seedLabels, ...
                                                         'sigma', 0.05, ...
                                                         'flagNormalize', false);
                                                     
% [segMask, ...
%  pixToLabelProbabilityMap] = random_walker_segmentation( imAdjusted, seedIndices, seedLabels, ...
%                                                          'sigma', 0.05, ...
%                                                          'flagNormalize', false);

%[segMask, pixToLabelProbabilityMap] = random_walker( imInput, seedIndices, seedLabels );
%segMask(~imThresh) = 0;

timeElapsed = toc

segMaskRGB = label2rgbND(segMask);

% display result
objSeedMaskRGB = cat(3, objSeedMask > 0, zeros([size(objSeedMask), 2]));
imseriesmaskshowrgb(imInput, {segMaskRGB, objSeedMaskRGB} );
set( gcf, 'Name', 'Random Walker Segmentaion Result -- mask, object-seed-mask' );

imseriesmaskshowrgb(imGradMag, {segMaskRGB, objSeedMaskRGB} );

% if ndims(imInput) == 2
%     
%     segMaskRep = label2rgbND( repmat(segMask, [ones(1,ndims(segMask)), numLabels]) );
%     objSeedMaskRep = label2rgbND( repmat(objSeedMask, [ones(1,ndims(objSeedMask)), numLabels]) );
% 
%     imseriesmaskshowrgb(pixToLabelProbabilityMap * 100, {segMaskRep, objSeedMaskRep} ); 
%     set( gcf, 'Name', 'Random Walker Segmentaion Result -- pixelLabelProbabilityMap' );
%     
% end