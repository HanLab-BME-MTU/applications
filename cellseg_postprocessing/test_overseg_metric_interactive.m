clc
clear 
close all;

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


% detect seed points
cellDiameterRange = [8, 12];

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
                                                                    'debugMode', false );
imCellSeedPoints(~imThresh) = 0;

objSeedMask = double(imdilate(imCellSeedPoints > 0, streldisknd([3,3])));

imseriesmaskshow( imResponse, objSeedMask );

imseriesmaskshow( imAdjusted, objSeedMask );
set( gcf, 'Name', 'Cell Seed Point Detection Result' );

% get line end points to examine
[x, y] = getpts();
close(gcf);

v1 = [x(1), y(1)];
v2 = [x(2), y(2)];

vmag = norm( v2 - v1 );

t = linspace(0, 1.0, vmag);
ptLine = [x(1) + t(:) * (x(2) - x(1)), y(1) + t(:) * (y(2) - y(1))];

% read intensity profile between points
figure;

subplot(3,1,1), plot( mat2gray( interp2( imInput, ptLine(:,1), ptLine(:,2) ) ) );
ylim( [0, 1] );
title( 'input image profile ' );

subplot(3,1,2), plot( mat2gray( interp2( imAdjusted, ptLine(:,1), ptLine(:,2) ) ) );
ylim( [0, 1] );
title( 'denoised image profile ' );

subplot(3,1,3), plot( mat2gray( interp2( imResponse, ptLine(:,1), ptLine(:,2) ) ) );
ylim( [0, 1] );
title( 'multi-scale LoG profile' );


