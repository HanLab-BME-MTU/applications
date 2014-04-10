clc
clear 
close all

fgMeanVar = [ 200, 20 ];
bgMeanVar = [ 180, 20 ];
cellRadius = 10 * [1, 2];

im = zeros(200,200);
fgGmObj = gmdistribution( fgMeanVar(:,1), reshape(fgMeanVar(:,2), [1,1,size(fgMeanVar,1)]) );
bgGmObj = gmdistribution( bgMeanVar(:,1), reshape(bgMeanVar(:,2), [1,1,size(bgMeanVar,1)]) );

imsize = size(im);
fgMask = zeros( size(im) );
[X,Y] = meshgrid(1:imsize(1),1:imsize(2));

im(:) = random(bgGmObj,numel(im));
xc = 0.5 * imsize(2);
yc = 0.5 * imsize(1);
pts = [ X(:) -  xc, Y(:) - yc ];
curEllipseInd = find( (pts(:,1).^2 / (cellRadius(2))^2)  + (pts(:,2).^2 / (cellRadius(1))^2) - 1 <= 0 );    
im( curEllipseInd ) = random(fgGmObj,numel(curEllipseInd)); 
fgMask( curEllipseInd ) = 1;

% compare different methods for cells seed point detection   
meanCellDiameter = 2 * mean(cellRadius);
cellDiameterRange = 2 * ( mean(cellRadius) + [-5,5] );
imCellSeedDetectionResult = {};
strSeedDetectionAlgorithm = {};

    % Local Maxima in LoG Filtered Image
    fprintf( '\n\nDetecting seed points as local maxima in LoG Filtered Image ...\n\n' );

    strSeedDetectionAlgorithm{end+1} = 'FixedScaleLoG';
    imCellSeedDetectionResult{end+1} = detectBlobsUsingLoG( im, ...
                                                            meanCellDiameter, ...
                                                            'debugMode', true );

    % Local Maxima in a Multiscale LoG Filtered Image
    fprintf( '\n\nDetecting seed points as local maxima in Multiscale LoG Filtered Image ...\n\n' );

    strSeedDetectionAlgorithm{end+1} = 'MultiScaleLoG';    
    imCellSeedDetectionResult{end+1} = detectBlobsUsingMultiscaleLoG( im, ...
                                                                      cellDiameterRange, ...
                                                                      'debugMode', true );
                                                        
% suppress seed points outside the thresholded foreground region
for i = 1:numel( imCellSeedDetectionResult )
    imCellSeedDetectionResult{i}( ~fgMask ) = 0;
end

% display result
imSeedMaskForDisplay = imCellSeedDetectionResult;
for i = 1:numel( imSeedMaskForDisplay )
    imSeedMaskForDisplay{i} = imdilate( imSeedMaskForDisplay{i}, strel('disk', 3) );
end
imseriesmaskshow( im, imSeedMaskForDisplay );
set( gcf, 'Name', [ 'Cell Seed Point Detection Results', ...
                    sprintf( ' -- %s', strSeedDetectionAlgorithm{:} ) ] );
