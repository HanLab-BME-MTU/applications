% 
% The following script tests the steerableDetector on a synthetic 2D image
% with lines drawn in random orientations
% 
% Author: Deepak Roy Chittajallu (Created on Aug 31, 2012)
%

clc
clear
close all

%*************************************************************************
%                               PARAMETERS
%*************************************************************************

    numLines = 5;
    
    meanLineWidth = 2;
    stdLineWidth = 1;
    
    orderSteerableDetector = 4;

    FeatureSpaceType = 'perpvec';
    
%*************************************************************************

im = zeros( 512, 512 );
imsize = size( im );

% normal distributons for foreground and background pixels
fgMeanVar = [ 200, 20 ];
bgMeanVar = [ 180, 20 ];
fgGmObj = gmdistribution( fgMeanVar(:,1), reshape(fgMeanVar(:,2), [1,1,size(fgMeanVar,1)]) );
bgGmObj = gmdistribution( bgMeanVar(:,1), reshape(bgMeanVar(:,2), [1,1,size(bgMeanVar,1)]) );
im(:) = random(bgGmObj,numel(im));

% generate lines randomly
fprintf( '\nGenerating %d Lines in Random Orientations: \n', numLines );

[X,Y] = meshgrid(1:imsize(2), 1:imsize(1));

try
    ptIm = gpuArray( [ X(:), Y(:) ] );
catch err
    ptIm = [ X(:), Y(:) ]; % you probably dont have a gpu, the code will be slow
end

lineMask = false( imsize );
lineWidthVec = zeros( numLines, 1 );
pixelClusterMap = zeros( imsize );

for i = 1:numLines   
    
    fprintf( ' %.3d ', i );
    
    if mod( i, 15 ) == 0
        fprintf( '\n' );
    end
    
    % generate a random point in the image
    ptRandInd = floor( rand * numel(im) );
    ptRef = ptIm( ptRandInd, : );
    
    % generate random orientation vector
    thetax = rand * pi;
    v = [ cos(thetax), sin(thetax) ];
    
    % generate random line width
    randLineWidth = meanLineWidth + stdLineWidth * abs( randn );
    
    if randLineWidth <= 0 
        randLineWidth = meanLineWidth;
    end
    
    lineWidthVec(i) = randLineWidth;
    
    % generate line
    refvec = ptIm - repmat( ptRef, numel(im), 1 );
    perp = refvec - (refvec * v') * v;
    sqDist = sum( perp .* perp, 2 );
    flagCurLineMask = sqDist <= randLineWidth^2;
    lineMask( flagCurLineMask ) = true;
    pixelClusterMap( flagCurLineMask ) = i;
    
end

fprintf( '\n' );

im( lineMask ) = random(fgGmObj, numel( find( lineMask ) ));
%im( lineMask ) = 1;
[imTrueRGBLineMask, clusterColorMap] = label2rgbND( pixelClusterMap );

% Run steerable detector to enhance the lines
fprintf( '\nRunning steerable detector to enhance curves on %d x %d sized image ...\n', imsize(2), imsize(1) );
sigmaTrialValues = (meanLineWidth + stdLineWidth * (0:0.5:2.5));
sigmaTrialValues( sigmaTrialValues <= 0 ) = [];

for i = 1:numel( sigmaTrialValues )
    
    fprintf( '\n\t%d/%d: Trying sigma value of %.2f ... ', i, numel( sigmaTrialValues ), sigmaTrialValues(i) );   
    
    tic
    [curRes, curTheta, curNms, rotations] = steerableDetector(im, orderSteerableDetector, sigmaTrialValues(i));
    timeElapsed = toc;
    
    fprintf( 'It took %.2f seconds\n', timeElapsed );   
    
    curRes = sigmaTrialValues(i)^2 * curRes; % scale normalization
    
    if i == 1        
        res = curRes;
        nms = curNms;
        theta = curTheta;
        pixelScaleMap = ones( size(res) );
    else
        indBetter = curRes > res;
        res(indBetter) = curRes(indBetter);
        nms(indBetter) = curNms(indBetter);
        theta(indBetter) = curTheta(indBetter);
        pixelScaleMap(indBetter) = i;
    end
end

theta = pi/2.0 + theta;
imThreshRes = res > thresholdOtsu(res);
imLineSegMask = imThreshRes > 0 & nms > 0;
imLineSegRGBMask = zeros( [size(imLineSegMask), 3] );
imLineSegRGBMask(:,:,1) = imLineSegMask;

pixelScaleMap( ~imLineSegMask ) = 0;
imScaleRGB = label2rgb( pixelScaleMap, 'jet', 'k' ) / 255.0;

% display
imseriesshow( im );
set( gcf, 'Name', 'Image with Randomly Generated Lines' );

imseriesmaskshowrgb( res, {imLineSegRGBMask, imScaleRGB, imTrueRGBLineMask} );
colorbar;
set( gcf, 'Name', 'Response of Steerable Detector with Otsu Line Mask and Pixel Scale Map' );

imseriesmaskshowrgb( nms, {imLineSegRGBMask, imScaleRGB, imTrueRGBLineMask} );
colorbar;
set( gcf, 'Name', 'Result of Non-maximal suppression with Otsu Line Mask and Pixel Scale Map' );

thetaDisplay = round( 1 + rad2deg(theta) );
thetaDisplay( ~imThreshRes ) = 0;
thetaDisplayRGB = label2rgb( thetaDisplay, 'hot', 'k' ) / 255.0;
imseriesmaskshowrgb( im, thetaDisplayRGB );
set( gcf, 'Name', 'Orientation Map' );

figure, imshow( res, [] );
hold on;
thetaSmooth = medfilt2( theta, [5,5] );
vx = cos(thetaSmooth);
vy = sin(thetaSmooth);
linind = find( imLineSegMask );
[yind, xind] = ind2sub( imsize, linind );
quiver( xind, yind, vx( linind ), vy( linind ) );
hold off;
set( gcf, 'Name', 'Orientation Vectors at Ridge Maxima' );

% w = 5;
% pixelClusterMap = zeros( size(im) );
% loc = [ 100:100:400 ];
% for i = 1:numel( loc )
% 
%     lineMask1 = false( size(im) );
%     lineMask1( loc(i):(loc(i)+w-1), : ) = true;
%     pixelClusterMap(lineMask1) = 2 * i - 1;
% 
%     lineMask2 = false( size(im) );
%     lineMask2( :, loc(i):(loc(i)+w-1) ) = true;
%     pixelClusterMap(lineMask2) = 2 * i;
%     
%     im( lineMask1 | lineMask2 ) = 1;
%     
% end
% 
% [imTrueRGBLineMask, clusterColorMap] = label2rgbND( pixelClusterMap );
% 
% [res, theta, nms, rotations] = steerableDetector(im, orderSteerableDetector, w / sqrt(2));
% 
% imThreshRes = res > thresholdOtsu(res);
% imLineSegMask = imThreshRes & nms > 0;
% imLineSegRGBMask = zeros( [size(imLineSegMask), 3] );
% imLineSegRGBMask(:,:,1) = imLineSegMask;
% 
% % display
% imseriesshow( im );
% set( gcf, 'Name', 'Image with Randomly Generated Lines' );
% 
% imseriesmaskshowrgb( res, {imLineSegRGBMask, imTrueRGBLineMask} );
% set( gcf, 'Name', 'Response of Steerable Detector with Otsu Line Mask and Pixel Scale Map' );
% 
% imseriesmaskshowrgb( nms, {imLineSegRGBMask, imTrueRGBLineMask} );
% set( gcf, 'Name', 'Result of Non-maximal suppression with Otsu Line Mask and Pixel Scale Map' );
% 
% imseriesshow( theta );
% set( gcf, 'Name', 'Orientation Map' );
% 
% figure, imshow( res, [] );
% hold on;
% vx = cos(theta);
% vy = sin(theta);
% [yind, xind] = ind2sub( imsize, find( imLineSegMask ) );
% quiver( xind, yind, vx( xind ), vy( yind ) );
% hold off;
% set( gcf, 'Name', 'Orientation Vectors at Ridge Maxima' );

% map line mask points to feature space
flagLineMask = imThreshRes > 0;
[yind, xind] = ind2sub( size( imLineSegMask ), find( flagLineMask ) );
ptLineMask = [ xind, yind];
lineOrientationVec = [ cos(theta(flagLineMask)), sin(theta(flagLineMask)) ];
perpvec = ptLineMask - repmat( sum( ptLineMask .* lineOrientationVec, 2 ), 1, 2 ) .* lineOrientationVec;

switch FeatureSpaceType
    
    case 'perpvec'
        
        ptLineFeature = perpvec + 5 * rand( size( perpvec ) );
        % ptLineFeature = perpvec; 

        figure;
        hold on;
        for i = 1:size( ptLineMask, 1 )
            curPtClusterId = pixelClusterMap( round( ptLineMask(i,2) ), round( ptLineMask(i,1) ) );
            if ~curPtClusterId
                continue;
            end
            plot( ptLineFeature(i,1), ptLineFeature(i,2), '.', 'Color', clusterColorMap(curPtClusterId, :) );
        end
        hold off;
        
        set( gcf, 'Name', 'Curvy Points in Feature Space' );
        title( 'Curvy Points in Feature Space' );
        xlabel( 'x-coordinate of perpendicular vector', 'FontWeight', 'Bold' );
        ylabel( 'y-coordinate of perpendicular vector', 'FontWeight', 'Bold' );
        axis( [0, 512, 0, 512] );
        
    case 'parametric'

        perpvec = perpvec + 5 * rand( size( perpvec ) );
        % perpdist = sqrt( sum( perpvec .* perpvec, 2 ) );
        ptLineFeature = [ atan2( perpvec(:,2), perpvec(:,1) ), perpdist ];
        
        figure;
        hold on;
        for i = 1:size( ptLineMask, 1 )
            curPtClusterId = pixelClusterMap( round( ptLineMask(i,2) ), round( ptLineMask(i,1) ) );
            if ~curPtClusterId
                continue;
            end
            plot( ptLineFeature(i,1), ptLineFeature(i,2), '.', 'Color', clusterColorMap(curPtClusterId, :) );
        end
        hold off;
        
        set( gcf, 'Name', 'Curvy Points in Feature Space' );
        title( 'Curvy Points in Feature Space' );
        xlabel( 'Angle between perpendicular line and x-axis', 'FontWeight', 'Bold' );
        ylabel( 'Perpendicular distance to origin', 'FontWeight', 'Bold' );
        axis( [-pi, pi, 0, norm( size(im) )] );
        
end


% % run mean shift
% bandwidth = 5;
% [clusterInfo, pointToClusterMap] = MeanShiftClustering(ptLineFeature, bandwidth, ... 
%                                                          'flagDebug', true, ...
%                                                          'kernel', 'gaussian', ...
%                                                          'method', 'optimized', ...
%                                                          'minClusterDistance', bandwidth, ...
%                                                          'flagUseKDTree', true );
% 
% % display segmentation result
% figure;
% hold on;
% for i = 1:max(imAdhesionSegLabel(:))
%    ptClusterCenter = clusterInfo(i).ptClusterCenter;
%    plot( ptClusterCenter(1), ptClusterCenter(2), 'o', 'Color', [0,0,0], 'MarkerFaceColor', labelColorMap(i,:), 'LineWidth', 2.0 ); 
% end
% hold off;
% set( gcf, 'Name', 'Clustering Result' );
