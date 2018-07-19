
clc
clear all
close all

im = zeros(200,200);
cellRadius = [10, 20];
overlap = -0.2;
theta = [0:45:359];
divergenceThreshold = 1.0;

imsize = size(im);
minObjectDiameter = 2 * min(cellRadius);

fgMeanVar = [ 200, 20 ];
bgMeanVar = [ 180, 20 ];

fgGmObj = gmdistribution( fgMeanVar(:,1), reshape(fgMeanVar(:,2), [1,1,size(fgMeanVar,1)]) );
bgGmObj = gmdistribution( bgMeanVar(:,1), reshape(bgMeanVar(:,2), [1,1,size(bgMeanVar,1)]) );

[X,Y] = meshgrid(1:imsize(1),1:imsize(2));

im(:) = random(bgGmObj,numel(im));
fgMask = zeros( size(im) );
%im(:) = bgMeanVar(1);

for i = 1:numel(theta)    
    xc = 0.5 * imsize(2) + max(cellRadius) * (1 - overlap) * cosd( theta(i) );
    yc = 0.5 * imsize(1) + max(cellRadius) * (1 - overlap) * sind( theta(i) );
    
    pts = [ X(:) -  xc, Y(:) - yc ];
    pts = pts * [ cosd(theta(i)) -sind(theta(i)); sind(theta(i)), cosd(theta(i)) ];
    
    curEllipseInd = find( (pts(:,1).^2 / (cellRadius(2))^2)  + (pts(:,2).^2 / (cellRadius(1))^2) - 1 <= 0 );    
    
    % randomly select a foreground distrubution and draw samples from it
    im( curEllipseInd ) = random(fgGmObj,numel(curEllipseInd)); 
    %im( curEllipseInd ) = 1; 
    fgMask( curEllipseInd ) = 1;
    
end

[imCellSeedPoints] = detect_cell_seeds_gvf_blurred_image( im, minObjectDiameter, ... 
                                                          'debugMode', true, ... 
                                                          'maxGVFIterations', 50 );
imCellSeedPoints( ~fgMask ) = 0;

imseriesmaskshow( im, fgMask );
hold on;
[yind,xind] = ind2sub( size(im), find( imCellSeedPoints > 0 ) );
plot( xind, yind, 'g+', 'MarkerSize', 10.0 );
hold off;
set( gcf, 'Name', 'Cell Seed Point Detection Result' );

% figvec = [];
% 
% imseriesmaskshow( im, fgMask );
% set( gcf, 'Name', 'Test Image' );
% figvec(end+1) = gcf;
% 
% imBlurred = filterGauss2D(im, 3);
% [gx,gy] = gradient( imBlurred );
% gmag = sqrt( gx .* gx + gy .* gy );
% 
% imFeature = gmag;
% 
% imseriesmaskshow( imFeature, fgMask );
% [fx,fy] = gradient( imFeature );
% fgmag = sqrt( fx .* fx + fy .* fy );
% hold on, quiver(fx,fy), hold off
% set( gcf, 'Name', 'Initial Vector Field' );
% figvec(end+1) = gcf;
% 
% profile on;
% %[vx,vy] = GVFND(imFeature, 0.15, 500, 'flagDebugMode', true);
% [vx,vy] = gvfc(imFeature, 0.15, 500);
% profile off;
% profile viewer;
% 
% vmag = sqrt( vx .* vx + vy .* vy );
% imseriesmaskshow( imFeature, fgMask );
% hold on, quiver(vx ./ vmag, vy ./ vmag), hold off
% set( gcf, 'Name', 'Diffused Vector Field' );
% figvec(end+1) = gcf;
% 
% imdiv = divergence(vx ./ vmag, vy ./ vmag);
% imseriesmaskshow( imdiv, double( abs(imdiv) > 0.1 ) );
% set( gcf, 'Name', 'Divergence of Diffused Vector Field' );
% figvec(end+1) = gcf;
% 
% % locate local intensity maxima in divergence image
% MaximaSuppressionSize = round( minObjectDiameter/1.5 );
% evenind = (mod( MaximaSuppressionSize, 2 ) == 0);
% MaximaSuppressionSize( evenind ) = MaximaSuppressionSize( evenind ) + 1;    
% imLocalMax = locmax2d( abs(imdiv), MaximaSuppressionSize);    
% imLocalMax(~fgMask) = 0;
% 
% for i = 1:numel(figvec)
%     
%     figure( figvec(i) );
%     hold on;
%     [yind,xind] = ind2sub( size(im), find( imLocalMax > 0 ) );
%     plot( xind, yind, 'r+', 'MarkerSize', 20.0, 'LineWidth', 2.0 );
% 
%     [yind,xind] = ind2sub( size(im), find( imLocalMax > 0 & imdiv > divergenceThreshold ) );
%     plot( xind, yind, 'go', 'MarkerSize', 20.0, 'LineWidth', 2.0 );
%     hold off;
%     
% end
% 
% imseriesmaskshow( im, fgMask );
% hold on;
% [yind,xind] = ind2sub( size(im), find( imLocalMax > 0 & imdiv > divergenceThreshold ) );
% plot( xind, yind, 'g+', 'MarkerSize', 20.0, 'LineWidth', 2.0 );
% plot( xind, yind, 'bo', 'MarkerSize', 20.0, 'LineWidth', 2.0 );
% hold off;
% set( gcf, 'Name', 'Cell Seed Point Detection Result' );
