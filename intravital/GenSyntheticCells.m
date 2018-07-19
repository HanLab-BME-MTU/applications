
clc
clear 
close all

im = zeros(500,500);

majorAxisRadiis = [ 20:25 ];
majorMinorRatios = [0.55:0.02:0.65];
cellOrientations = [0:20:360];
numCells = 20;

fgMeanVar = [ 240, 20 ; 220, 20 ];
bgMeanVar = [ 180, 20 ];

bgGmObj = gmdistribution( bgMeanVar(:,1), reshape(bgMeanVar(:,2), [1,1,size(bgMeanVar,1)]) );

imsize = size(im);

[X,Y] = meshgrid(1:imsize(2),1:imsize(1));

im(:) = random(bgGmObj,numel(im));
mask = zeros(size(im));

for i = 1:numCells
    
    pixind = randsample(prod(imsize),1);
    [yc, xc] = ind2sub( imsize, pixind );
    
%     yc = round(0.5 * imsize(1));
%     xc = round(0.5 * imsize(2));
    
    majorRad = randsample(majorAxisRadiis,1);
    minorRad = round( majorRad * randsample(majorMinorRatios,1) );
    theta = randsample(cellOrientations,1);
            
    pts = [ X(:) - xc, Y(:) - yc ];
    pts = pts * [ cosd(theta) -sind(theta); sind(theta), cosd(theta) ];
    
    curEllipseInd = find( (pts(:,1).^2 / majorRad^2)  + (pts(:,2).^2 / minorRad^2) - 1 <= 0 );
    
    % randomly select a foreground distrubution and draw samples from it
    fgdistind = randsample(size(fgMeanVar,1),1);
    fgGmObj = gmdistribution( fgMeanVar(fgdistind,1), fgMeanVar(fgdistind,2) );    
    im( curEllipseInd ) = random(fgGmObj,numel(curEllipseInd)); 
    %im( curEllipseInd ) = 255; 
    mask( curEllipseInd ) = 1;
end

mask = imdilate( mask, strel('disk',3) );

imseriesmaskshow( im, mask );

% imCellSeedPoints = GradientBasedCellCenterDetecter( im, [9,13], mask );
% imseriesshow( im );
% imseriesmaskshow_multichannel( { double(edge(im)), imCellSeedPoints } );

imCellSeedPoints = detect_cell_seeds_radial_symmetry( im, [15,25], ...
                                                      'roiMask', mask, ...                  
                                                      'numRadiusSamples', 10 );                                                        

imseriesmaskshow( im, imdilate(imCellSeedPoints, ones(3,3)) );
