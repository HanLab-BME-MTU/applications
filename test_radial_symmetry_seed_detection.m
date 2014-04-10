clc
clear 
close all

im = zeros(200,200);
cellRadius = [10, 20];
overlap = 0.0;
theta = [0:90:359];

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

roiMask = imdilate( fgMask, strel('disk',3) );
imseriesmaskshow( im, roiMask );

imCellSeedPoints = detect_cell_seeds_radial_symmetry( im, [15 15], ...
                                                      'roiMask', roiMask, ...                  
                                                      'numRadiusSamples', 1 );                                                        

imseriesmaskshow( im, imdilate(imCellSeedPoints, ones(3,3)) );
