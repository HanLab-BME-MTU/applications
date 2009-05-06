function [mask] = makeCellMaskDetections(positions,closureRadius,dilationRadius,doFill,imSize,plotMask,overlayImage);

% makeCellMaskDetections calculates mask based on positions of particles by
% dilating, closing, and filling these particles
%
% INPUT:    positions = matrix where first column contains the x positions
%               of all points and second column contains the y positions

%           closureRadius =   radius for disk used in closure

%           dilationRadius = radius for dialation

%           doFill = 1 to fill mask, 0 to not fill

%           imSize = size of frame, [y x]

%           plotMask = 1 to display mask, 0 to not display

%           overlayImage = image on which to overlay mask; can be obtained
%               using imread(imagePath)
% OUTPUT
%           mask
%
% Uses:
%
%
% Daniel Nunez, updated May 04, 2009


%INTERPRET INPUTS

%POSITIONS
%positions must have 3 or more points
if length(positions) < 3
    error('must have more than three particles')
end
%get size of positions matrix
sizePositions = size(positions);
%since positions should include only x and y positions, one element of size must be equal
%to two; find which direction of the matrix is equal to two
direction = find(sizePositions == 2);
if direction == 1
    positions = positions';
end


%MAKE MASK OUT OF POINTS
mask = zeros(imSize(1),imSize(2));
%fill in pixels that have a point with a one
xpos = nonzeros(ceil(positions(:,1)));
ypos = nonzeros(ceil(positions(:,2)));
for ipoint = 1:length(xpos)
    mask(ypos(ipoint),xpos(ipoint)) = 1;
end

%pad image so that closure functions properly near edges
mask = padarray(mask,[2*closureRadius 2*closureRadius],0);

%PERFORM DILATION
seDilate = strel('disk',dilationRadius);
mask = imdilate(mask,seDilate);

%PERFORM CLOSURE
seClose = strel('disk',closureRadius);
mask = imclose(mask,seClose);

%GRAB LARGEST CONNECTED AREA AS MASK
[mask,labelNum] = bwlabeln(mask);
for ilabel = 1:labelNum
    labelSize(ilabel) = length(find(mask == ilabel));
end
maxArea = find(labelSize == max(labelSize));
mask(mask ~= maxArea) = 0;

%PERFORM FILL
if doFill
    mask = imfill(mask);
end

%remove padding from image edge
mask(1:2*closureRadius,:) = [];
mask(:,1:2*closureRadius) = [];
mask(1+imSize(1):end,:) = [];
mask(:,1+imSize(2):end) = [];


%PLOT MASK
if plotMask
    boundaries = bwboundaries(mask);
    for iter = 1:2
        figure
        if exist(overlayImage,'var') && ~isempty(overlayImage)
            imagesc(overlayImage)
            colormap(gray)
        end
        hold on
        axis equal
        if iter == 1
            plot(positions(:,1),positions(:,2),'r.','MarkerSize',0.2)
        end
        for ibound = 1:length(boundaries)
            bounds = boundaries{ibound};
            plot(bounds(:,2),bounds(:,1),'g')
        end
    end
end

end %of function