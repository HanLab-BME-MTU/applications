function [reg_grid,xvec,yvec,gridSpacing]=createRegGridFromDisplField(displField,mag)
% For Benedikt's software to work, the displacement field has to be
% interpolated on a rectangular grid, with an even number of grid points
% along each edge. Furthermore, one has to make sure that the noisy data 
% has not to be extrapolated. This may happen along the edges. To prevent
% this, extract the corner of the displacement grid, calculate how often 
% (even number) the optimal gridspacing fits into each dimension, then 
% place the regular grid centered to the orignal bounds. Thereby make sure 
% that the edges have been eroded to a certain extend.
if nargin < 2
    mag=1;
end
% First get the corners (here we take field no.1 but that is arbitrary):
LeftUpperCorner(1:2) = [min(displField(1).pos(:,1)), min(displField(1).pos(:,2))];
RightLowerCorner(1:2) = [max(displField(1).pos(:,1)), max(displField(1).pos(:,2))];
oldSize=RightLowerCorner-LeftUpperCorner;

% Next get the average inter bead distance (here we take field no.1 but that is arbitrary):
numBeads=length(displField(1).pos);
areaImg=prod(RightLowerCorner-LeftUpperCorner);
interBeadDist=sqrt(areaImg/numBeads);

%display('Interbead Distance has been upscaled in createRegGridFromDisplField. This has to be taken out again:')
%interBeadDist=sqrt(2)*interBeadDist;

% See how often the optimal grid spacing fits into the bounds:
gridSpacing=ceil(interBeadDist/mag);
numPieces=floor((RightLowerCorner-LeftUpperCorner)/gridSpacing);

% Make sure that we get an odd Number of pieces which leads to an even 
% number of grid points later:
numPieces=numPieces-(1-mod(numPieces,2));

% The size of the new grid is then given by:
newSize=numPieces*gridSpacing;

% If there is not much difference to the old size, we need to further erode
% the edges (this part is optional and could be skipped):

edgeErode=1;
if edgeErode==1
    sizeDiff=oldSize-newSize;
    if sizeDiff(1)<gridSpacing
        numPieces(1)=numPieces(1)-2;
    end
    if sizeDiff(2)<gridSpacing
        numPieces(2)=numPieces(2)-2;
    end
    newSize=numPieces*gridSpacing;
end

% Now place the new grid in the center of the old bounds:
sizeDiff=oldSize-newSize;
newLeftUpperCorner=LeftUpperCorner+floor(sizeDiff/2);
newRightLowerCorner=newLeftUpperCorner+numPieces*gridSpacing;

xvec=newLeftUpperCorner(1):gridSpacing:newRightLowerCorner(1);
yvec=newLeftUpperCorner(2):gridSpacing:newRightLowerCorner(2);

% Now we have everything needed for making the new grid:
[X,Y] = meshgrid(xvec, yvec);
reg_grid(:,:,1) = X';
reg_grid(:,:,2) = Y';

xvec=reshape(reg_grid(:,:,1),[],1);
yvec=reshape(reg_grid(:,:,2),[],1);