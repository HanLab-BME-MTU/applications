function  [spotsidx, mask] = discernspots(cordList,mSize,dataProperties);
%DISCERNSPOTS distinguish between partially overlapping spots and isolated spots
%
% SYNOPSIS [spotsidx, mask] = discernspots(cordList,mSize);
%
% INPUT cordList : nx3 matrix (nspots in 3 dimsension)
%            mSze     : Size of mask [x y z]
%
% OUTPUT mask      :  binary mask matrix
%                spotdsidx :  index into cordList to overlapping spots

% c: 10/10/01	dT

%Define constants
PIXELSIZE_XY=dataProperties.PIXELSIZE_XY;
PIXELSIZE_Z=dataProperties.PIXELSIZE_Z;
FT_SIGMA=dataProperties.FT_SIGMA;

%take first check distance to others 
%dm=distMat(cordList,diag([PIXELSIZE_XY^2 PIXELSIZE_XY^2 PIXELSIZE_Z^2]));
dmXY=distMat(cordList(:,1:2),diag([PIXELSIZE_XY^2 PIXELSIZE_XY^2]));
dmZ=distMat(cordList(:,3),diag([PIXELSIZE_Z^2]));

spotsidx=rec_find(dmXY,dmZ,1,PIXELSIZE_XY, PIXELSIZE_Z, FT_SIGMA);

%create mask
mask=zeros(mSize);
% default sphere
ellipsoid=fillSphereData(5*FT_SIGMA(1),1);
for i = spotsidx
    cen=round(cordList(i,:));
    %      center=[cen(2) cen(1) cen(3)];
    mask=mask | pushstamp3d(mask,ellipsoid,cen);
end;


function spidx=rec_find(dMatrixXY,dMatrixZ,curIdxList,PIXELSIZE_XY, PIXELSIZE_Z, FT_SIGMA)
% recursive search for points which are separated by a distance smaller than MIN_DIST


MIN_DIST_XY=7*FT_SIGMA(1)*PIXELSIZE_XY;
MIN_DIST_Z=7*FT_SIGMA(3)*PIXELSIZE_Z;
%if smaller than d add spot
spidx=curIdxList;
for i=1:size(dMatrixXY,2)
    if (i~=curIdxList(end) & dMatrixXY(curIdxList(end),i)<MIN_DIST_XY & dMatrixZ(curIdxList(end),i)<MIN_DIST_Z)
        if (isempty(spidx) | ~any(spidx==i))
            spidx(end+1)=i;
            spidx=rec_find(dMatrixXY,dMatrixZ,spidx,PIXELSIZE_XY, PIXELSIZE_Z, FT_SIGMA);
        end;
    end;
end;
