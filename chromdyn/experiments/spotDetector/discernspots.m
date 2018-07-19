function  [spotsidx, mask] = discernspots(cordList,mSize,dataProperties)
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
FILTERPRM=dataProperties.FILTERPRM;
mask = [];

%take first check distance to others
%dm=distMat(cordList,diag([PIXELSIZE_XY^2 PIXELSIZE_XY^2 PIXELSIZE_Z^2]));
dmXY=distMat(cordList(:,1:2),diag([PIXELSIZE_XY^2 PIXELSIZE_XY^2]));
dmZ=distMat(cordList(:,3),diag(PIXELSIZE_Z^2));

% jonas,4/09: adjust filterprm for sigmaCorrection. Otherwise, so many 
% spots are being fitted jointly in mammalian metaphase that there is no 
% hope of ever finishing.
if isfield(dataProperties,'sigmaCorrection') && ~isempty(dataProperties.sigmaCorrection)
    FILTERPRM(1:3) = FILTERPRM(1:3)./dataProperties.sigmaCorrection([1,1,2]);
end
spotsidx=rec_find(dmXY,dmZ,1,PIXELSIZE_XY, PIXELSIZE_Z, FILTERPRM);

if nargout > 1
    % make mask for intensity-fit: An ellipsoid with radius 5*sigma. This
    % replaces Dom's old (and slightly wrong) fillSphereData
    ellipsoidRadius = 5*FILTERPRM(1:3);
    maskSize = ceil(ellipsoidRadius);
    % for every pixel: calculate the distance from the origin as a function of
    % the radius of the ellipsoid
    [xx,yy,zz] = ndgrid(-maskSize(1):maskSize(1),...
        -maskSize(2):maskSize(2),...
        -maskSize(3):maskSize(3));
    distance = xx.^2/ellipsoidRadius(1)^2 + ...
        yy.^2/ellipsoidRadius(2)^2 + ...
        zz.^2/ellipsoidRadius(3)^2;
    % every pixel whose center is less than 1 radius from the origin will be
    % counted
    ellipsoid = distance < 1;



    % %create mask 
    % jonas@march10,2009 : use logical mask to save memory
    mask=false(mSize);
    % % default sphere
    % ellipsoid=fillSphereData(5*FT_SIGMA(1),1);

    for i = spotsidx
        cen=round(cordList(i,:));
        %      center=[cen(2) cen(1) cen(3)];
        mask=mask | pushstamp3d(mask,ellipsoid,cen);
    end;

end

function spidx=rec_find(dMatrixXY,dMatrixZ,curIdxList,PIXELSIZE_XY, PIXELSIZE_Z, FILTERPRM)
% recursive search for points which are separated by a distance smaller than MIN_DIST

% use 7 sigma as minimum distance cutoff 
MIN_DIST_XY=7*FILTERPRM(1)*PIXELSIZE_XY;
MIN_DIST_Z=7*FILTERPRM(3)*PIXELSIZE_Z;
%if smaller than d add spot
spidx=curIdxList;
for i=1:size(dMatrixXY,2)
    if (i~=curIdxList(end) && dMatrixXY(curIdxList(end),i)<MIN_DIST_XY && dMatrixZ(curIdxList(end),i)<MIN_DIST_Z)
        if (isempty(spidx) || ~any(spidx==i))
            spidx(end+1)=i; %#ok<AGROW>
            spidx=rec_find(dMatrixXY,dMatrixZ,spidx,PIXELSIZE_XY, PIXELSIZE_Z, FILTERPRM);
        end;
    end;
end;
