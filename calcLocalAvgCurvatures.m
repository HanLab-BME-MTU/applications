function curvOut = calcLocalAvgCurvatures(maskProp,radius,dsSurf,dsSamp)

%calculates local averages of the various surface properties within the
%specified radius.
%
%

%Hunter Elliott, 1/2013

nFaces = size(maskProp.SmoothedSurface.faces,1);

%Get the face center positions, since we have curv data at the per-face
%level.

if dsSurf > 1
    useFaces = randperm(nFaces,round(nFaces/dsSurf));
else
    useFaces = 1:nFaces;
end
nFaceUse = numel(useFaces);
faceCenters = zeros(nFaceUse,3);

for j = 1:nFaceUse
    faceCenters(j,:) = mean(maskProp.SmoothedSurface.vertices(maskProp.SmoothedSurface.faces(useFaces(j),:),:),1);
end

%Find points within the averaging radius of each surface face

iClosest = KDTreeBallQuery(faceCenters,faceCenters(1:dsSamp:end,:),radius);


%Now average all these.
maxAbsPC = max(abs(real([maskProp.CurvaturePC1 maskProp.CurvaturePC2])),[],2);
for j = 1:numel(iClosest)
    curvOut.LocMeanCurvaturePC1(j,1) = mean(maskProp.CurvaturePC1(iClosest{j}));
    curvOut.LocMeanCurvaturePC2(j,1) = mean(maskProp.CurvaturePC2(iClosest{j}));
    curvOut.LocMeanMeanCurvature(j,1) = mean(maskProp.MeanCurvature(iClosest{j}));
    curvOut.LocMeanGaussianCurvature(j,1) = mean(maskProp.GaussianCurvature(iClosest{j}));
    curvOut.LocMeanMaxAbsCurvature(j,1) = mean(maxAbsPC(iClosest{j}));
end



