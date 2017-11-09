function maskProp = analyze3DMaskGeometry(maskIn,smoothIter,roiInf,avgRad)
%ANALYZE3DMASKGEOMETRY calculates varios properties of the geometry of the objects in the input 3D mask
% 
% maskProp = analyze3DMaskGeometry(maskIn)
% maskProp = analyze3DMaskGeometry(maskIn,smoothIter)
% 
% This function analyzes the geometry of the input 3D mask and returns the
% calculated properties for each object in the mask. This is analagous to
% regionprops.m, but with additional properties which regionprops.m does
% not support for 3D matrices. See output description for list of
% calculated properties.
% 
%
% Input:
% 
%   maskIn - A 3D logical matrix containing the mask to analyze.
%
%   smoothIter - A positive integer scalar specifying the number of
%   iterations of mesh smoothing to perform on the mask surface prior to
%   surface curvature calculations. Due to the discretized nature of the
%   binary mask, the curvature information in the raw mask surface is
%   nearly meaningless, so smoothing is performed beforehand. If set to 0,
%   no smoothing is performed. Default is 3.
%   *NOTE* Insufficient smoothing may result in inaccurate surface
%   curvature values, or in undefined (NaN) surface curvature values.
%
%
% Output:
%
%   maskProp - A 1xM structure, where M is the number of objects in the
%   input mask, with fields containing the calculated properties. The
%   fields include:
%
%       maskProp.Volume - The number of voxels each object occupies.
%       maskProp.Centroid = The centroid of each object.
%       maskProp.SmoothedSurface = The polygon giving the smoothed surface of
%           the mask. Stored in the FV form used by patch, isosurface etc.
%       maskProp.SurfaceNorms = The normals of the smoothed surface at each
%           vertex.
%       maskProp.GaussianCurvature = The gaussian curvature at each face of the
%           smoothed surface.
%       maskProp.MeanCurvature = The mean curvature at each face of the
%           smoothed surface.
%       maskProp.curvaturePC1 = The first principle component of the surface
%           curvature at each face of the smoothed surface.
%       maskProp.curvaturePC2 = The second principle component of the surface
%           curvature at each face of the smoothed surface.
%       maskProp.PixelList = The list of pixels belonging to each object.
%       maskProp.CenterMost = The location of the pixel(s) within each object
%           which is furthest from the object edge. Depending on the object
%           geometry, multiple pixels may share this maximum distance.
%       maskProp.CenterMostDist = The distance of the centermost point(s)
%           from the object edge.
%       
%
% Hunter Elliott
% 3/2011
%
%
%% ----------------- Parameters --------------- %%

smoothMethod = 1;%1 = gaussian smoothing of mask, 2 = curvature smoothing of mesh
nObj = 1; %TEMP - only 1 object now!!!

%% -------------------------- Input ----------------------------- %%

if nargin < 1 || ~islogical(maskIn) || ndims(maskIn) ~= 3
    error('The first input must be a 3-dimensional logical matrix!')
end

if nargin < 2 || isempty(smoothIter)
    smoothIter = .25;%We use this for the isovalue when gaussian smoothing is used.
end

if nargin < 2 || isempty(roiInf)
    isROI = false;
else    
    isROI = true;
end

if nargin < 3 || isempty(avgRad)
    avgRad = 2;
end

%% --------------------------- Init ---------------------------- &&

%Label the mask so we can calculate object-specific properties
maskCC = bwconncomp(maskIn);
labelMask = labelmatrix(maskCC);
[~,objInd] = sort(cellfun(@numel,maskCC.PixelIdxList),'descend');%Sort the objects by size.


%Initialize mask property structure
maskProp(1:nObj) = struct('SmoothedSurface',[],...
                          'SurfaceNorms',[],...
                          'GaussianCurvature',[],...
                          'CurvaturePC1',[],...
                          'CurvaturePC2',[],...
                          'PixelList',[],...
                          'Centroid',[],...
                          'Volume',[],...
                          'CenterMostDist',[],...
                          'CenterMost',[]);

[M,N,P] = size(maskIn);                      

%% -------------- Property Calculation ------------------------- %%


%Loop through the objects and get their properties
for iObj = 1:nObj
    
    
    if smoothMethod == 1
        
        %Smooth the mask
        maskSmooth = fastGauss3D(labelMask == objInd(iObj),1,[5 5 5],2);
        
        
        if isROI
            %If this is an ROI, we do it this way so the cooridantes are still
            %correct in the full image, and so the surface will be properly
            %truncated by the boundaries of the ROI.
            [X,Y,Z] = meshgrid(roiInf.cropY(1):roiInf.cropY(2),roiInf.cropX(1):roiInf.cropX(2),roiInf.cropZ(1):roiInf.cropZ(2));
            maskProp(iObj).SmoothedSurface = isosurface(X,Y,Z,maskSmooth(roiInf.cropX(1):roiInf.cropX(2),roiInf.cropY(1):roiInf.cropY(2),roiInf.cropZ(1):roiInf.cropZ(2)),smoothIter);
            maskProp(iObj).SurfaceNorms = isonormals(X,Y,Z,maskSmooth(roiInf.cropX(1):roiInf.cropX(2),roiInf.cropY(1):roiInf.cropY(2),roiInf.cropZ(1):roiInf.cropZ(2)),maskProp(iObj).SmoothedSurface.vertices);
        else
            %Ge the isosurface and normals                
            maskProp(iObj).SmoothedSurface = isosurface(maskSmooth,smoothIter);
            maskProp(iObj).SurfaceNorms = isonormals(maskSmooth,maskProp(iObj).SmoothedSurface.vertices);        
            
        end
        
        
    else
        
        %Get the mask surface for this object
        maskProp(iObj).SmoothedSurface = isosurface(labelMask == objInd(iObj),0);

        if smoothIter > 0    
            %Smooth the mask surface
            maskProp(iObj).SmoothedSurface = ...
                smoothpatch(maskProp(iObj).SmoothedSurface,0,smoothIter);    
        end    

        %Get the vertex normals of this surface for curvature calculation
        [~,maskProp(iObj).SurfaceNorms] = surfaceNormals(maskProp(iObj).SmoothedSurface);

    end
        
    %Calculate local gaussian and mean curvature of surface
    [maskProp(iObj).GaussianCurvature,maskProp(iObj).MeanCurvature] =...
        surfaceCurvature(maskProp(iObj).SmoothedSurface,maskProp(iObj).SurfaceNorms);

    %Calculate principal curvatures at each face:
    maskProp(iObj).CurvaturePC1 = maskProp(iObj).MeanCurvature + ...
        sqrt(maskProp(iObj).MeanCurvature .^2 - maskProp(iObj).GaussianCurvature);
    maskProp(iObj).CurvaturePC2 = maskProp(iObj).MeanCurvature - ...
        sqrt(maskProp(iObj).MeanCurvature .^2 - maskProp(iObj).GaussianCurvature);
    
    %Calculate locally-averaged curvature data
    maskProp(iObj).locAvgCurv = calcLocalAvgCurvatures(maskProp,avgRad,1,1);

    %Get distance transform of interior of this object
    distX = bwdist(~(labelMask==objInd(iObj)));

    %Get x-y-z coordinates of each point in this object.
    objCoord = [];
    [objCoord(:,2),objCoord(:,1),objCoord(:,3)] = ...
        ind2sub(size(maskIn),maskCC.PixelIdxList{objInd(iObj)});
    
    maskProp(iObj).PixelList = maskCC.PixelIdxList{objInd(iObj)};
    maskProp(iObj).Centroid(iObj,:) = mean(objCoord,1);
    maskProp(iObj).Volume(iObj) = numel(maskCC.PixelIdxList{objInd(iObj)});
    
    %Find centermost point of this object
    currDistX = zeros(size(maskIn));
    currDistX(maskCC.PixelIdxList{objInd(iObj)}) = distX(maskCC.PixelIdxList{objInd(iObj)});    
    maskProp(iObj).CenterMostDist(iObj) = max(currDistX(:));
    %Get x-y-z coordinates of this point
    [maskProp(iObj).CenterMost(:,2),maskProp(iObj).CenterMost(:,1),...
        maskProp(iObj).CenterMost(:,3)] = ind2sub(size(maskIn),...
        find(currDistX == maskProp(iObj).CenterMostDist(iObj))); %In case this is a non-unique point, get all of them
                 
    
end

