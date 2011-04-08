function maskProp = analyze3DMaskGeometry(maskIn,smoothIter)
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




%% -------------------------- Input ----------------------------- %%

if nargin < 1 || ~islogical(maskIn) || ndims(maskIn) ~= 3
    error('The first input must be a 3-dimensional logical matrix!')
end

if nargin < 2 || isempty(smoothIter)
    smoothIter = 3;
elseif numel(smoothIter)> 1 || round(abs(smoothIter)) ~= smoothIter
    error('The smoothIter parameter must be a postive integer scalar!')
end


%% --------------------------- Init ---------------------------- &&

%Label the mask so we can calculate object-specific properties
maskCC = bwconncomp(maskIn);
labelMask = labelmatrix(maskCC);
nObj = maskCC.NumObjects;


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


%% -------------- Property Calculation ------------------------- %%


%Loop through the objects and get their properties
for iObj = 1:nObj
    
    %Get the mask surface for this object
    maskProp(iObj).SmoothedSurface = isosurface(labelMask == iObj,0);
    
    if smoothIter > 0    
        %Smooth the mask surface
        maskProp(iObj).SmoothedSurface = ...
            smoothpatch(maskProp(iObj).SmoothedSurface,0,smoothIter);    
    end    

    %Get the vertex normals of this surface for curvature calculation
    [~,maskProp(iObj).SurfaceNorms] = surfaceNormals(maskProp(iObj).SmoothedSurface);

    %Calculate local gaussian and mean curvature of surface
    [maskProp(iObj).GaussianCurvature,maskProp(iObj).MeanCurvature] =...
        surfaceCurvature(maskProp(iObj).SmoothedSurface,maskProp(iObj).SurfaceNorms);

    %Calculate principal curvatures at each face:
    maskProp(iObj).CurvaturePC1 = maskProp(iObj).MeanCurvature + ...
        sqrt(maskProp(iObj).MeanCurvature .^2 - maskProp(iObj).GaussianCurvature);
    maskProp(iObj).CurvaturePC2 = maskProp(iObj).MeanCurvature - ...
        sqrt(maskProp(iObj).MeanCurvature .^2 - maskProp(iObj).GaussianCurvature);

    %Get distance transform of interior of this object
    distX = bwdist(~(labelMask==iObj));

    %Get x-y-z coordinates of each point in this object.
    objCoord = [];
    [objCoord(:,2),objCoord(:,1),objCoord(:,3)] = ...
        ind2sub(size(maskIn),maskCC.PixelIdxList{iObj});
    
    maskProp(iObj).PixelList = maskCC.PixelIdxList{iObj};
    maskProp(iObj).Centroid(iObj,:) = mean(objCoord,1);
    maskProp(iObj).Volume(iObj) = numel(maskCC.PixelIdxList{iObj});
    
    %Find centermost point of this object
    currDistX = zeros(size(maskIn));
    currDistX(maskCC.PixelIdxList{iObj}) = distX(maskCC.PixelIdxList{iObj});    
    maskProp(iObj).CenterMostDist(iObj) = max(currDistX(:));
    %Get x-y-z coordinates of this point
    [maskProp(iObj).CenterMost(:,2),maskProp(iObj).CenterMost(:,1),...
        maskProp(iObj).CenterMost(:,3)] = ind2sub(size(maskIn),...
        find(currDistX == maskProp(iObj).CenterMostDist(iObj))); %In case this is a non-unique point, get all of them
                 
    
end

