function maskProp = analyze3DMaskGeometry(maskIn,varargin)
%ANALYZE3DMASKGEOMETRY calculates varios properties of the geometry of the objects in the input 3D mask
% 
% maskProp = analyze3DMaskGeometry(maskIn)
% maskProp = analyze3DMaskGeometry(maskIn,'OptionName1',optionValue1,'OptionName2',optionValue2,...)
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
%   'OptionName',optionValue - A string with an option name followed by the
%   value for that option. 
%
%     Possible option name/value pairs:
%
%       ('OptionName'->possible values)
%   
%       ('SmoothSigma'->positive scalar) The sigma of the gaussian filter
%       to use to smooth the mask prior calculating surface curvature
%       properties. If set to 0, no smoothing is done.
%       Default is 1.
%
%       ('IsoValue'->positive scalar 0 > x < 1) Specifies the value to
%       threshold the smoothed mask image at to obtain the smoothed mask.
%       Default is .5
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

[smoothSig,isoVal] = parseInput(varargin);

if isempty(smoothSig)
    smoothSig = 1;
end

if isempty(isoVal)
    isoVal = .5;
elseif isoVal >= 1 || isoVal <= 0
    error('The IsoValue input must be greater than 0 and less than 1!')
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
    
    %Smooth each object independently to avoid object merging via
    %smoothing.    
    if smoothSig > 0    
        %Smooth the shit out of the mask
        smMask = fastGauss3D(double(labelMask==iObj),smoothSig);
    else
        smMask = double(labelMask==iObj);
    end    

    %Extract isosurface of this object
    maskProp(iObj).SmoothedSurface = isosurface(smMask,isoVal);
    maskProp(iObj).SurfaceNorms = isonormals(smMask,maskProp(iObj).SmoothedSurface.vertices);

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

function [smoothSig,isoVal] = parseInput(argArray)
%Sub-function for parsing variable input arguments


smoothSig = [];
isoVal = [];


if isempty(argArray)
    return
end

nArg = length(argArray);


%Make sure there is an even number of arguments corresponding to
%optionName/value pairs
if mod(nArg,2) ~= 0
    error('Inputs must be as optionName/ value pairs!')
end

for i = 1:2:nArg
    
   switch argArray{i}                     
                         
       case 'SmoothSigma'
           
           smoothSig = argArray{i+1};
           
       case 'IsoValue'
           
           isoVal = argArray{i+1};
           
       otherwise
                  
           imviewArgs  = [imviewArgs argArray{i:i+1}]; %#ok<AGROW>
           
   end
               
      
   
end
