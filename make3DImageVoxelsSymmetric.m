function image = make3DImageVoxelsSymmetric(image,pixXY,pixZ)
%MAKE3DIMAGEVOXELSSYMMETRIC scales the input 3D image so that the voxels are symmetric 
% 
% image = make3DImageVoxelsSymmetric(image,pixXY,pixZ)
% 
% This function scales the input 3D image matrix so that the resulting
% voxels are symmetric. By symmetric it is meant that the physical
% dimensions of the voxels in each dimension. This is done based on the
% input pixel sizes in X-Y and Z. It is assumed that the last dimension of
% the matrix is the Z axis, and that the pixel sizes in the X and Y
% directions are equal. The scaling is done in the Z-direction, so that the
% resulting matrix is the same size in the first two dimensions, but larger
% in the third dimension.
%   This function is intended to be used to pre-process 3D image data where
% the spacing between z-planes is not equal to the pixel size in X-Y. The
% units of the input pixel sizes are unimportant, as only the relative
% sizes are used in scaling.
%
% Input:
% 
%   image - 3D image matrix.
% 
%   pixXY - Positive scalar specifying the size of the pixels in the X-Y
%   plane (dimensions 1,2)
% 
%   pixZ - Positive scalar specifying the size of the pixels in the z
%   direction (dimension 3)
% 
% 
% Output:
%
%   image - 3D matrix containing symmetric voxels, whose sizes in the x,y
%   and z directions are all equal to the input pixel size in the x-y
%   direction, pixXY.
%
% Hunter Elliott
% 4/2011
%

if nargin < 3 || isempty(image) || isempty(pixXY) || isempty(pixZ)
    error('You must input an image, and x-y pixel size and a z-pixel size!')
end

if ndims(image) ~=3
    error('Input image must be 3-dimensional!')
end

if any([pixXY pixZ]) <= 0 || numel(pixXY)>1 || numel(pixZ)>1
    error('The input pixel sizes must be positive scalars!');
end

%Get the image class
ogClass = class(image);
%Always do interpolation as double. interp3 doesnt support uint classes.
image = double(image);


%Factor to scale z dimension by.
scFact = pixZ/pixXY;


%Set up the interpolation points
[M,N,P] = size(image);
[X,Y,Z] = meshgrid(1:N,1:M,1:P);
[Xi,Yi,Zi] = meshgrid(1:N,1:M,linspace(1,P,P*scFact));

%Interpolate the image.
image = interp3(X,Y,Z,image,Xi,Yi,Zi);

if ~strcmp(ogClass,'logical')    
    %Restore the original class if it has changed
    if ~strcmp(ogClass,'double')
        image = cast(image,ogClass);
    end
else
    %Special case for logical - a cast to logical is equivalent to
    %ceil(image), whereas we want round(image)
    image = image >= .5;
    
end
