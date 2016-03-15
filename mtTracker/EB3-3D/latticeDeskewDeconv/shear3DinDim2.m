function [ output ] = shear3DinDim2(inArray, angle, bReverse, dz, xypixelsize, trans, fillVal)
%   shear3DinDim2 -- deskew inArray about the 1st dimension (y in image terms)
%   x in sheet-scan terms)
%   inArray -- input 3D array
%   angle -- skewing angle
%   bReverse -- is skewing in reverse direction ? (like in Bi-Chang's data)
%   dz -- sample-scan step size (in microns)
%   xypixelsize -- in microns
%   trans -- extra left-right translation (positive number translates the result toward
%   the right)
%   fillVal -- use this value to fill the new void

center = (size(inArray)+1)/2;
trans1 = [1 0 0 0
    0 1 0 0
    0 0 1 0
    -center 1];

rot = [1 0 0 0
    0 1 0 0
    0 cos(angle*pi/180)*dz/xypixelsize 1 0
    0 0 0 1];

if bReverse
    rot(3, 2) = -rot(3, 2);
end

% widenBy = ceil(size(inArray, 2) * xypixelsize /(cos(angle*pi/180)*dz));
widenBy = ceil(size(inArray, 3)*cos(angle*pi/180)*dz/xypixelsize);


trans2 = [1 0 0 0
    0 1 0 0
    0 0 1 0
    center+[0 widenBy/2+trans 0] 1];

T = trans1*rot*trans2;
tform = maketform('affine', T);
R = makeresampler('cubic', 'fill');
output = tformarray(inArray, tform, R, [1 2 3], [1 2 3], size(inArray)+[0 widenBy 0] , [], fillVal);

end
