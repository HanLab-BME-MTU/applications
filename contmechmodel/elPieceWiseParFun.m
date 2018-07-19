function val = elPieceWiseParFun(x,y,A,roi,V,varargin)
%elPieceWiseParFun : Define any parameter in the elastic equation by assigning 
%                    values to a set of subregions.
%
% SYNOPSIS : 
%    val = elImgParFun(x,y,A,roi,V)
%    This function can be used to define any parameter in the elastic equation
%    which has different values on a set of subregions defined by 'roi'.
%
% INPUT :
%    x,y : Points where the parameter is evaluated.
%    A   : The base value for the background region.
%    roi : A cell array of polygons that define the set of subregions. Each 
%          element of the cell array is a two-column vector [x y] that specifies
%          the vertices of the polygon.
%    V   : A vector that specify the values for the set of subregions.

if nargin < 5
   error('Not enough input arguments.');
end

numSubRegions = length(roi);

if V == 0
   val = A*ones(size(x));
   return;
end

xySZ = size(x);

x = reshape(x,length(x(:)),1);
y = reshape(y,length(y(:)),1);

val = A*ones(size(x));
for k = 1:numSubRegions
   in = inpolyon(x,y,roi{k}(:,1),roi{k}(:,2));
   val(find(in==1)) = V(k);
end

val = reshape(val,xySZ);
