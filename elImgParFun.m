function val = elImgParFun(x,y,imI,A,V,threshold)
%elImgParFun : Define any parameter in the elastic equation based on
%              the intensity of an image.
%
% SYNOPSIS : val = elImgParFun(x,y,imI,A,V,threshold)
%    This function can be used to define any parameter in the elastic equation
%    whose values equal the base value 'A' plus the variation 'V' multiplied
%    by the intensity image 'imI'. The intensity of the image is scaled to be
%    between 0 and 1.
%
% INPUT :
%    x,y : Points where the parameter is evaluated.
%    imI : The intensity matrix of the image. It can be of class double, uint8
%       or uint16.
%    A : The base value of the parameter.
%    V : The variation of the parameter.
%    threshold : A trhreshold in term of the percentage of the maximum
%       intensity can also be specified so that any intensity that
%       is below the threshold will be set to the base value 'A'. Pass [] to
%       ignore it.

if V == 0
   val = A*ones(size(x));
   return;
end

xySZ = size(x);

x = reshape(x,length(x(:)),1);
y = reshape(y,length(y(:)),1);

[sp,val] = imInterp(imI,[x y],[],'spap2',5);
%val = imInterp(imI,[y x],[],'Gaussian',10);

if ~isempty(threshold)
   val(find(val<max(abs(val(:)))*threshold)) = 0;
end

val = A + V*reshape(val,xySZ);
