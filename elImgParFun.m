function val = elImgParFun(varargin)
%elImgParFun : Define any parameter in the elastic equation based on
%              the intensity of an image.
%
% SYNOPSIS : 
%    val = elImgParFun(x,y,imI,A,V)
%    This function can be used to define any parameter in the elastic equation
%    whose values equal the base value 'A' plus the variation 'V' multiplied
%    by the intensity image 'imI'. The intensity of the image is scaled to be
%    between 0 and 1.
%
%    val = elImgParFun(x,y,imI,A,V,threshold)
%    A threshold in term of the percentage of the maximum
%    intensity can also be specified so that any intensity that
%    is below the threshold will be set to the base value 'A'. Pass [] to
%    ignore it.
%
%    val = elImgParFun(x,y,imI,A,V,threshold,method,par1,par2,...,parm)
%    The image can also be filtered with the following methods:
%
%    method : A string that specifies the interpolation method to use.
%       'spap2'    : Least square B-Spline.
%       'spaps'    : Smoothing spline.
%       'Gaussian' : Convolution with Gausian kernel.
%       'ppform'   : The pp- or B-form of the spline interpolation of the image
%                   have already been given.
%    par1,...,parm : Parameters that are needed for the various interpolation
%       methods.
%        method        (par1,...,parm)
%       ------------------------------------
%       'spap2'    : ({knots1,knots2}, [k1 k2])
%       'spaps;    : (tol)
%       'Gaussian' : (corLen)
%       'ppform'   : The pp- or B-form of the spline interpolation.
%
%    See help imInterp.
%
% INPUT :
%    x,y : Points where the parameter is evaluated.
%    imI : The intensity matrix of the image. It can be of class double, uint8
%          or uint16.
%    A   : The base value of the parameter.
%    V   : The variation of the parameter.

if nargin < 5
   error('Not enough input arguments.');
end

x   = varargin{1};
y   = varargin{2};
imI = double(varargin{3});
A   = varargin{4};
V   = varargin{5};

if V == 0
   val = A*ones(size(x));
   return;
end

xySZ = size(x);

x = reshape(x,length(x(:)),1);
y = reshape(y,length(y(:)),1);

if nargin > 6
   if strcmp(varargin{7},'ppform') == 1
      val = fnval(varargin{7},[y x].');
   else
      val = imInterp(imI,[y x],varargin{7:end});
   end
else
   val = imInterp(imI,[y x]);
end

threshold = 0;
if nargin >= 6 & ~isempty(varargin{6})
   threshold = varargin{6};
end
val = (val-threshold)*V;
val(find(val<0)) = 0;

val = A + reshape(val,xySZ);
