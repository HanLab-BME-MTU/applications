function vargout = imInterp(vargin)
%IMINTERP Interpolation of image intensity to subpixel points and return the 
%         intensity as double values.
%
% SYNOPSIS :
%    outI = imInterp(img,XX,method,par1,...,parm,range)
%    The function interpolates the intensity of an image to points specified
%    in 'XX'. The output of the interpolated values are of class double and 
%    are linearly scaled so that the minimum intensity is mapped to 'low' and 
%    the maximum is mapped to 'high' where 'low' and 'high' are specified in 
%    'range = [low high]'. When the input, 'range' is ommitted, the default
%    value is [0 1].
%
%    [sp,outI] = imInterp(img,XX,method,par1,...,parm,range)
%    When 'method' is the B-Spline interpolation, the B-form of the
%    interpolation, can also be output in the variable, 'sp'. 
%
%    sp = imInterp(img,[],...)
%    When no interpolation point is specified, the output is the B-form of the
%    interpolation. In this case, 'method' can only be one of the B-Spline
%    interpolation methods.
%
% INPUT :
%    img : The image intensity matrix. It can be of class uint8, uint16, or
%       double. It can also be a string that specifies the name of the image
%       file.
%    XX  : An m-by-2 matrix that specifies the coordinates of m points where
%       we want to find the values by interpolation. It is given in the form
%       [x y] where 'x' is the x coordinates and 'y' is the y coordinates. The
%       unit of the coordinates are in pixels. For points that are outside of
%       the image, the output is NaN.
%    method : A string that specifies the interpolation method to use.
%       'spap2'   : Least square B-Spline.
%       'spaps'   : Smoothing spline.
%       'Gausian' : Convolution with Gaunsian kernel.
%    par1,...,parm : Parameters that are needed for the various interpolation
%       methods.
%        method        (par1,...,parm)
%       ------------------------------------
%       'spap2'   : ({knots1,knots2}, [k1 k2])
%       'spaps;   : (tol)
%       'Gausian' : (corLen)
%    range : The range of the values to which the image intensity is
%       interpolated (or mapped) to. It is given in the form [low high] where
%       the lowest intensity is mapped to 'low' and the highest intensity is
%       mapped to 'high'. The default value is [0 1].
%
% OUTPUT :
%    outI : The scaled intensities at the points given in 'XX'. The values are
%       of class double even if the input image is of class uint8, or uint16.
%    sp : The B-form of the spline interpolation when 'method' is one of the
%       spline interpolation methods.

%Set the default range of values to which the intensities are mapped.
range = [0 1];

if nargin < 2
   error('There have to be at least two input arguments.');
elseif nargin == 2
   img = vargin{1};
   XX  = vargin{2};

   if isempty(XX)
      error(['The interpolation points have to be specified ' ...
         'when there are only two input arguments.']);
   end

   %The default interpolation method.
   method = 'Gausian';
   corLen = 1;
else
   img    = vargin{1};
   XX     = vargin{2};
   method = vargin{3}
end

%Check the legitimacy of the inputs, 'img' and 'XX'.
if ischar(img)
   imI = imread(img);
else
   imI = img;
end

if isrgb(imI)
   error('RGB images are not supported. Call RGB2GRAY first.');
elseif isind(imI)
   error('Indexed images are not supported. Call IND2GRAY first.');
elseif ~isgray(imI)
   error('Not an image intensity matrix.');
end

numPixelsX = size(imI,2);
numPixelsY = size(imI,1);

%Parse the rest of the arguments according to 'method'.
if strcmp(method,'Gaussian') == 1
   corLen = 1; %Default correlation length.

   if nargin > 5 
      error(['Too many input arguments for the specified ' ...
         'interpolation method.']);
   elseif nargin == 5 & ~isempty(vargin{5})
      range  = vargin{5};
   end

   if nargin >= 4 & ~isempty(vargin{4})
      corLen = vargin{4};
   end

   %Check if 'corLen' is set correctly.
   if ~isnumeric(corLen) | (isnumeric(corLen) & length(corLen) > 2)
      error(['The correlation length should be provided as a ' ...
         'numerical array of 1 or 2 numbers.']);
   elseif min(corLen) < 0
      error('The correlation length has to be positive.');
   end

   if length(corLen) == 1
      corLenX = corLen;
      corLenY = corLen;
   else
      corLenX = corLen(1);
      corLenY = corLen(2);
   end
elseif strcmp(method,'spap2') == 1
   order = 4; %Default order of spline interpolation.

   if nargin > 6
      error(['Too many input arguments for the specified ' ...
         'interpolation method.']);
   elseif nargin == 6 & ~isempty(vargin{6})
      range = vargin{6};
   end

   if nargin >= 5 & ~isempty(vargin{5})
      order = vargin{5};
   end

   %Check if 'order' is correctly defined.
   if ~isnumeric(order) | (isnumeric(order) & length(order) > 2)
      error(['The order of the spline interpolation should be ' ...
         'provided as a numerical array of 1 or 2 positive integers.']);
   elseif min(order) < 0
      error('The order of the spline interpolation has to be positive.');
   end

   if length(order) == 1
      orderX = order;
      orderY = order;
   else
      orderX = order(1);
      orderY = order(2);
   end

   if floor(orderX) ~= orderX | floor(orderY) ~= orderY
      error('The order of the spline interpolation has to be integer.');
   end

   %If no knot sequence is provided, the default knot sequence with
   % intervals of appoximately every other pixel is created and the default
   % order is 4. For even number of pixels, the end interval has a length
   % of 3 pixels.
   if rem(numPixelsX/2) == 0 
      knotsX = augknt([1:2:numPixelsX-3 numPixelsX],orderX);
   else
      knotsY = augknt([1:2:numPixelsX],orderX);
   end

   if rem(numPixelsY/2) == 0 
      knotsX = augknt([1:2:numPixelsY-3 numPixelsY],orderY);
   else
      knotsY = augknt([1:2:numPixelsY],orderY);
   end

   if nargin >= 4 & ~isempty(vargin{4})
      if iscell(vargin{4}) & length(vargin{4}) == 2
         knotsX = vargin{4}{1};
         knotsY = vargin{4}{2};
      else
         error(['The two knot sequences should be provided ' ...
            'in a cell array of two elements.']);
      end
   end
elseif strcmp(method,'spaps') == 1
else
   error('The specified interpolation method is not recogonized.');
end

%Check if 'range' is correctly defined.
if ~isnumeric(range) | (isnumeric(range) & length(range) > 2) | ...
   (isnumeric(range) & range(1) == range(2)
   error(['The range of values to which the intensity is ' ...
      'mapped to should be provided in a numerical array ' ...
      'two different numbers.']);
end

%Start the interpolation.
if strcmp(method,'Gaussian') == 1
end
