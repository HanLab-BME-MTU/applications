function M = vectorFieldAnimate(varargin)
%vectorFieldAnimate  Show the animation of the given vector field over time.
%
% SYNOPSIS :
%    M = vectorFieldAnimate(XY,V,scale,'name1',value1,...)
%    M = vectorFieldAnimate(XY,V,...)
%
% INPUT :
%    XY : Points where the vector field is shown. It is an m-by-2 matrix in
%       the form [x y] where 'm' is the number of points.
%    V : The stack of vector fields over time. It is an m-by-2-by-n
%       multidimentional array where 'm' equals the number of points stored in
%       'XY' and 'n' is the number of time steps (or frames). The two colums
%       of the 2nd dimention has the form [vx vy] where 'vx' and 'vy' are the
%       x and y components of the vector respectively.
%    scale : A non-negative numerical value. When it is zero, the vectors
%       will be automatically scaled. The scaling is the same for all the
%       frames in order to display the dynamical change of the vector length.
%       The default value is zero (automatical scaling).
%
%    The following properties can also be specified as optional inputs.
%    'bgImg' : It can be either one image or a cell array of
%       images whose length equals the number of time steps.
%    'vc' : The color of the vector. Default: 'b' (blue).
%
% OUTPUT :
%    M : A matlab movie is returned.

%Parse the inputs.
if nargin < 2
   error('Not enought input arguments.');
end

XY = varargin{1};
V  = varargin{2};

%Check the legitimacy of 'XY' and 'V'
if ~isnumeric(XY) | ndims(XY) > 2 | size(XY,2) ~= 2
   error(['The first argument should be an m-by-2 numerical matrix that ' ...
      'specifies the coordinates of the points where the vector field ' ...
      'is shown.']);
end

if ~isnumeric(V) | ndims(V) > 3 | size(V,2) ~= 2 | ...
   size(V,1) ~= size(XY,1)
   error(['The second argument that specifies the vector field is not ' ...
      'defined correctly.']);
end

numFrames = size(V,3);

scale  = 0; %Default scale.
pStart = 3; %Where the properties start.
if nargin >= 3 & isnumeric(varargin{3})
   scale = varargin{3};

   if length(scale) ~= 1
      error(['The scale specified in the third argument is ' ...
         'a single numerical value.']);
   end
   if scale < 0
      error(['The scale specified in the third argument can not ' ...
         'be negative.']);
   end
   pStart = 4;
end

%Default property values.
bgImg = [];
vc    = 'b';

if nargin >= pStart
   if rem(nargin-pStart+1,2) ~= 0
      error('The properties and values specified do not match in pair.');
   end

   numProperties = (nargin-pStart+1)/2;
   inProperty    = cell(numProperties,1);
   inValues      = cell(numProperties,1);
   for k = 1:numProperties
      inProperty{k} = varargin{pStart+2*k-2};
      inValue{k}    = varargin{pStart+2*k-1};

      switch inProperty{k}
      case 'bgImg'
         if iscell(inValue{k})
            if length(inValue{k}) ~= 1 | ...
               length(inValue{k}) ~= numFrames
               error(['The number of background images can either be ' ...
                  'one or the number of frames.']);
            end

            if ischar(inValue{k}{1})
               for jj = 1:length(inValue{k})
                  bgImg{jj} = imread(inValue{k}{jj});
               end
            else
               bgImg = inValue{k};
            end
         else
            if ischar(inValue{k})
               bgImg{1} = imread(inValue{k});
            else
               bgImg{1} = inValue{k};
            end
         end
      case 'vc'
         vc = inValue{k};
      otherwise
         error([inProperty{k} ' is not an acceptable property.']);
      end
   end
end

%Display the vector field and make the movie.
figure(gcf); hold off;
for k = 1:numFrames
   if ~isempty(bgImg)
      if length(bgImg) == 1
         imshow(bgImg{1},[]); hold on;
      else
         imshow(bgImg{k},[]); hold on;
      end
   end

   quiver(XY(:,1),XY(:,2),V(:,1,k)*scale,V(:,2,k)*scale,0,vc); hold off;
   M(k) = getframe;
end
