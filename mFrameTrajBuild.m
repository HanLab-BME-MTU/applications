function MFT = mFrameTrajBuild(varargin)
%mFrameTrajBuild Create trajectories of speckles over multiframes.
%
% SYNOPSIS : 
%    MFT = mFrameSpeckTrack(M,corLen,YX)
%    MFT = mFrameSpeckTrack(M,...)
%
% INPUT :
%    M : The vector field of the displacements of speckles created by Aaron's 
%       Speckle Tracking Software. It is an m-by-4-by-n multi-dimentional array
%       where 'm' is roughly (due to the death of old speckles and the birth of 
%       new speckles) the number of speckles and 'n+1' is the total number of 
%       frames stored in 'M'. The 4 columns have the form [y0 x0 y x] where 
%       (y0,x0) is the base and (y,x) is the tip of the displacement vector.
%    corLen : The correlation length used for the interpolation of the vector
%       field by convolution. If it is not given, the default value is 2
%       (pixels).
%    YX : The coordinates of the points whose trajectories are to be tracked.
%       It is an m-by-2 matrix in the form [y x] where 'm' is the number of 
%       points. If it is not given, the speckles in the first frame of 'M' are
%       tracked.
%
% OUTPUT :
%    MFT : An m-by-2*(n+1) matrix where 'm' is the number of points to be
%       tracked and 'n+1=size(M,3)+1' is the total number of frames stored in 
%       'M'. The '2*(n+1)' columns have the form [y1 x1 ... yi xi ...]
%       where (yi,xi) is the tracked position of the points in the i-th 
%       frame.

%Parse and check the inputs.
[M,corLen,YX] = parse_inputs(varargin{:});

%Get the total number of frames in 'M'.
numFrames = size(M,3)+1;

%Number of points to be tracked.
numPoints = size(YX,1);

MFT = zeros(numPoints,2*numFrames);

frameID    = 1;
MFT(:,1:2) = YX;
for jj = 1:2:2*numFrames-3
   MFT(:,jj:jj+3) = vectorFieldInterp( ...
      M(find(M(:,1,frameID)~=0 & M(:,3,frameID)~=0),:,frameID), ...
      MFT(:,jj:jj+1),corLen,[]);
   frameID = frameID+1;
end

%%%%% End of the Main Function %%%%%%%%%%%%%%%%%%

function [M,corLen,YX] = parse_inputs(varargin)
%Parse and check the inputs.

if nargin < 1
   error('Not enough input arguments.');
elseif nargin > 3
   error('Too many input arguments.');
end

M = varargin{1};

if nargin == 3
   YX = varargin{3};
else
   YX = M(find(M(:,1,1)~=0 & M(:,3,1)~=0),1:2,1);
end

if nargin >= 2
   corLen = varargin{2};
else
   corLen = 2;
end

%Check the legitimacy of the input arguments.
if ~isnumeric(corLen) | length(corLen) ~= 1 | ...
   corLen <= 0
   error('The correlation length is a positive numerical value.');
end

if ~isnumeric(YX) | ndims(YX) > 2 | size(YX) ~= 2
   error(['The coordinates of the points to be tracked is not ' ...
      'correctly defined.']);
end
