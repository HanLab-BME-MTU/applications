function img=fsmPrepPrepareImage(img,IcorrFactor,W,sigma)
% fsmPrepPrepareImage prepares the image to be subjected to further analysis in LTA
%
% Images in fsm are loaded and prepared in fsmPrepMain.m and in fsmBuildSaveSpeckleArray.m.
% For this reason, this function has been written in order to assure that no changes
% in the images are introduced by wrong handling
%
% SYNOPSIS      img=fsmPrepPrepareImage(img,IcorrFactor,W,sigma)
%
% INPUT         img        : loaded input image
%               IcorrFactor: factor to compensate the overall intensity variation
%                            due to bleaching / focusing effects
%               W          : matrix defining the grid overlaid to the image
%                            (see function FRAMEWORK for more details)
%               sigma      : sigma for image low-pass filtering (optional, default = 1).
%
% OUTPUT        img        : prepared image
%
% DEPENDENCES   fsmPrepPrepareImage uses { }
%               fsmPrepPrepareImage is used by { fsmPrepMain, fsmBuildSaveSpeckleArray } 
%
% Aaron Ponti, October 4th, 2002

% Check input
if nargin<3 | nargin>4
    error('Three or four input parameters expected');
end

% Sigma has a default value of 1.
if nargin==3
    sigma=1;
end

% Check for image class
if ~isa(img,'double')
    img=double(img);
end

% Multiply image with corresponding factor
img=IcorrFactor*img;

% Cut images (defined by grid)
img=img(W(1,1):W(end,3),W(1,2):W(end,4));

% Filter with correction of border effect (if needed)
if sigma~=0
    bck=mean(img(:))*ones(size(img));
    [img,M]=Gauss2D(img,sigma);
    dy=fix(size(M,1)/2);
    dx=fix(size(M,2)/2);
    bck(dy+1:end-dy,dx+1:end-dx)=img(dy+1:end-dy,dx+1:end-dx); 
    img=bck;
end