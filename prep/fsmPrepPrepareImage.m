function img=fsmPrepPrepareImage(img,IcorrFactor,W)
% fsmPrepPrepareImage prepares the image to be subjected to further analysis in LTA
%
% Images in fsm are loaded and prepared in fsmPrepMain.m and in fsmBuildSaveSpeckleArray.m.
% For this reason, this function has been written in order to assure that no changes
% in the images are introduced by wrong handling
%
% SYNOPSIS      img=fsmPrepPrepareImage(img,IcorrFactor,W)
%
% INPUT         img        : loaded input image
%               IcorrFactor: factor to compensate the overall intensity variation
%                            due to bleaching / focusing effects
%               W          : matrix defining the grid overlaid to the image
%                            (see function FRAMEWORK for more details)
%
% OUTPUT        img        : prepared image
%
% DEPENDENCES   fsmPrepPrepareImage uses { }
%               fsmPrepPrepareImage is used by { fsmPrepMain, fsmBuildSaveSpeckleArray } 
%
% Aaron Ponti, October 4th, 2002

% Check for image class
if ~isa(img,'double')
    img=double(img);
end

% Multiply image with corresponding factor
img=IcorrFactor*img;

% Cut images (defined by grid)
img=img(W(1,1):W(end,3),W(1,2):W(end,4));

% Filter with correction of border effect
bck=mean(img(:))*ones(size(img));
img=Gauss2D(img,1);
bck(3:end-2,3:end-2)=img(3:end-2,3:end-2); 
img=bck;
