function [lMax]=validateSpeckles(img,lMax,info,parameters)

% SYNOPSIS   [lMax]=validateSpeckles(img,lMax,info,parameters)
%
% INPUT      img        : image
%            lMax       : local max map of image
%            info       : Imax and Imin for every speckles obtained 
%                         from analyzeSpeckles (Delaunay triangulation)
%            parameters : [k sDN sP]
%                             k: parameter
%                             sDN : dark noise
%                             sP  : Poisson noise (dependent on I)%                              
% OUTPUT     lMax       : the loc max map with only those loc max which have been 
%                         selected as speckles
%

% Parameters
k=parameters(1);
sigmaD=parameters(2);
PoissonNoise=parameters(3);

% Correct intensity at the four edges 
img([1 size(img,1) (size(img,2)-1)*size(img,1)+1 size(img,1)*size(img,2)])=mean(img(:));

% Create addressing matrix
A=zeros(size(img));

% Calculate difference Imax-Imin for every speckle
for c1=1:size(info,3)
   
   % Read values for Imax and Imin from the image at the coordinates specified by info
   Imax=img(info(1,1,c1),info(1,2,c1));
   Imin=mean([img(info(2,1,c1),info(2,2,c1)) img(info(3,1,c1),info(3,2,c1)) img(info(4,1,c1),info(4,2,c1))]);
   
   % Calculate error on Imax (noise model)
   sigmaMax=sigmaD+sqrt(PoissonNoise*Imax);
   
   % Calculate error on Imin (noise model)
   sigmaMin=sigmaD+sqrt(PoissonNoise*Imin);
   
   % Calculate difference and error
   deltaI=Imax-Imin;
   sigmaDiff=sqrt(sigmaMax^2+sigmaMin^2);
   
   % Check for the validity of the speckle
   if deltaI>=k*sigmaDiff;
      A(info(1,1,c1),info(1,2,c1))=1;   % Loc max accepted as speckle
   end
end

% Removing false speckles from loc max map
lMax=A.*lMax;
