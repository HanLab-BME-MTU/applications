function [Inew,Imaxima,nmb]=fsmPrepSubstructMaxima(IG,Imax,SIG,cands)

% fsmPrepSubstructMaxima takes the filtered image IG and substructs Gaussian Kernels 
% at the positions of the local maxima
%
%
% SYNOPSIS   Inew=fsmPrepSubstructMaxima(IG,Imax,SIG,cands)
%
% INPUT      IG         :  filtered image
%            Imax       :  local maxima map
%            SIG        :  sigma of the GKs
%            cands      :  cands strucure
%
% OUTPUT     Inew       :  the resulting image after substaction
%            Imaxima    :  the synthetic image to be substracted from the real
%            nmb        :  number of gaussians kernels/speckles to be substracted
%        
%
% DEPENDENCES   fsmPrepSubstructMaxima uses { Gauss2d1 } 
%               fsmPrepSubstructMaxima is used by { fsmPrepMainSecondarySpeckles }
%
% Alexandre Matov, November 7th, 2002

Imax=zeros(size(Imax)); % prepearing Imax for substruction 

% setting the value of Imax to the values from delataI, i.e. taking into account the background
s=length(cands);

% counter
nmb=0;

for i=1:s
    if cands(i).status==1 & cands(i).deltaI>0 % status flag - a local maximun is significant or not
        Imax(cands(i).Lmax(1),cands(i).Lmax(2))=cands(i).deltaI; % the intensity at the speckle position is set to "delta I"
        nmb=nmb+1;
    end

end

% the masked with a GK local maxima points with intensity delta.I for the raw data image
Imaxima=Gauss2D1(Imax,SIG);

% substruction  
Inew=IG-Imaxima;