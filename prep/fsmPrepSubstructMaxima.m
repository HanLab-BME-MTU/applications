function [Inew,Imaxima]=fsmPrepSubstructMaxima(IG,SIG,cands)

% fsmPrepSubstructMaxima takes the filtered image IG and substracts Gaussian
% kernels at the positions of the local maxima
%
%
% SYNOPSIS   Inew=fsmPrepSubstructMaxima(IG,Imax,SIG,cands)
%
% INPUT      IG         :  filtered image
%            SIG        :  sigma of the GKs
%            cands      :  cands strucure
%
% OUTPUT     Inew       :  the resulting image after substaction
%            Imaxima    :  the synthetic image to be substracted from the
%                          real
%
% Alexandre Matov, November 7th, 2002
% Modified by Sylvain Berlemont, 2010

Imax=zeros(size(IG)); % prepearing Imax for substraction

% Replace loop
validCands = ([cands(:).status] == 1) & ([cands(:).deltaI] > 0);
if any(validCands)
    validLmax = vertcat(cands(validCands).Lmax);
    validDeltaI = [cands(validCands).deltaI];
    validIdx = sub2ind(size(IG), validLmax(:,1), validLmax(:,2));
    Imax(validIdx) = validDeltaI;
end
% % setting the value of Imax to the values from delataI, i.e. taking into account the background
% s=length(cands);
%
% % counter
% nmb=0;
%
% for i=1:s
%     if cands(i).status==1 && cands(i).deltaI>0 % status flag - a local maximun is significant or not
%         Imax(cands(i).Lmax(1),cands(i).Lmax(2))=cands(i).deltaI; % the intensity at the speckle position is set to "delta I"
%         nmb=nmb+1;
%     end
%
% end

% the masked with a GK local maxima points with intensity delta.I for the raw data image
Imaxima=filterGauss2D(Imax,SIG);

% substruction
Inew=IG-Imaxima;