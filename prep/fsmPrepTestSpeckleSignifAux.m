function [Imin,deltaI,k,sigmaDiff,sigmaMax,sigmaMin,status,A]=fsmPrepTestSpeckleSignifAux(aImg,img,locMax,maskR,stdImg,A,parameters)
% fsmPrepTestSpeckleSignificanceAux tests the significance of a local maximum in the case of failed Delaunay
%
% SYNOPSIS   [Imin,deltaI,k,sigmaDiff,sigmaMax,sigmaMin,status,A]=fsmPrepTestSpeckleSignificanceAux(aImg,img,locMax,maskR,stdImg,A,parameters)
%
%

% Assign noise parameters
k=parameters(1);
sigmaD=parameters(2);
PoissonNoise=parameters(3);
I0=parameters(4);

% Positions of all local maxima
y=locMax(1);x=locMax(2);

% Add maskR offset to the local maximum coordinates
ym=maskR+y; xm=maskR+x;

% Read mask from image
mask=aImg(ym-maskR:ym+maskR,xm-maskR:xm+maskR);

% Intensity of the local maximum
Imax=img(y,x);

% Check - this check is no longer correct since now we work on sutracted images
% if mask(maskR+1,maskR+1)~=Imax
% 	error('Mask is not centered around the local maximum');
% end

% Intensity of the background - mean intensity of an area of size [2maskR,2maskR] around loc max
%      from which one removes the central part (corresponding to the diffracted local maximum peak)
v=maskR-2;
ver=mask(:,[1:v end-(v-1):end]);ver=ver(:);
hor=mask([1:v end-(v-1):end],[v+1:end-v]);hor=hor(:);

% Calculate background intensity from the region selected
Imin=mean(cat(1,ver,hor));

% Calculate
deltaI=Imax-Imin;

% Check for the statistical significance of the local maximum
if deltaI>0.1*stdImg
	if size(A)~=[1 1] % Function called by validateSpeckle2 and not saveSpeckleArray
		A(y,x)=1;
	end
	status=1;
else
	status=0;
end

% Calculate all other values needed (noise model)
if (Imax-I0)<0
	Imax=I0;
	%warning('validateSpeckles: deltaI<0');
end
sigmaMax=sqrt(sigmaD^2+PoissonNoise*(Imax-I0));

% Calculate error on Imin (noise model) - sigmaMin is multiplied by 1/sqrt(3) 
%   because Imin is the mean of 3 intensity values
if (Imin-I0)<0
	Imin=I0;
	%warning('validateSpeckles: deltaI<0');
end
sigmaMin=(1/sqrt(3))*sqrt(sigmaD^2+PoissonNoise*(Imin-I0));

% Calculate error on difference
sigmaDiff=sqrt(sigmaMax^2+sigmaMin^2);
