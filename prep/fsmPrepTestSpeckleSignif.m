function [Imin,deltaI,k,sigmaDiff,sigmaMax,sigmaMin,status,A]=testSpeckleSignificance(locMaxPos,Imax,Imin,k,sigmaD,PoissonNoise,I0,A)
% fsmPrepTestSpeckleSignificance tests the significance of a local maximum in the normal case of successful Delaunay
%
% SYNOPSIS   [Imin,deltaI,k,sigmaDiff,sigmaMax,sigmaMin,status,A]=testSpeckleSignificance(locMaxPos,Imax,Imin,k,sigmaD,PoissonNoise,I0,A)
%
%

% Calculate error on Imax (noise model)
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

% Calculate difference and error
deltaI=Imax-Imin;
sigmaDiff=sqrt(sigmaMax^2+sigmaMin^2);

% Check for the validity of the speckle
if deltaI>=k*sigmaDiff;
	if size(A)~=[1 1] % Function called by validateSpeckle2 and not saveSpeckleArray
		A(locMaxPos(1),locMaxPos(2))=1;   % Loc max accepted as speckle
	end
	status=1;                         
else
	status=0;                         % Loc max rejected as speckle
end
