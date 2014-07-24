function [Imin,deltaI,k,sigmaDiff,sigmaMax,sigmaMin,status]=testSpeckleSignificance(Imax,Imin,k,sigmaD,PoissonNoise,I0)
% testSpeckleSignificance tests the significance of a local maximum in the normal case of successful Delaunay
%
% SYNOPSIS   [Imin,deltaI,k,sigmaDiff,sigmaMax,sigmaMin,status]=testSpeckleSignificance(Imax,Imin,k,sigmaD,PoissonNoise,I0)
%
% INPUT 
%        Imax : local maximum peak intensity
%        Imin : local background intensity
%        k: t-test confidence probability
%        sigmaD :
%        poissonNoise:
%        I0   : average background intensity
%
% OUTPUT
%
% References:
% A. Ponti et al., Biophysical J., 84 336-3352, 2003.

% Sebastien Besson, Aug 2011
% Adapated from fsmPrepTestSpeckleSignif.m

% Calculate error on Imax (noise model)
Imax = max(Imax,I0);
sigmaMax=sqrt(sigmaD^2+PoissonNoise*(Imax-I0));

% Calculate error on Imin (noise model) - sigmaMin is multiplied by 1/sqrt(3) 
%   because Imin is the mean of 3 intensity values
if (Imin-I0)<0, Imin=I0; end
sigmaMin=(1/sqrt(3))*sqrt(sigmaD^2+PoissonNoise*(Imin-I0));

% Calculate difference and error
deltaI=Imax-Imin;
sigmaDiff=sqrt(sigmaMax^2+sigmaMin^2);

% Check for the validity of the speckle
status = deltaI >= k * sigmaDiff;