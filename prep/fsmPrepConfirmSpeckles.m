function [cands,triMin,pMin]=fsmPrepConfirmSpeckles(IG,Imin,noiseParam,userROIbw)

% fsmPrepConfirmSpeckles uses statistical tests to confirm the significance
% of detected speckles
%
% SYNOPSIS   [cands,triMin,pMin]=fsmPrepConfirmSpeckles(IG,Imin,noiseParam)
%
% INPUT      IG         :  filtered image
%            Imin       :  local minima
%            noiseParam :  noise parameters for statistical speckle selection
%            userROIbw  :  (optional) User-defined black-and-white mask to select speckles
%
% OUTPUT     cands      :  cands structure (see fsmPrepTestLocalMaxima.m)
%            triMin     :  set of Delaunay triangles
%            pMin       :  attached dSet (vertex coordinates) to lMin
%
%
%
% DEPENDENCES   fsmPrepConfirmSpeckles uses { fsmPrepBkgEstimationDelaunay, fsmPrepTestLocalMaxima, locmax2d }
%               fsmPrepConfirmSpeckles is used by { fsmPrepMainSecondarySpeckles }
%
% Alexandre Matov, November 7th, 2002
% Sylvain Berlemont, Jan 2010

% Find the local maxima  
Imax=locmax2d(IG,[5,5]);

if nargin == 4 && ~isempty(userROIbw) 
    Imax=Imax.*userROIbw;
end

% Finds 3 loc min around each loc max
[cands,triMin,pMin]=fsmPrepBkgEstimationDelaunay(Imax,Imin);

% analyze speckles - validate, locmax, locmin...
cands = fsmPrepTestLocalMaxima(IG,cands,noiseParam,IG);  


