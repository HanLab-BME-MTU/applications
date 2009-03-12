function [yi,xi,y,x,Imax,cands,triMin,pMin]=fsmPrepConfirmSpeckles(IG,Imin,noiseParam,userROIbw)

% fsmPrepConfirmSpeckles uses statistical tests to confirm the significance of detected speckles
%
%
% SYNOPSIS   [yi,xi,y,x,Imax,cands]=fsmPrepConfirmSpeckles(IG,Imin,noiseParam)
%
% INPUT      IG         :  filtered image
%            Imin       :  local minima
%            noiseParam :  noise parameters for statistical speckle selection
%            userROIbw  :  (optional) User-defined black-and-white mask to select speckles
%
% OUTPUT     yi         :  initial local maxima positions before the test
%            xi         :  initial local maxima positions befire the test
%            y          :  positions of the significant local maxima
%            x          :  positions of the significant local maxima
%            Imax       :  the matrix of [y,x]
%            cands      :  cands structure (see fsmPrepTestLocalMaxima.m)
%            triMin     :  set of Delaunay triangles
%            pMin       :  attached dSet (vertex coordinates) to lMin
%
%
%
% DEPENDENCES   fsmPrepConfirmSpeckles uses { fsmPrepBkgEstimDelauNoEnh, fsmPrepTestLocalMaxima, locmax2d }
%               fsmPrepConfirmSpeckles is used by { fsmPrepMainSecondarySpeckles }
%
% Alexandre Matov, November 7th, 2002

if nargin==4
    userROIbw=[]; % Set the optional value userROI to 0
end

% Find the local maxima  
Imax=locmax2d(IG,[5,5]);

if ~isempty(userROIbw)
    % Mask Imax
    Imax=Imax.*userROIbw;
end

% Find the coordinates/positions of the initial local maxima before significance test (for comparision)
[yi,xi]=find(ne(Imax,0));

[cands,triMin,pMin]=fsmPrepBkgEstimDelauNoEnh(size(IG),Imax,Imin); % Finds 3 loc min around each loc max

% analyze speckles - validate, locmax, locmin...
[Imax,cands]=fsmPrepTestLocalMaxima(IG,Imax,cands,noiseParam,IG);  

% find the coordinates/positions of the local maxima after selecting only the significant local maxima/speckles
[y,x]=find(ne(Imax,0));

