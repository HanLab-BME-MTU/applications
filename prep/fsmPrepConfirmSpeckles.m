function [yi,xi,y,x,Imax,cands,triMin,pMin]=fsmPrepConfirmSpeckles(IG,Imin,noiseParam,enhTriang)

% fsmPrepConfirmSpeckles uses statistical tests to confirm the significance of detected speckles
%
%
% SYNOPSIS   [yi,xi,y,x,Imax,cands]=fsmPrepConfirmSpeckles(IG,Imin,noiseParam)
%
% INPUT      IG         :  filtered image
%            Imin       :  local minima
%            noiseParam :  noise parameters for statistical speckle selection
%            enhTriang  :  turns on enhanced triangulation for Matlab Version < 6.5
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
% DEPENDENCES   fsmPrepConfirmSpeckles uses { fsmPrepBkgEstimationDelaunay, fsmPrepTestLocalMaxima, locmax2d }
%               fsmPrepConfirmSpeckles is used by { fsmPrepMainSecondarySpeckles }
%
% Alexandre Matov, November 7th, 2002

% find the local maxima  
Imax=locmax2d(IG,[5,5]);

% find the coordinates/positions of the initial local maxima before significance test (for comparision)
[yi,xi]=find(ne(Imax,0));

% Use Delaunay triangulation to assign 3 local minima to each local maximum
matlabVersion=ver('MATLAB');

if str2num(matlabVersion.Version)<6.1
    % Use enhanced triangulation
    [cands,triMin,pMin]=fsmPrepBkgEstimationDelaunay(size(IG),Imax,Imin,enhTriang); % Finds 3 loc min around each loc max
else
    % Use delaunay triangulation of Matlab version >=6.1
    [cands,triMin,pMin]=fsmPrepBkgEstimDelauNoEnh(size(IG),Imax,Imin); % Finds 3 loc min around each loc max
end

% analyze speckles - validate, locmax, locmin...
[Imax,cands]=fsmPrepTestLocalMaxima(IG,Imax,cands,noiseParam,IG);  

% find the coordinates/positions of the local maxima after selecting only the significant local maxima/speckles
[y,x]=find(ne(Imax,0));

