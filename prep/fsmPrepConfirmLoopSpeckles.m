function cands=fsmPrepConfirmLoopSpeckles(Inew,noiseParam,triMin,pMin,IG,userROIbw)

% fsmPrepConfirmLoopSpeckles uses statistical tests to confirm the significance of detected speckles
% of higher than one hierarchical level (in the main loop)
%
% SYNOPSIS   cands=fsmPrepConfirmLoopSpeckles(Inew,noiseParam,triMin,pMin,IG)
%
% INPUT      Inew       :  substracted image
%            noiseParam :  noise parameters for statistical speckle selection
%            triMin     :  set of Delaunay triangles
%            pMin       :  attached dSet (vertex coordinates) to lMin
%            IG         :  original (filtered) image
%
% OUTPUT     cands      :  cands structure (see fsmPrepTestLocalMaxima.m)
%
%
%
% DEPENDENCES   fsmPrepConfirmLoopSpeckles uses { locmax2d, tsearch,fsmPrepTestLocalMaxima }
%               fsmPrepConfirmLoopSpeckles is used by { fsmPrepMainSecondarySpeckles }
%
% Alexandre Matov, April 2nd, 2003
% Sylvain Berlemont, Jan 2010

% find the local maxima  
Imax=locmax2d(Inew,[5,5]);

if nargin == 6 && ~isempty(userROIbw)
    % Mask Imax
    Imax=Imax.*userROIbw;
end

% find the coordinates/positions of the initial local maxima before
% significance test (for comparision)
[yi,xi]=find(ne(Imax,0));
pMax=[yi,xi];

% Assign local maxima to local minimum triangles
triangles=tsearch(pMin(:,1),pMin(:,2),triMin,pMax(:,1),pMax(:,2));

% Store information into cands structure
validTri = 1:numel(triangles);
validTri = validTri(~isnan(triangles));
n = ones(numel(validTri), 1);

cands = struct(...
    'Lmax', mat2cell(pMax(validTri,:), n), ...
    'Bkg1', mat2cell(pMin(triMin(triangles(validTri),1),:), n),...
    'Bkg2', mat2cell(pMin(triMin(triangles(validTri),2),:), n),...
    'Bkg3', mat2cell(pMin(triMin(triangles(validTri),3),:), n));

% analyze speckles - validate, locmax, locmin...
cands = fsmPrepTestLocalMaxima(Inew,cands,noiseParam,IG);  


