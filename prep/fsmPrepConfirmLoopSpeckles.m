function [yi,xi,y,x,Imax,cands]=fsmPrepConfirmLoopSpeckles(Inew,noiseParam,enhTriang,triMin,pMin,IG,userROIbw)

% fsmPrepConfirmLoopSpeckles uses statistical tests to confirm the significance of detected speckles
% of higher than one hierarchical level (in the main loop)
%
% SYNOPSIS   [yi,xi,y,x,Imax,cands]=fsmPrepConfirmLoopSpeckles(Inew,noiseParam,enhTriang,triMin,pMin,IG)
%
% INPUT      Inew       :  substracted image
%            noiseParam :  noise parameters for statistical speckle selection
%            enhTriang  :  turns on enhanced triangulation for Matlab Version < 6.5
%            triMin     :  set of Delaunay triangles
%            pMin       :  attached dSet (vertex coordinates) to lMin
%            IG         :  original (filtered) image
%
% OUTPUT     yi         :  initial local maxima positions before the test
%            xi         :  initial local maxima positions befire the test
%            y          :  positions of the significant local maxima
%            x          :  positions of the significant local maxima
%            Imax       :  the matrix of [y,x]
%            cands      :  cands structure (see fsmPrepTestLocalMaxima.m)
%
%
%
% DEPENDENCES   fsmPrepConfirmLoopSpeckles uses { locmax2d, tsearch,fsmPrepTestLocalMaxima }
%               fsmPrepConfirmLoopSpeckles is used by { fsmPrepMainSecondarySpeckles }
%
% Alexandre Matov, April 2nd, 2003

if nargin==6
    userROIbw=[];
end

% find the local maxima  
Imax=locmax2d(Inew,[5,5]);

if ~isempty(userROIbw)
    % Mask Imax
    Imax=Imax.*userROIbw;
end

% find the coordinates/positions of the initial local maxima before significance test (for comparision)
[yi,xi]=find(ne(Imax,0));
pMax=[yi,xi];

% Assign local maxima to local minimum triangles
triangles=tsearch(pMin(:,1),pMin(:,2),triMin,pMax(:,1),pMax(:,2));

% Store information into cands structure
for i=1:size(triangles,1)
   if ~isnan(triangles(i))  % If NaN -> no triangle found
      cands(i).Lmax=pMax(i,:);
      % Read local maxima positions into the 3x2 matrix Bkg
      Bkg=[pMin([triMin(triangles(i),:)],1) pMin([triMin(triangles(i),:)],2)];
      % Store positions into the cands structure
      cands(i).Bkg1=Bkg(1,:);
      cands(i).Bkg2=Bkg(2,:);
      cands(i).Bkg3=Bkg(3,:);
   else   % Mark failed Delaunay triangulation
      cands(i).Lmax=pMax(i,:);
      cands(i).Bkg1=[-1 -1];
      cands(i).Bkg2=[-1 -1];
      cands(i).Bkg3=[-1 -1];
   end
end;

% analyze speckles - validate, locmax, locmin...
[Imax,cands]=fsmPrepTestLocalMaxima(Inew,Imax,cands,noiseParam,IG);  

% find the coordinates/positions of the local maxima after selecting only the significant local maxima/speckles
[y,x]=find(ne(Imax,0));

