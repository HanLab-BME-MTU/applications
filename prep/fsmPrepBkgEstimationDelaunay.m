function [cands,triMin,pMin]=fsmPrepBkgEstimationDelaunay(lMax,lMin)
% fsmPrepBkgEstimationDelaunay uses Delaunay triangulation to assign 3 local minima to every local maximum
%
% SYNOPSIS   [cands,triMin,pMin]=fsmPrepBkgEstimationDelaunay(lMax,lMin)
%
% INPUT      lMax      :   local max map (the output of the locMax2D function)
%            lMin      :   local min max (the output of the locMin2D function)
%
% OUTPUT     triMin    :   Delaunay triangulation
%            pMin      :   attached dSet (vertex coordinates) to lMin
%            cands     :   structure containing statistical information for each local maximum
%                          (see help for detail) - fsmPrepBkgEstimationDelaunay only stores 
%                          local maximum and local minimum coordinates
%
% DEPENDENCES
%
% Aaron Ponti, October 23th, 2002
% Sylvain Berlemont, Jan 2010

% Collect lMax coordinates into matrics
[y x]=find(lMax);
pMax=[y x];

if isempty(pMax)
   return;
end

% Collect lMin coordinates into matrices
[y x]=find(lMin);
pMin=[y x];

% Delaunay triangulation
% (New delaunay function (MATLAB 7), 'Qt' -> triangulated output)
triMin=delaunay(pMin(:,1),pMin(:,2),{'Qt'});

% Search triangles
triangles=tsearch(pMin(:,1),pMin(:,2),triMin,pMax(:,1),pMax(:,2));

% Store information
validTri = 1:numel(triangles);
validTri = validTri(~isnan(triangles));
n = ones(numel(validTri), 1);

cands = struct(...
    'Lmax', mat2cell(pMax(validTri,:), n), ...
    'Bkg1', mat2cell(pMin(triMin(triangles(validTri),1),:), n),...
    'Bkg2', mat2cell(pMin(triMin(triangles(validTri),2),:), n),...
    'Bkg3', mat2cell(pMin(triMin(triangles(validTri),3),:), n));
