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

% Initialize cands structure
% SB: this is useless
% cands=struct(...
%     'Lmax',[0 0],...                 % Local maximum position - [y x]
%     'Bkg1',[0 0],...                 % First local minimum position - [y x]
%     'Bkg2',[0 0],...                 % Second local minimum position - [y x]
%     'Bkg3',[0 0],...                 % Third local minimum position - [y x]
%     'ILmax',0,...                    % Local maximum intensity
%     'IBkg',0,...                     % Mean background intensity
%     'deltaI',0,...                   % Intensity difference: ILmax-IBkg
%     'deltaICrit',0,...               % Critical intensity difference as calculated with the noise model
%     'sigmaLmax',0,...                % Error on local maximum intensity
%     'sigmaBkg',0,...                 % Error on background intensity 
%     'status',0);                     % Significance of the local maximum: 1, speckle; 0, weak local maximum

% Collect lMax coordinates into matrics
[y x]=find(lMax);
pMax=[y x];

if isempty(pMax)
   return;
end

% SB: We need to remove dSet since virtual points on the cell edge has been
% added.

% Setting vertex coordinates
%dSet=[1,1;imgSize(1),1;imgSize(1),imgSize(2);1,imgSize(2)];

% Collect lMin coordinates into matrices
% SB: check whether the pixels with -1000 value are caught by find. GOOD
[y x]=find(lMin);
pMin=[y x];

% Attach dSet to lMin
%pMin=cat(1,pMin,dSet);

% Delaunay triangulation
triMin=delaunay(pMin(:,1),pMin(:,2),{'Qt'}); % New delaunay function (MATLAB 7), 'Qt' -> triangulated output

% Search triangles
triangles=tsearch(pMin(:,1),pMin(:,2),triMin,pMax(:,1),pMax(:,2));

% SB: the case triangles(i) == NaN should now happens only at the image
% border and not anymore at the cell edge. => Discard this case (no Bkg ==
% [-1 -1] anymore.

% Store information

% SB: Replace the loop
validTri = 1:numel(triangles);
validTri = validTri(~isnan(triangles));
n = ones(numel(validTri), 1);

cands = struct(...
    'Lmax', mat2cell(pMax(validTri,:), n), ...
    'Bkg1', mat2cell(pMin(triMin(triangles(validTri),1),:), n),...
    'Bkg2', mat2cell(pMin(triMin(triangles(validTri),2),:), n),...
    'Bkg3', mat2cell(pMin(triMin(triangles(validTri),3),:), n));

% for i=1:size(triangles,1)
%    if ~isnan(triangles(i))  % If NaN -> no triangle found
%       cands(i).Lmax=pMax(i,:);
%       % Read local maxima positions into the 3x2 matrix Bkg
%       Bkg=[pMin(triMin(triangles(i),:),1) pMin(triMin(triangles(i),:),2)];
%       % Store positions into the cands structure
%       cands(i).Bkg1=Bkg(1,:);
%       cands(i).Bkg2=Bkg(2,:);
%       cands(i).Bkg3=Bkg(3,:);
%    else   % Mark failed Delaunay triangulation
%       cands(i).Lmax=pMax(i,:);
%       cands(i).Bkg1=[-1 -1];
%       cands(i).Bkg2=[-1 -1];
%       cands(i).Bkg3=[-1 -1];
%    end
%end;
