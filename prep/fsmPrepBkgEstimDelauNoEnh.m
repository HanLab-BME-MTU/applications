function [cands,triMin,pMin]=fsmPrepBkgEstimDelauNoEnh(imgSize,lMax,lMin)
% fsmPrepBkgEstimDelauNoEnh uses Delaunay triangulation to assign 3 local minima to every local maximum
%
% SYNOPSIS   [cands,triMin,pMin]=fsmPrepBkgEstimationDelaunay(imgSize,lMax,lMin)
%
% INPUT      imgSize   :   image size [y x]
%            lMax      :   local max map (the output of the locMax2D function)
%            lMin      :   local min max (the output of the locMin2D function)
%
% OUTPUT     triMin    :   Delaunay triangulation
%            pMin      :   attached dSet (vertex coordinates) to lMin
%            cands     :   structure containing statistical information for each local maximum
%                          (see help for detail) - fsmPrepBkgEstimationDelaunay only stores 
%                          local maximum and local minimum coordinates
%
% Note about previous 'Enhanced triangulation' option:
%
% As of Matlab 6.1, the delaunay function has been enhanced and does no longer accept the 
% 'fuzz' parameter.
%
% DEPENDENCES
%
% Aaron Ponti, October 23th, 2002

% DEBUG
showInfo=0;

% Initialize cands structure
cands=struct(...
    'Lmax',[0 0],...                 % Local maximum position - [y x]
    'Bkg1',[0 0],...                 % First local minimum position - [y x]
    'Bkg2',[0 0],...                 % Second local minimum position - [y x]
    'Bkg3',[0 0],...                 % Third local minimum position - [y x]
    'ILmax',0,...                    % Local maximum intensity
    'IBkg',0,...                     % Mean background intensity
    'deltaI',0,...                   % Intensity difference: ILmax-IBkg
    'deltaICrit',0,...               % Critical intensity difference as calculated with the noise model
    'sigmaLmax',0,...                % Error on local maximum intensity
    'sigmaBkg',0,...                 % Error on background intensity 
    'status',0);                     % Significance of the local maximum: 1, speckle; 0, weak local maximum


    
% Define 'd' as the maximum side length for a triangle to be kept
d=30;

% Setting vertex coordinates
dSet=[1,1;imgSize(1),1;imgSize(1),imgSize(2);1,imgSize(2)];

% Collect coordinates into matrices
[y x]=find(lMax);
pMax=[y x];
if isempty(pMax)
   info=[];
   return;
end
[y x]=find(lMin);
pMin=[y x];

% Attach dSet to lMin
pMin=cat(1,pMin,dSet);
y=cat(1,y,dSet(:,1));
x=cat(1,x,dSet(:,2));

% Delaunay triangulation
triMin=delaunay(pMin(:,1),pMin(:,2),{'Qt'}); % New delaunay function (MATLAB 7), 'Qt' -> triangulated output

% Search triangles
triangles=tsearch(pMin(:,1),pMin(:,2),triMin,pMax(:,1),pMax(:,2));

if showInfo==1
	% Plot loc min
	figure;
	%imshow(img),hold on;
	%plot(pMin(:,2),pMin(:,1),'r.');
	hold on;
	% Plot loc max
	handle=plot(pMax(:,2),pMax(:,1),'b.');
    set(handle,'MarkerSize',20);
	% Draw selected Delaunay triangles
	for i=1:size(triangles,1)
		if ~isnan(triangles(i))
			fVecX = pMin([triMin(triangles(i),:),triMin(triangles(i),1)],2);
			fVecY = pMin([triMin(triangles(i),:),triMin(triangles(i),1)],1);
			handle=plot(fVecX,fVecY,'k-');
			set(handle,'LineWidth',2);
            set(handle,'Marker','.')
            set(handle,'MarkerSize',20);
            set(handle,'MarkerSize',15);
            set(handle,'MarkerEdgeColor','r');
            set(handle,'MarkerFaceColor','r');
			%pause;
		end
	end;
	pause;
	hold off;
	close;
end

% Store information
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
