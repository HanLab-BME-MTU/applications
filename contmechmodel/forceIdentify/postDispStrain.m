%This is a script file that displays the colormap of strain tensor
% overlayed on the original cell image.

%The time step to display.
dispTimeStep = 1;

%Get the dimension of the cell image.
pixelX = [1:size(cellImg,2)];
pixelY = [1:size(cellImg,1)];

%Calculate the intensity map of adhesion force. The intensity is between 0 and
% 1.
strain11 = fnval(ppStrain{1,dispTimeStep},{pixelY pixelX});
strain12 = fnval(ppStrain{2,dispTimeStep},{pixelY pixelX});
strain22 = fnval(ppStrain{3,dispTimeStep},{pixelY pixelX});
compStrnMap  = strain11+strain22;
shearStrnMap = strain12;

%Scale so that the max strain is one
strainScale  = max(max(abs(compStrnMap(:))),max(abs(shearStrnMap(:))));
compStrnMap  = (compStrnMap/strainScale+1)/2;
shearStrnMap = (shearStrnMap/strainScale+1)/2;

cMap     = colormap('jet');
cMapRang = size(cMap,1)-1;
[compStrnImg,compPix]   = applyColorMap(cellImg,compStrnMap*cMapRang, ...
   [0 cMapRang],cMap,1);
[shearStrnImg,shearPix] = applyColorMap(cellImg,shearStrnMap*cMapRang, ...
   [0 cMapRang],cMap,1);

cutOffI = floor(cMapRng/4);
%cMap(1:cutOffI,:) = ones(cutOffI,1)*cMap(cutOffI,:);
%cMap(end-cutOffI+1:end,:) = ones(cutOffI,1)*cMap(end-cutOffI+1,:);
postplot(fem,'tridata','-(u1x+u2y)','tribar','on','trimap',cMap);

