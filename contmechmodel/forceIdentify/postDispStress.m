%This is a script file that displays the colormap of stress tensor
% overlayed on the original cell image.

%The time step to display.
dispTimeStep = 1;

%Get the dimension of the cell image.
pixelX = [1:size(cellImg,2)];
pixelY = [1:size(cellImg,1)];

%Calculate the intensity map of adhesion force. The intensity is between 0 and
% 1.
stress11 = fnval(ppStress{1,dispTimeStep},{pixelY pixelX});
stress12 = fnval(ppStress{2,dispTimeStep},{pixelY pixelX});
stress22 = fnval(ppStress{3,dispTimeStep},{pixelY pixelX});
compStrsMap  = stress11+stress22;
shearStrsMap = stress12;

%Scale so that the max stress is one
scale  = max(max(abs(compStrsMap(:))),max(abs(shearStrsMap(:))));
compStrsMap  = compStrsMap/scale;
shearStrsMap = shearStrsMap/scale;

cMap    = colormap('jet');
cMapRng = size(cMap,1)-1;
[compStrsImg,compPix] = applyColorMap(cellImg,compStrsMap*cMapRng, ...
   [0 cMapRng],cMap,1);
[shearStrsImg,shearPix] = applyColorMap(cellImg,shearStrsMap*cMapRng, ...
   [0 cMapRng],cMap,1);
%imshow(compressImg,[]);
%imshow(shearImg,[]);

cutOffI = floor(cMapRng/4);
%cMap(1:cutOffI,:) = ones(cutOffI,1)*cMap(cutOffI,:);
%cMap(end-cutOffI+1:end,:) = ones(cutOffI,1)*cMap(end-cutOffI+1,:);
%postplot(fem,'tridata','-2*(lambda+mu)*(u1x+u2y)*stressScale', ...
%   'tribar','on','trimap',cMap);
%postplot(fem,'tridata','-mu*(u1y+u2x)/2*stressScale','tribar','on', ...
%   'trimap',cMap);
