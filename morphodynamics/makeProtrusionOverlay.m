function [] = makeProtrusionOverlay(smoothedEdge,img,minIm,maxIm,colorUnit )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[nFrame tmp]=size(smoothedEdge);
cmap=jet(nFrame);
%Load the image(s).

h=figure;    
imgRGB(:,:,1)=img;imgRGB(:,:,2)=img;imgRGB(:,:,3)=img;
imgRGB = (double(imgRGB)-minIm)/(maxIm-minIm);
imgRGB(imgRGB>1)=1;
imgRGB(imgRGB<0)=0;
imagesc(imgRGB), axis image, axis off;
caxis([1 nFrame]);
hc=colorbar;
xlabel(hc,colorUnit);
hold on;
    
for iFrame = 1:nFrame

    plot(smoothedEdge{iFrame}(:,1),smoothedEdge{iFrame}(:,2),'Color',cmap(iFrame,:))
   
end
hold off