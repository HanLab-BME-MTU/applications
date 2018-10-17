function [h] = makeProtrusionOverlay(smoothedEdge,img,minIm,maxIm,timeInterval)
%function [] =
%makeProtrusionOverlay(smoothedEdge,img,minIm,maxIm,colorUnit) creates the
%colored edge overlay on top of the single image, img.
%   input:  smoothedEdge: edges from the protrusion process output.
%           img: a single image
%           timeInterval: time interval in sec
%   output: h: figure handle
% Look also: makeMovieProtrusionOverlay(MD)
% Kwonmoo Lee, modified by Sangyoon Han. Oct 2018
if nargin<5
    timeInterval=[];
end

nFrame = size(smoothedEdge,1);
cmap=parula(nFrame);
%Load the image(s).

h=figure;    
imgRGB(:,:,1)=img;imgRGB(:,:,2)=img;imgRGB(:,:,3)=img;
imgRGB = (double(imgRGB)-minIm)/(maxIm-minIm);
imgRGB(imgRGB>1)=1;
imgRGB(imgRGB<0)=0;
imagesc(imgRGB), axis image, axis off;
hold on;
    
for iFrame = 1:nFrame

    plot(smoothedEdge{iFrame}(:,1),smoothedEdge{iFrame}(:,2),'Color',cmap(iFrame,:))
   
end

if isempty(timeInterval)
    caxis([1 nFrame]);
    colorUnit='frames';
else
    caxis([0 (nFrame-1)*timeInterval/60]);
    colorUnit='min';
end
hc=colorbar;
xlabel(hc,colorUnit);

hold off