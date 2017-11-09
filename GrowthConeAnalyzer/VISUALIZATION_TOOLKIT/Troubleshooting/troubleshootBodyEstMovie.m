function [ output_args ] = troubleshootBodyEstMovie(movieData,channelIndex,analInfo,backboneInfo,frame,saveDir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% if have a problematic structure it is helpful to troubleshoot the body
% estimation 

% make a troubleshoot folder 
% if ~isdir([saveDir filesep 'bodyEstTrouble']) 
%     mkdir([saveDir filesep 'bodyEstTrouble']) 
% end 
nImTot = numel(analInfo); 
  %fmt = ['%0' num2str(ceil(log10(nImTot))) 'd'];

%img = double(imread(analInfo(frame).imgPointer)); 
img =movieData.channels_(channelIndex).loadImage(frame); 
img = double(img); 
[ny,nx] = size(img); 
setFigure(nx,ny,'off');


imshow(-img,[]); 


analInfo1 = analInfo(frame); 
bodyEst = analInfo1.bodyEst; 


toPlot = backboneInfo(frame).maxNMSLarge; 
hold on 
spy(toPlot,'b'); % maxNMSAll % don't know if making this additive really helps...could try.. right now only looking at one scael

spy(bodyEst.backbone,'r',5); 


erodForBody = analInfo(frame).bodyEst.erodForBody; 

cellfun(@(x) plot(x(:,2),x(:,1),'g'),erodForBody);
% overlay skeleton of thick body regions
skel = analInfo(frame).bodyEst.thinnedBodyAll; 
spy(skel,'m');

% get the final neuriteEdge 
mask = analInfo(frame).masks.neuriteEdge; 
roiYX  = bwboundaries(mask); 
cellfun(@(x) plot(x(:,2),x(:,1),'y'),roiYX); %final 

pixels = 10/0.216; 
plotScaleBar(pixels,pixels/20,'Color',[0,0,0]);




saveas(gcf,[saveDir filesep 'troubleshootVeilStem' num2str(frame,'%03d') '.fig']); 
saveas(gcf, [saveDir filesep 'troubleshootVeilStem' num2str(frame,'%03d') '.png']); 

close gcf

end



