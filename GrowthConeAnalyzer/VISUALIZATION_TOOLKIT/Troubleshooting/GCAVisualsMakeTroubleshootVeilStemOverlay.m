function [ h ] = GCAVisualsMakeTroubleshootVeilStemOverlay(img,analInfo,backboneInfo)
% GCAVisualsMakeTroubleshootVeilStemOverlay: 
% INPUT: 
% img 
% analInfo : (single frame) structure with fields 
% backboneInfo: (single frame) structure with fields
% Note: was troubleshootBodyEst until 20150321
[ny,nx] = size(img); 
h = setFigure(nx,ny,'off');

imshow(img,[]); 

bodyEst = analInfo.bodyEst; 

toPlot = backboneInfo.maxNMSLarge; 
hold on 
spy(toPlot,'b'); % maxNMSAll % don't know if making this additive really helps...could try.. right now only looking at one scael

spy(bodyEst.backbone,'r',5); 

erodForBody = analInfo.bodyEst.erodForBody; 

cellfun(@(x) plot(x(:,2),x(:,1),'g'),erodForBody);
% overlay skeleton of thick body regions
skel = analInfo.bodyEst.thinnedBodyAll; 
spy(skel,'m');

% get the final neuriteEdge 
mask = analInfo.masks.neuriteEdge; 
roiYX  = bwboundaries(mask); 
cellfun(@(x) plot(x(:,2),x(:,1),'y'),roiYX); %final 



end



