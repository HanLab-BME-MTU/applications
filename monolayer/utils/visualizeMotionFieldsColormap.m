function [] = visualizeMotionFieldsColormap(I,ROI,dxs,dys,patchSize,reduceMotionResolution,clim,outFname)
%% 
% Input: 
%   I - image
%   ROI - segmentaion mask (1 - cellular, 0 - background)
%   dx,dy - motion fiels
%   patchSize - in pixels
%   reduceMotionResolution - reduced resolution for motion fields (0.5 --> 2x2 patches)
%   clim - magnitude limits to classify the vectors
%
% Output: writes figure to outFname
%
% Assaf Zaritsky, May 2016

% % one channel
% if size(I,3) > 1
%     I = I(:,:,1);
% end

close all;
 
addpath(genpath('/home2/azaritsky/code/common/graphics/quiverColormap.m'));

dxs1 = imresize(dxs,reduceMotionResolution);
dys1 = imresize(dys,reduceMotionResolution);
dxs2 = imresize(dxs1,size(I(:,:,1)));
dys2 = imresize(dys1,size(I(:,:,1)));

MASK1 = ROI & ~isnan(dxs2) & ~isnan(dys2);

reducedResolutionPatchSize = round(patchSize * 1/reduceMotionResolution);
xs = (1+round(reducedResolutionPatchSize/2)) : reducedResolutionPatchSize : (size(I,2)-round(reducedResolutionPatchSize/2));
ys = (1+round(reducedResolutionPatchSize/2)) : reducedResolutionPatchSize : (size(I,1)-round(reducedResolutionPatchSize/2));
MASK2 = false(size(MASK1));
MASK2(ys,xs) = true;
MASK = MASK1 & MASK2;

[yInds,xInds] = ind2sub(size(MASK),find(MASK));

us = dxs2(MASK);
vs = dys2(MASK);


%% Figure
% [left bottom width height]
FPosition = [0 0 350 350];
APosition = [0.2 0.2 0.75 0.75]; 

% % % Arrange colormap
% % magnitude = sqrt(us.^2+vs.^2);
% % nColors = 60;
% % cmap = colormap(hsv(nColors));
% % vColor=floor(scaleContrast(magnitude,clim,[1 nColors]));
% % vColor(vColor<1)=1;
% % vColor(vColor>nColors)=nColors;
% % vIndex= unique(vColor);


% I = zeros(size(I));

h = figure;
imagesc(I);
colormap(gray);
hold on;


%%


quiver(xInds',yInds',us',vs',1.5,'-r','LineWidth',1,'MarkerFaceColor','r'); % ,,'MarkerEdgeColor','k' 
% adjust_quiver_arrowhead_size(qh, 1.5);
% % % Create array of quiverplots
% % for i=1:numel(vIndex)
% %     idx = find(vColor==vIndex(i));
% %     curColor = cmap(i,:);    
% %     quiver(xInds(idx)',yInds(idx)',us(idx)',vs(idx)',1.5,'-','LineWidth',1,'MarkerFaceColor',curColor,'MarkerEdgeColor',curColor); 
% % end

%%



% qh = quiverColormap(xInds',yInds',us',vs','LineWidth',1,'CLim',clim);% 1.5,'-','MarkerFaceColor','r',
% qh = quiverColormapMonolayerUtils(xInds',yInds',us',vs','LineWidth',1,'CLim',clim);% 1.5,'-','MarkerFaceColor','r',
% qh = quiver(xInds',yInds',us',vs',2,'-r','LineWidth',2,'MarkerFaceColor','r');



%%
% children=get(qh,'children'); % retrieve the plot-children - second element are the arrow tips
% 
% XData=get(children(2),'XData'); % retrieve the coordinates of the tips
% YData=get(children(2),'YData');
% 
% delete(children(2))  % delete old arrow tips
% 
% for l=1:4:length(XData)-3   % paint new arrow tips, skipping the NaN-values
%     ArrowTips((l-1)/4+1)=fill(XData(l:l+2),YData(l:l+2),'r');
% end
%%

% set(qh,'MarkerSize',10)
haxes = get(h,'CurrentAxes');
set(haxes,'XTick',[]);
set(haxes,'YTick',[]);
set(haxes,'XTickLabel',[]);
set(haxes,'YTickLabel',[]);
set(h,'Color','w');
% set(h,'Position',FPosition,'PaperPositionMode','auto');
% axisHandle= findobj(h,'type','axes');
% set(axisHandle,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
% set(get(axisHandle,'XLabel'),'FontSize',fontsize); set(get(axisHandle,'YLabel'),'FontSize',fontsize);
axis equal; % no image stretching
% set(haxes,'box','off');
axis off;
axis tight;
hold off;
export_fig(outFname);
end