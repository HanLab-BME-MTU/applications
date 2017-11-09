function [h] = visualizeMotionFields(I,ROI,dxs,dys,patchSize,reduceMotionResolution)
%% 
% Input: 
%   I - image
%   ROI - segmentaion mask (1 - cellular, 0 - background)
%   dx,dy - motion fiels
%   patchSize - in pixels
%   reduceMotionResolution - reduced resolution for motion fields (0.5 --> 2x2 patches)
%   outFname - to save the file
%
% Output: handle to a figure visualization of motion fields overlayed on the image
%
% Assaf Zaritsky, December 2014

% % one channel
% if size(I,3) > 1
%     I = I(:,:,1);
% end

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

h = figure;
imagesc(I);
colormap(gray);
hold on;
qh = quiver(xInds',yInds',us',vs',1.5,'-r','LineWidth',1,'MarkerFaceColor','r');
% qh = quiver(xInds',yInds',us',vs',2,'-r','LineWidth',2,'MarkerFaceColor','r');

% adjust_quiver_arrowhead_size(qh, 2);

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
set(haxes,'XTickLabel',[]);
set(haxes,'YTickLabel',[]);
set(h,'Color','none');
hold off;

end