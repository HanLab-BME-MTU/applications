function [h] = visualizeStressEllipses(I,ROI,sigmaMax,sigmaMin,stressAngle,patchSize,reduceMotionResolution)
%% 
% Input: 
%   I - image
%   ROI - segmentaion mask (1 - cellular, 0 - background)
%   sigmaMax,sigmaMin - max principle stress and shear
%   patchSize - in pixels
%   reduceMotionResolution - reduced resolution for motion fields (0.5 --> 2x2 patches)
%   outFname - to save the file
%
% Output: handle to a figure visualization of motion fields overlayed on the image
%
% Assaf Zaritsky, February 2015

% one channel
% if size(I,3) > 1
%     I = I(:,:,1);
% end

sigmaMax1 = imresize(sigmaMax,reduceMotionResolution);
sigmaMin1 = imresize(sigmaMin,reduceMotionResolution);
stressAngle1 = imresize(stressAngle,reduceMotionResolution);
sigmaMax2 = imresize(sigmaMax1,size(I(:,:,1)));
sigmaMin2 = imresize(sigmaMin1,size(I(:,:,1)));
stressAngle2 = imresize(stressAngle1,size(I(:,:,1)));

MASK1 = ROI & ~isnan(sigmaMax2) & ~isnan(sigmaMin2) & ~isnan(stressAngle2);

reducedResolutionPatchSize = round(patchSize * 1/reduceMotionResolution);
xs = (1+round(reducedResolutionPatchSize/2)) : reducedResolutionPatchSize : (size(I,2)-round(reducedResolutionPatchSize/2));
ys = (1+round(reducedResolutionPatchSize/2)) : reducedResolutionPatchSize : (size(I,1)-round(reducedResolutionPatchSize/2));
MASK2 = false(size(MASK1));
MASK2(ys,xs) = true;
MASK = MASK1 & MASK2;

[yInds,xInds] = ind2sub(size(MASK),find(MASK));

sigma_max = sigmaMax2(MASK);
sigma_min = sigmaMin2(MASK);
stress_angles = stressAngle2(MASK);

h = figure;
imagesc(I);
colormap(gray);
hold on;
for i = 1 : length(xInds)
    hEllipse = ellipsedraw(sigma_max(i)/10,sigma_min(i)/10,xInds(i),yInds(i),stress_angles(i),'-b');
    %     set(hEllipse,'LineWidth',2);
end

% set(qh,'MarkerSize',10)
haxes = get(h,'CurrentAxes');
set(haxes,'XTickLabel',[]);
set(haxes,'YTickLabel',[]);
set(h,'Color','none');
hold off;

end