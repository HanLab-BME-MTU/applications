function [outEdgeL1,outSpeedL1] = showMapsFlowAndProtrusion(curPath,verbose)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

s=load([curPath '/WindowingPackage/window_sampling/Protrusion Speed map - channel 1.mat']);
if (verbose)
    figure, imagesc(squeeze(s.samples.avg(:,1,:))),
    axis xy
    GreenBlackRedColorMap; caxis([-900 900]),colorbar
    title('Protrusion')
end
% flow speed
s2 = load([curPath '/WindowingPackage/window_sampling/Speed map - channel 1.mat']);
if verbose
    figure, imagesc(squeeze(s2.samples.avg(:,1,:))),
    axis xy
    colormap jet; caxis([0 2000]),colorbar
    title('Actin speed')
end
% output two maps
outSpeedL1 = squeeze(s2.samples.avg(:,1,:));
outEdgeL1=squeeze(s.samples.avg(:,1,:));
end

