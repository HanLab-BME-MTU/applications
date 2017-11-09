function plotTrackClassesCond(data)
% P. Roudot 2016
tmp=arrayfun(@(d) load([d.source 'Analysis' filesep 'lifetimeData.mat']),data);

plotTrackClasses(vertcat(tmp.catIdx));