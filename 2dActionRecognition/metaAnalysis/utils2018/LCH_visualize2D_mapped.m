function [] = LCH_visualize2D_mapped(mappedFeats,labels,outFname)
inds = find(~cellfun(@isempty,labels));
mappedFeats = mappedFeats(inds,:);
labels = labels(inds);
h = figure;
hold on;
gscatter(mappedFeats(:,1), mappedFeats(:,2), labels');
% legend(unique(labels));
hold off;
saveas(h,outFname);

% saveas(h,[outFname(1:end-4) '.eps']);
% save([outFname(1:end-4) '.mat'],'mappedFeats','labels');

close all;
end