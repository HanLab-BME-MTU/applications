function [ h, rps, all, areas ] = plotFaceStatistics( skeletons, property, scale )
%plotFaceStatistics
%
% skeletons - an array of lamins.classes.Skeletons
% scale - normalization factor to scale statistic by

rowHeight = 0.75;
rowSpacing = 1 - rowHeight;
rowToRow = rowHeight + rowSpacing;

if(nargin < 2)
    property = 'Area';
end

if(nargin < 3)
    scale = 1000;
end


h = figure;
hold on;

rps = arrayfun(@(S) regionprops(S.faces,property),skeletons,'UniformOutput',false);

all.rp = vertcat(rps{:});
numFaces = numel(all.rp);
[all.area.data] = [all.rp.(property)]/scale;
all.area.median = median(all.area.data(:));
all.area.std = std(all.area.data(:));
all.area.mean = mean(all.area.data(:));
all.area.top = prctile(all.area.data(:),95);
all.area.cutoff = all.area.median + 3*all.area.std;
% all.area.cutoff = 1;
all.area.outliers = all.area.data > all.area.cutoff;
[all.area.N,all.area.edges] = histcounts(all.area.data(~all.area.outliers));
all.area.cutoff = max(all.area.cutoff,all.area.edges(end));
all.area.histmax = max(all.area.N);
all.area.centers = ( all.area.edges(1:end-1) + all.area.edges(2:end) )/2;


plot([all.area.edges(1) all.area.centers all.area.edges(end) all.area.cutoff([ 1 1])], ...
    [0 all.area.N/all.area.histmax*rowHeight 0 0 sum(all.area.outliers)/all.area.histmax]);
% scatter(all.area.data(all.area.outliers),zeros(sum(all.area.outliers),1));
% line(all.area.cutoff([ 1 1]),[0 sum(all.area.outliers)/all.area.histmax],'Marker','+');
line(all.area.median([ 1 1]),[0 rowHeight],'Color','r');
line(all.area.mean([ 1 1]),[0 rowHeight],'Color','g');


areas(length(skeletons)) = struct();

for i=1:length(skeletons)
    baseline = i*rowToRow;
    [areas(i).data] = [rps{i}.(property)]/scale;
    areas(i).median = median(areas(i).data(:));
    areas(i).std = std(areas(i).data(:));
    areas(i).mean = mean(areas(i).data(:));
    areas(i).top = prctile(areas(i).data(:),95);
    areas(i).outliers = areas(i).data > all.area.cutoff;
    [areas(i).N,areas(i).edges] = histcounts(areas(i).data(~areas(i).outliers),all.area.edges);
    areas(i).histmax = max(areas(i).N);
    areas(i).centers = ( areas(i).edges(1:end-1) + areas(i).edges(2:end) )/2;
    
    line([0 all.area.cutoff],[baseline baseline],'Color',[0.5 0.5 0.5]);
    plot([areas(i).edges(1) areas(i).centers areas(i).edges(end) all.area.cutoff([1 1])], ...
        [0 areas(i).N/areas(i).histmax*rowHeight 0 0 sum(areas(i).outliers)/areas(i).histmax]+baseline);
    % scatter(areas(i).data(areas(i).outliers),zeros(sum(areas(i).outliers),1));
%     line(all.area.cutoff([ 1 1]),[0 sum(areas(i).outliers)/areas(i).histmax]+baseline,'Marker','+');
    line(areas(i).median([ 1 1]),[0 rowHeight]+baseline,'Color','r');
    line(areas(i).mean([ 1 1]),[0 rowHeight]+baseline,'Color','g');
    
end

set(gca,'YTick',0:length(skeletons));
set(gca,'YTickLabel',[{'Combined'}; num2cell(num2str((1:length(skeletons))'),2)]);

% set(gca,'YTickLabel',[{'Combined'}; cellfun(@(MD) strrep(MD.movieDataFileName_,'_','\_'), ML.movies_,'Unif',false)'])
% xlabel('Face Area (\mum^2)')

end

