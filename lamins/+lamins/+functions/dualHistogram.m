function [ h ] = dualHistogram( data, labels, edges)
%plotFaceStatistics
%
% skeletons - an array of lamins.classes.Skeletons
% scale - normalization factor to scale statistic by

if(nargin < 3)
    edges = {};
else
    edges = { edges };
end

% data = {data1 data2};

rowHeight = 0.75;
rowSpacing = 1 - rowHeight;
rowToRow = rowHeight + rowSpacing;

h = figure;
hold on;

all.data = [ data{:} ];
all.median = median(all.data);
all.std = std(all.data);
all.mean = mean(all.data);
all.cutoff = all.median + 5*all.std;
all.outliers = all.data > all.cutoff;
[all.N , all.edges] = histcounts(all.data(~all.outliers),edges{:});
all.cutoff = max(all.cutoff, all.edges(end));
all.histmax = max(all.N);
all.centers = ( all.edges(1:end-1) + all.edges(2:end) )/2;

% plot([all.edges(1) all.centers all.edges(end) all.cutoff([ 1 1])], ...
%     [0 all.N/all.histmax*rowHeight 0 0 sum(all.outliers)/all.histmax]);
% % scatter(all.data(all.outliers),zeros(sum(all.outliers),1));
% % line(all.cutoff([ 1 1]),[0 sum(all.outliers)/all.histmax],'Marker','+');
% line(all.median([ 1 1]),[0 rowHeight],'Color','r');
% line(all.mean([ 1 1]),[0 rowHeight],'Color','g');

% baseline = 0;

cm = lines(7);
cm = cm([5 7],:);

components(length(data)) = struct();

for i=1:length(data)
%     baseline = i*rowToRow;
    baseline = (length(data)-i)*rowToRow;
    [components(i).data] = [data{i}];
    components(i).median = median(components(i).data(:));
    components(i).std = std(components(i).data(:));
    components(i).mean = mean(components(i).data(:));
    components(i).top = prctile(components(i).data(:),95);
    components(i).outliers = components(i).data > all.cutoff;
    [components(i).N,components(i).edges] = histcounts(components(i).data(~components(i).outliers),all.edges);
    components(i).histmax = max(components(i).N);
    components(i).centers = ( components(i).edges(1:end-1) + components(i).edges(2:end) )/2;
    
    line([0 all.cutoff],[baseline baseline],'Color',[0.5 0.5 0.5]);
    if(builtin('all',round(all.centers) == all.centers))
        stems = [components(i).N/components(i).histmax*rowHeight sum(components(i).outliers)/components(i).histmax];
        stems = [zeros(size(stems)); stems; NaN(size(stems))];
        stems = stems + baseline;
        X = [components(i).centers all.cutoff; components(i).centers all.cutoff];
        X(3,:) = NaN;
        X = joinColumns(X);
        line(X,stems(:),'Color',cm(mod(i-1,2)+1,:));
%         stem([components(i).edges(1) components(i).centers components(i).edges(end) all.cutoff([1 1])], ...
%             [0 components(i).N/components(i).histmax*rowHeight 0 0 sum(components(i).outliers)/components(i).histmax]+baseline,'Color',cm(mod(i-1,2)+1,:),'BaseValue',baseline);
    else
        plot([components(i).edges(1) components(i).centers components(i).edges(end) all.cutoff([1 1])], ...
            [0 components(i).N/components(i).histmax*rowHeight 0 0 sum(components(i).outliers)/components(i).histmax]+baseline,'Color',cm(mod(i-1,2)+1,:));
        % scatter(areas(i).data(areas(i).outliers),zeros(sum(areas(i).outliers),1));
    %     line(all.cutoff([ 1 1]),[0 sum(areas(i).outliers)/areas(i).histmax]+baseline,'Marker','+');
        line(components(i).median([ 1 1]),[0 rowHeight]+baseline,'Color','k','LineStyle','-');
        line(components(i).mean([ 1 1]),[0 rowHeight]+baseline,'Color','k','LineStyle','--');
    end
    
end

ylabel('Frequency');
ticks = (0:length(data)) + rowHeight/2;
ticks(2,:) = (0:length(data)) + rowHeight;
set(gca,'YTick',ticks(:));
histMaxima = round([components.histmax] ./ [arrayfun(@(x) length(x.data),components)],3);
allLabels = [labels ; num2cell(histMaxima)];
allLabels = fliplr(allLabels);
set(gca,'YTickLabel',allLabels(:));
set(gcf,'Position',[369 33 290 1064]);
hold off;
% set(gca,'YTickLabel',[{'Combined'}; num2cell(num2str((1:length(data))'),2)]);



% set(gca,'YTickLabel',[{'Combined'}; cellfun(@(MD) strrep(MD.movieDataFileName_,'_','\_'), ML.movies_,'Unif',false)'])
% xlabel('Face Area (\mum^2)')

end

