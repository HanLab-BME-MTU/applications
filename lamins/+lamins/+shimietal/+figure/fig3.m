%% Initialize

import lamins.plot.*;

hist.area.data = {  stats(21).area(1).all.data, stats(21).area(2).all.data, ...
                    stats(22).area(1).all.data, stats(22).area(2).all.data, ...
                    stats(25).area(1).all.data, stats(25).area(2).all.data, ...
                    stats(27).area(1).all.data, stats(27).area(2).all.data, ...
                    stats(28).area(1).all.data, stats(28).area(2).all.data, ...
                    stats(23).area(1).all.data, stats(23).area(2).all.data};

hist.edgeLength.data = {    stats(21).edgeLength(1).all.data, stats(21).edgeLength(2).all.data, ...
                            stats(22).edgeLength(1).all.data, stats(22).edgeLength(2).all.data, ...
                            stats(25).edgeLength(1).all.data, stats(25).edgeLength(2).all.data, ...
                            stats(27).edgeLength(1).all.data, stats(27).edgeLength(2).all.data, ...
                            stats(28).edgeLength(1).all.data, stats(28).edgeLength(2).all.data, ...
                            stats(23).edgeLength(1).all.data, stats(23).edgeLength(2).all.data};
                
hist.edgesPerVertex.data = {    stats(21).edgesPerVertex(1).all.data, stats(21).edgesPerVertex(2).all.data, ...
                                stats(22).edgesPerVertex(1).all.data, stats(22).edgesPerVertex(2).all.data, ...
                                stats(25).edgesPerVertex(1).all.data, stats(25).edgesPerVertex(2).all.data, ...
                                stats(27).edgesPerVertex(1).all.data, stats(27).edgesPerVertex(2).all.data, ...
                                stats(28).edgesPerVertex(1).all.data, stats(28).edgesPerVertex(2).all.data, ...
                                stats(23).edgesPerVertex(1).all.data, stats(23).edgesPerVertex(2).all.data};
                            
hist.eccentricity.data = {    stats(21).eccentricity(1).all.data, stats(21).eccentricity(2).all.data, ...
                                stats(22).eccentricity(1).all.data, stats(22).eccentricity(2).all.data, ...
                                stats(25).eccentricity(1).all.data, stats(25).eccentricity(2).all.data, ...
                                stats(27).eccentricity(1).all.data, stats(27).eccentricity(2).all.data, ...
                                stats(28).eccentricity(1).all.data, stats(28).eccentricity(2).all.data, ...
                                stats(23).eccentricity(1).all.data, stats(23).eccentricity(2).all.data};
                            
hist.edgesPerFace.data = {    stats(21).edgesPerFace(1).all.data, stats(21).edgesPerFace(2).all.data, ...
                                stats(22).edgesPerFace(1).all.data, stats(22).edgesPerFace(2).all.data, ...
                                stats(25).edgesPerFace(1).all.data, stats(25).edgesPerFace(2).all.data, ...
                                stats(27).edgesPerFace(1).all.data, stats(27).edgesPerFace(2).all.data, ...
                                stats(28).edgesPerFace(1).all.data, stats(28).edgesPerFace(2).all.data, ...
                                stats(23).edgesPerFace(1).all.data, stats(23).edgesPerFace(2).all.data};

% calculated in fig 6
% for i=1:length(stats)
%     for j=1:length(stats(i).perimeter)
%         stats(i).avgEdgeLength(j).all.data = stats(i).perimeter(j).all.data ./ stats(i).edgesPerFace(j).all.data;
%         stats(i).circularity(j).all.data = 4*pi*stats(i).area(j).all.data ./ stats(i).perimeter(j).all.data.^2;
%         for k=1:length(stats(i).perimeter(j).movie)
%             stats(i).circularity(j).movie(k).data = 4*pi*stats(i).area(j).movie(k).data ./ stats(i).perimeter(j).movie(k).data.^2;
%         end
%     end
% end
hist.avgEdgeLength.data = {    stats(21).avgEdgeLength(1).all.data, stats(21).avgEdgeLength(2).all.data, ...
                                stats(22).avgEdgeLength(1).all.data, stats(22).avgEdgeLength(2).all.data, ...
                                stats(25).avgEdgeLength(1).all.data, stats(25).avgEdgeLength(2).all.data, ...
                                stats(27).avgEdgeLength(1).all.data, stats(27).avgEdgeLength(2).all.data, ...
                                stats(28).avgEdgeLength(1).all.data, stats(28).avgEdgeLength(2).all.data, ...
                                stats(23).avgEdgeLength(1).all.data, stats(23).avgEdgeLength(2).all.data};
                
hist.labels = {    'LA'  'LB1'  ...
                   'LA'  'LB2'  ...
                   'LB2' 'LB1'  ...
                   'LC'  'LB1'  ...
                   'LC'  'LB2'  ...
                   'LC'  'LA' };
       
% data2 = reshape(data,2,[]);
       
%% Area
h = lamins.functions.dualHistogram(hist.area.data,hist.labels);
set(h,'PaperPosition',[0 0 3 11]);
print(h,'fig3_area_hist','-dpdf');

[h , hq , qh] = fig3_qq(hist.area.data,hist.labels);
% lim.x = [0 max(cellfun(@(x) x(end),arrayfun(@(x) xlim(x),hq,'Unif',false)))];
% lim.y = [0 max(cellfun(@(x) x(end),arrayfun(@(x) ylim(x),hq,'Unif',false)))];
lim.x = [0 8];
lim.y = [0 8];
for i=1:6
    xlim(hq(i),lim.x);
    ylim(hq(i),lim.y);
end
set(h,'PaperPosition',[0 0 3 11]);
print(h,'fig3_area_qq','-dpdf');

[h , hq , qh] = fig3_qq(hist.area.data,hist.labels,[0 0.6]);
set(h,'PaperPosition',[0 0 3 11]);
print(h,'fig3_area_qq_zoom','-dpdf');


%% Edge Length
h = lamins.functions.dualHistogram(hist.edgeLength.data,hist.labels,[0:40]/sqrt(1000));
set(h,'PaperPosition',[0 0 3 11]);
print(h,'fig3_edgeLength_hist','-dpdf');

[h , hq , qh] = fig3_qq(hist.edgeLength.data,hist.labels);
lim.x = [0 max(cellfun(@(x) x(end),arrayfun(@(x) xlim(x),hq,'Unif',false)))];
lim.y = [0 max(cellfun(@(x) x(end),arrayfun(@(x) ylim(x),hq,'Unif',false)))];
for i=1:6
    xlim(hq(i),lim.x);
    ylim(hq(i),lim.y);
end
set(h,'PaperPosition',[0 0 3 11]);
print(h,'fig3_edgeLength_qq','-dpdf');

[h , hq , qh] = fig3_qq(hist.edgeLength.data,hist.labels,[0 1]);
set(h,'PaperPosition',[0 0 3 11]);
print(h,'fig3_edgeLength_qq_zoom','-dpdf');
% 
% h = figure;
% for i=1:6
%     hq(i) = subplot(6,1,i);
%     lamins.plot.qqplot_enhanced(hist.edgeLength.data{i*2-1},hist.edgeLength.data{i*2},[0:10:99 99.01:0.01:100]);
%     xlabel(labels{i*2-1});
%     ylabel(labels{i*2});
%     text(0.05,0.9,labels{i*2},'Units','normalized')
%     text(0.9,0.1,labels{i*2-1},'Units','normalized')
% end
% lim.x = [0 max(cellfun(@(x) x(end),arrayfun(@(x) xlim(x),hq,'Unif',false)))];
% lim.y = [0 max(cellfun(@(x) x(end),arrayfun(@(x) ylim(x),hq,'Unif',false)))];
% for i=1:6
%     xlim(hq(i),lim.x);
%     ylim(hq(i),lim.y);
% end
% set(h,'Position',[369 33 290 1064]);
% set(h,'PaperPosition',[0 0 3 11]);
% % set(h,'PaperUnits','inches');
% set(h,'PaperPositionMode','manual');
% print(h,'fig3_edgeLength_qq','-dpdf');

%% Edges per Vertex

h = lamins.functions.dualHistogram(hist.edgesPerVertex.data,hist.labels);
set(h,'PaperPosition',[0 0 3 11]);
print(h,'fig3_edgesPerVertex_hist','-dpdf');

[h , hq , qh] = fig3_qq(hist.edgesPerVertex.data,hist.labels);
set(h,'PaperPosition',[0 0 3 11]);
print(h,'fig3_edgesPerVertex_qq','-dpdf');

[h , hq , qh] = fig3_qq(hist.edgesPerVertex.data,hist.labels,[0 5]);
set(h,'PaperPosition',[0 0 3 11]);
print(h,'fig3_edgesPerVertex_qq_zoom','-dpdf');

% h = figure;
% for i=1:6
%     hq(i) = subplot(6,1,i);
%     lamins.plot.qqplot_enhanced(hist.edgesPerVertex.data{i*2-1},hist.edgesPerVertex.data{i*2},[0:10:99 99.01:0.01:100]);
%     xlabel(labels{i*2-1});
%     ylabel(labels{i*2});
%     text(0.05,0.9,labels{i*2},'Units','normalized')
%     text(0.9,0.1,labels{i*2-1},'Units','normalized')
% end
% lim.x = [0 max(cellfun(@(x) x(end),arrayfun(@(x) xlim(x),hq,'Unif',false)))];
% lim.y = [0 max(cellfun(@(x) x(end),arrayfun(@(x) ylim(x),hq,'Unif',false)))];
% for i=1:6
%     xlim(hq(i),lim.x);
%     ylim(hq(i),lim.y);
% end
% set(h,'Position',[369 33 290 1064]);
% set(h,'PaperPosition',[0 0 3 11]);
% print(h,'fig3_edgesPerVertex_qq','-dpdf');

%% Edges per Face

h = lamins.functions.dualHistogram(hist.edgesPerFace.data,hist.labels);
set(h,'PaperPosition',[0 0 3 11]);
print(h,'fig3_edgesPerFace_hist','-dpdf');

[h , hq , qh] = fig3_qq(hist.edgesPerFace.data,hist.labels);
set(h,'PaperPosition',[0 0 3 11]);
print(h,'fig3_edgesPerFace_qq','-dpdf');

[h , hq , qh] = fig3_qq(hist.edgesPerFace.data,hist.labels,[0 10]);
set(h,'PaperPosition',[0 0 3 11]);
print(h,'fig3_edgesPerFace_qq_zoom','-dpdf');

%% Average Edge Length
[h , hq , qh] = fig3_qq(hist.avgEdgeLength.data,hist.labels);
set(h,'PaperPosition',[0 0 3 11]);
print(h,'fig3_avgEdgeLength_qq','-dpdf');

%% Eccentricity

h = lamins.functions.dualHistogram(hist.eccentricity.data,hist.labels);
set(h,'PaperPosition',[0 0 3 11]);
print(h,'fig3_eccentricity_hist','-dpdf');

% [h , hq , qh] = fig3_qq(hist.eccentricity.data,labels);


    h = figure;
%     pvec = 100 - 0.5.^(1:10);
%     pvec = 100 - [0.5.^(0:10)*10 0];
    pvec = 100 - [0.5.^(0:10)*10 0];
    pvec_low = [0.955.^(0:1:10)*1 0];


    for i=1:6
        hq(i) = subplot(6,1,i);
        qh(i) = qqplot_enhanced(data{i*2+[-1 0]},[0:10:90]);
        hold on;
        scatter_prctile(data{i*2+[-1 0]},pvec_low,'r','Marker','+');
        scatter_prctile(data{i*2+[-1 0]},pvec,'g','Marker','+');
        xlabel(hist.labels{i*2-1});
        ylabel(hist.labels{i*2});
    %     ylim([0 8]);
        text(0.05,0.9,hist.labels{i*2},'Units','normalized')
        text(0.75,0.1,hist.labels{i*2-1},'Units','normalized')
        text(0.2,0.3,['m = ' num2str(qh(i).scaleFactor,3)],'Units','normalized');
        lim = [ 0 1];
            xlim(lim);
            ylim(lim);
    end

set(h,'Position',[369 33 290 1064])
set(h,'PaperPosition',[0 0 3 11]);
print(h,'fig3_eccentricity_qq','-dpdf');

h = figure;
    for i=1:6
        hq(i) = subplot(6,1,i);
        qh(i) = qqplot_enhanced(data{i*2+[-1 0]},[0:10:90]);
        hold on;
        scatter_prctile(data{i*2+[-1 0]},pvec_low,'r','Marker','+');
        scatter_prctile(data{i*2+[-1 0]},pvec,'g','Marker','+');
        xlabel(hist.labels{i*2-1});
        ylabel(hist.labels{i*2});
    %     ylim([0 8]);
        text(0.05,0.9,hist.labels{i*2},'Units','normalized')
        text(0.75,0.1,hist.labels{i*2-1},'Units','normalized')
        text(0.2,0.3,['m = ' num2str(qh(i).scaleFactor,3)],'Units','normalized');
        lim = [ 0 1e-2];
            xlim(lim);
            ylim(lim);
    end

set(h,'Position',[369 33 290 1064])
set(h,'PaperPosition',[0 0 3 11]);
print(h,'fig3_eccentricity_qq_zoom','-dpdf');

% [h , hq , qh] = fig3_qq(hist.eccentricity.data,labels,[0 1]);
% set(h,'PaperPosition',[0 0 3 11]);
% print(h,'fig3_eccentricity_qq_zoom','-dpdf');
% 
% h = figure;
% for i=1:6
%     hq(i) = subplot(6,1,i);
%     lamins.plot.qqplot_enhanced(hist.eccentricity.data{i*2-1},hist.eccentricity.data{i*2},[0:0.1:1 1.1:10:100]);
%     xlabel(labels{i*2-1});
%     ylabel(labels{i*2});
%     text(0.05,0.9,labels{i*2},'Units','normalized')
%     text(0.9,0.1,labels{i*2-1},'Units','normalized')
% end
% lim.x = [0 max(cellfun(@(x) x(end),arrayfun(@(x) xlim(x),hq,'Unif',false)))];
% lim.y = [0 max(cellfun(@(x) x(end),arrayfun(@(x) ylim(x),hq,'Unif',false)))];
% for i=1:6
%     xlim(hq(i),lim.x);
%     ylim(hq(i),lim.y);
% end
% set(h,'Position',[369 33 290 1064]);
% set(h,'PaperPosition',[0 0 3 11]);
% print(h,'fig3_eccentricity_qq','-dpdf');