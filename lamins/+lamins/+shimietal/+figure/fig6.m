%% Preliminaries

import lamins.plot.*;

%% Area
quantiles = 10:10:90;
% pvec = 100 - [(sqrt(2).^-(0:50))*100 0];
% pvec = 100 - [(sqrt(2).^-(0:50))*10 0];
pvec = 100 - [0.5.^(1:10)*10 0];



h = figure;
subplot(4,3,1);
qh(1) = lamins.plot.qqplot_enhanced(stats(20).area(2).all.data,stats(2).area(1).all.data,quantiles);
hold on;
scatter_prctile(stats(20).area(2).all.data,stats(2).area(1).all.data,pvec,'g','Marker','^');
ylabel('LAC null Area');
xlabel('WT Area');
lim = ylim;
lim(1) = 0;
ylim(lim);
text(0.4,0.3,['m = ' num2str(qh(1).scaleFactor,3)],'Units','normalized');
        text(0.05,0.9,'Lmna-/-','Units','normalized','FontAngle','italic')
        text(0.85,0.1,'wt','Units','normalized')

subplot(4,3,2);
qh(2) = lamins.plot.qqplot_enhanced(stats(20).area(1).all.data,stats(7).area(1).all.data,quantiles)
hold on;
scatter_prctile(stats(20).area(1).all.data,stats(7).area(1).all.data,pvec,'g','Marker','^');
ylabel('LB1 null Area');
xlabel('WT  Area');
lim = ylim;
lim(1) = 0;
ylim(lim);
text(0.4,0.3,['m = ' num2str(qh(2).scaleFactor,3)],'Units','normalized');
        text(0.05,0.9,'Lmnb1-/-','Units','normalized','FontAngle','italic')
        text(0.85,0.1,'wt','Units','normalized')


subplot(4,3,3);
% movies 3 and 7 did not acquire a good focal plane
goodLB2null = [stats(12).area(1).movie([1 2 4 5 6 8 9 10]).data];
qh(3) = lamins.plot.qqplot_enhanced(stats(20).area(1).all.data,goodLB2null,quantiles)
hold on;
scatter_prctile(stats(20).area(1).all.data,goodLB2null,pvec,'g','Marker','^');
ylabel('LB2 null Area');
xlabel('WT Area');
lim = ylim;
lim(1) = 0;
ylim(lim);
text(0.4,0.3,['m = ' num2str(qh(3).scaleFactor,3)],'Units','normalized');
        text(0.05,0.9,'Lmnb2-/-','Units','normalized','FontAngle','italic')
        text(0.85,0.1,'wt','Units','normalized')

%% average edge length

quantiles = 10:10:90;
pvec = 100 - [0.5.^(1:10)*10 0];
pvec_low = [0.955.^(0:1:10)*1 0];

for i=1:length(stats)
    for j=1:length(stats(i).perimeter)
        stats(i).avgEdgeLength(j).all.data = stats(i).perimeter(j).all.data ./ stats(i).edgesPerFace(j).all.data;
        stats(i).circularity(j).all.data = 4*pi*stats(i).area(j).all.data ./ stats(i).perimeter(j).all.data.^2;
        for k=1:length(stats(i).perimeter(j).movie)
            stats(i).avgEdgeLength(j).movie(k).data = stats(i).perimeter(j).movie(k).data ./ stats(i).edgesPerFace(j).movie(k).data;
            stats(i).circularity(j).movie(k).data = 4*pi*stats(i).area(j).movie(k).data ./ stats(i).perimeter(j).movie(k).data.^2;
        end
    end
end

subplot(4,3,4);
qh(4) = lamins.plot.qqplot_enhanced(stats(20).avgEdgeLength(2).all.data, stats(2).avgEdgeLength(1).all.data,quantiles);
hold on;
% scatter_prctile(stats(20).avgEdgeLength(2).all.data,stats(2).avgEdgeLength(1).all.data,pvec_low,'r','Marker','+');
scatter_prctile(stats(20).avgEdgeLength(2).all.data,stats(2).avgEdgeLength(1).all.data,pvec,'g','Marker','^');
ylabel('LAC null avgEdgeLength');
xlabel('WT avgEdgeLength');
lim = ylim;
lim(1) = 0;
% ylim([0 3]);
text(0.4,0.3,['m = ' num2str(qh(4).scaleFactor,3)],'Units','normalized');
        text(0.05,0.9,'Lmna-/-','Units','normalized','FontAngle','italic')
        text(0.85,0.1,'wt','Units','normalized')

subplot(4,3,5);
qh(5) = lamins.plot.qqplot_enhanced(stats(20).avgEdgeLength(1).all.data,stats(7).avgEdgeLength(1).all.data,quantiles);
hold on;
% scatter_prctile(stats(20).avgEdgeLength(1).all.data,stats(7).avgEdgeLength(1).all.data,pvec_low,'r','Marker','+');
scatter_prctile(stats(20).avgEdgeLength(1).all.data,stats(7).avgEdgeLength(1).all.data,pvec,'g','Marker','^');
ylabel('LB1 null avgEdgeLength');
xlabel('WT  avgEdgeLength');
lim = ylim;
lim(1) = 0;
% ylim([0 3]);
text(0.4,0.3,['m = ' num2str(qh(5).scaleFactor,3)],'Units','normalized');
        text(0.05,0.9,'Lmnb1-/-','Units','normalized','FontAngle','italic')
        text(0.85,0.1,'wt','Units','normalized')


subplot(4,3,6);
% movies 3 and 7 did not acquire a good focal plane
goodLB2null = [stats(12).avgEdgeLength(1).movie([1 2 4 5 6 8 9 10]).data];
qh(6) = lamins.plot.qqplot_enhanced(stats(20).avgEdgeLength(1).all.data,goodLB2null,quantiles);
hold on;
% scatter_prctile(stats(20).avgEdgeLength(1).all.data,goodLB2null,pvec_low,'r','Marker','+');
scatter_prctile(stats(20).avgEdgeLength(1).all.data,goodLB2null,pvec,'g','Marker','^');
ylabel('LB2 null avgEdgeLength');
xlabel('WT avgEdgeLength');
lim = ylim;
lim(1) = 0;
% ylim([0 3]);
text(0.4,0.3,['m = ' num2str(qh(6).scaleFactor,3)],'Units','normalized');
        text(0.05,0.9,'Lmnb2-/-','Units','normalized','FontAngle','italic')
        text(0.85,0.1,'wt','Units','normalized')

%% Edges Per Face

subplot(4,3,7);
qh(7) = lamins.plot.qqplot_enhanced(stats(20).edgesPerFace(2).all.data,stats(2).edgesPerFace(1).all.data,quantiles)
hold on;
% scatter_prctile(stats(20).edgesPerFace(2).all.data,stats(2).edgesPerFace(1).all.data,pvec_low,'r','Marker','+');
scatter_prctile(stats(20).edgesPerFace(2).all.data,stats(2).edgesPerFace(1).all.data,pvec,'g','Marker','^');
ylabel('LAC null edgesPerFace');
xlabel('WT edgesPerFace');
lim = ylim;
lim(1) = 0;
ylim([0 30]);
xlim([0 30]);
text(0.4,0.3,['m = ' num2str(qh(7).scaleFactor,3)],'Units','normalized');
        text(0.05,0.9,'Lmna-/-','Units','normalized','FontAngle','italic')
        text(0.85,0.1,'wt','Units','normalized')

subplot(4,3,8);
qh(8) = lamins.plot.qqplot_enhanced(stats(20).edgesPerFace(1).all.data,stats(7).edgesPerFace(1).all.data,quantiles)
hold on;
% scatter_prctile(stats(20).edgesPerFace(1).all.data,stats(7).edgesPerFace(1).all.data,pvec_low,'r','Marker','+');
scatter_prctile(stats(20).edgesPerFace(1).all.data,stats(7).edgesPerFace(1).all.data,pvec,'g','Marker','^');
ylabel('LB1 null edgesPerFace');
xlabel('WT edgesPerFace');
lim = ylim;
lim(1) = 0;
ylim(lim);
ylim([0 30]);
xlim([0 30]);
text(0.4,0.3,['m = ' num2str(qh(8).scaleFactor,3)],'Units','normalized');
        text(0.05,0.9,'Lmnb1-/-','Units','normalized','FontAngle','italic')
        text(0.85,0.1,'wt','Units','normalized')



subplot(4,3,9);
% movies 3 and 7 did not acquire a good focal plane
goodLB2null = [stats(12).edgesPerFace(1).movie([1 2 4 5 6 8 9 10]).data];
qh(9) = lamins.plot.qqplot_enhanced(stats(20).edgesPerFace(1).all.data,goodLB2null,quantiles)
hold on;
% scatter_prctile(stats(20).edgesPerFace(1).all.data,goodLB2null,pvec_low,'r','Marker','+');
scatter_prctile(stats(20).edgesPerFace(1).all.data,goodLB2null,pvec,'g','Marker','^');
ylabel('LB2 null edgesPerFace');
xlabel('WT edgesPerFace');
lim = ylim;
lim(1) = 0;
ylim([0 30]);
xlim([0 30]);
text(0.4,0.3,['m = ' num2str(qh(9).scaleFactor,3)],'Units','normalized');
        text(0.05,0.9,'Lmnb2-/-','Units','normalized','FontAngle','italic')
        text(0.85,0.1,'wt','Units','normalized')

%% Eccentricity
% figure;

% pvec = 1.0471.^(-100:100);
pvec = [0.955.^(0:1:10)*1 0];

subplot(4,3,10);
qh(10) = lamins.plot.qqplot_enhanced(stats(20).eccentricity(2).all.data,stats(2).eccentricity(1).all.data,quantiles)
hold on;
scatter_prctile(stats(20).eccentricity(2).all.data,stats(2).eccentricity(1).all.data,pvec,'r','Marker','s');
ylabel('LAC null Eccentricity');
xlabel('WT Eccentricity');
ylim([0 1]);
xlim([0 1]);
text(0.35,0.3,['m = ' num2str(qh(10).scaleFactor,3)],'Units','normalized');
        text(0.05,0.9,'Lmna-/-','Units','normalized','FontAngle','italic')
        text(0.85,0.1,'wt','Units','normalized')

subplot(4,3,11);
qh(11) = lamins.plot.qqplot_enhanced(stats(20).eccentricity(1).all.data,stats(7).eccentricity(1).all.data,quantiles)
hold on;
scatter_prctile(stats(20).eccentricity(1).all.data,stats(7).eccentricity(1).all.data,pvec,'r','Marker','s');
ylabel('LB1 null Eccentricity');
xlabel('WT Eccentricity');
ylim([0 1]);
xlim([0 1]);
text(0.35,0.3,['m = ' num2str(qh(11).scaleFactor,3)],'Units','normalized');
        text(0.05,0.9,'Lmnb1-/-','Units','normalized','FontAngle','italic')
        text(0.85,0.1,'wt','Units','normalized')



subplot(4,3,12);
% movies 3 and 7 did not acquire a good focal plane
goodLB2null = [stats(12).eccentricity(1).movie([1 2 4 5 6 8 9 10]).data];
qh(12) = lamins.plot.qqplot_enhanced(stats(20).eccentricity(1).all.data,goodLB2null,quantiles)
hold on;
scatter_prctile(stats(20).eccentricity(1).all.data,goodLB2null,pvec,'r','Marker','s');
ylabel('LB2 null Eccentricity');
xlabel('WT Eccentricity');
ylim([0 1]);
xlim([0 1]);
text(0.35,0.3,['m = ' num2str(qh(12).scaleFactor,3)],'Units','normalized');
        text(0.05,0.9,'Lmnb2-/-','Units','normalized','FontAngle','italic')
        text(0.85,0.1,'wt','Units','normalized')

%% Save figure
set(h,'Position',[680 39 809 1058]);
set(h,'PaperOrientation','portrait');
set(h,'PaperPositionMode','auto');
print(h,'fig6','-dpdf');

%%
pause;

%% Edge length

quantiles = 10:10:90;
pvec = 100 - [0.5.^(0:10)*10 0];
pvec_low = [0.955.^(0:1:10)*1 0];

h = figure;
subplot(2,3,1);
qh(1) = lamins.plot.qqplot_enhanced(stats(20).edgeLength(2).all.data,stats(2).edgeLength(1).all.data,quantiles);
hold on;
scatter_prctile(stats(20).edgeLength(2).all.data,stats(2).edgeLength(1).all.data,pvec_low,'r','Marker','+');
scatter_prctile(stats(20).edgeLength(2).all.data,stats(2).edgeLength(1).all.data,pvec,'g','Marker','+');
ylabel('LAC null edgeLength');
xlabel('WT edgeLength');
lim = ylim;
lim(1) = 0;
ylim([0 3]);
text(0.2,0.3,['m = ' num2str(qh(1).scaleFactor,3)],'Units','normalized');

subplot(2,3,2);
qh(2) = lamins.plot.qqplot_enhanced(stats(20).edgeLength(1).all.data,stats(7).edgeLength(1).all.data,quantiles);
hold on;
scatter_prctile(stats(20).edgeLength(1).all.data,stats(7).edgeLength(1).all.data,pvec_low,'r','Marker','+');
scatter_prctile(stats(20).edgeLength(1).all.data,stats(7).edgeLength(1).all.data,pvec,'g','Marker','+');
ylabel('LB1 null edgeLength');
xlabel('WT  edgeLength');
lim = ylim;
lim(1) = 0;
ylim([0 3]);
text(0.2,0.3,['m = ' num2str(qh(2).scaleFactor,3)],'Units','normalized');


subplot(2,3,3);
% movies 3 and 7 did not acquire a good focal plane
goodLB2null = [stats(12).edgeLength(1).movie([1 2 4 5 6 8 9 10]).data];
qh(3) = lamins.plot.qqplot_enhanced(stats(20).edgeLength(1).all.data,goodLB2null,quantiles);
hold on;
scatter_prctile(stats(20).edgeLength(1).all.data,goodLB2null,pvec_low,'r','Marker','+');
scatter_prctile(stats(20).edgeLength(1).all.data,goodLB2null,pvec,'g','Marker','+');
ylabel('LB2 null edgeLength');
xlabel('WT edgeLength');
lim = ylim;
lim(1) = 0;
ylim([0 3]);
text(0.2,0.3,['m = ' num2str(qh(3).scaleFactor,3)],'Units','normalized');

%% Edges Per Face

subplot(2,3,4);
qh(4) = lamins.plot.qqplot_enhanced(stats(20).edgesPerFace(2).all.data,stats(2).edgesPerFace(1).all.data,quantiles)
hold on;
scatter_prctile(stats(20).edgesPerFace(2).all.data,stats(2).edgesPerFace(1).all.data,pvec_low,'r','Marker','+');
scatter_prctile(stats(20).edgesPerFace(2).all.data,stats(2).edgesPerFace(1).all.data,pvec,'g','Marker','+');
ylabel('LAC null edgesPerFace');
xlabel('WT edgesPerFace');
lim = ylim;
lim(1) = 0;
ylim([0 40]);
xlim([0 40]);
text(0.2,0.3,['m = ' num2str(qh(4).scaleFactor,3)],'Units','normalized');

subplot(2,3,5);
qh(5) = lamins.plot.qqplot_enhanced(stats(20).edgesPerFace(1).all.data,stats(7).edgesPerFace(1).all.data,quantiles)
hold on;
scatter_prctile(stats(20).edgesPerFace(1).all.data,stats(7).edgesPerFace(1).all.data,pvec_low,'r','Marker','+');
scatter_prctile(stats(20).edgesPerFace(1).all.data,stats(7).edgesPerFace(1).all.data,pvec,'g','Marker','+');
ylabel('LB1 null edgesPerFace');
xlabel('WT edgesPerFace');
lim = ylim;
lim(1) = 0;
ylim(lim);
ylim([0 40]);
xlim([0 40]);
text(0.2,0.3,['m = ' num2str(qh(5).scaleFactor,3)],'Units','normalized');



subplot(2,3,6);
% movies 3 and 7 did not acquire a good focal plane
goodLB2null = [stats(12).edgesPerFace(1).movie([1 2 4 5 6 8 9 10]).data];
qh(6) = lamins.plot.qqplot_enhanced(stats(20).edgesPerFace(1).all.data,goodLB2null,quantiles)
hold on;
scatter_prctile(stats(20).edgesPerFace(1).all.data,goodLB2null,pvec_low,'r','Marker','+');
scatter_prctile(stats(20).edgesPerFace(1).all.data,goodLB2null,pvec,'g','Marker','+');
ylabel('LB2 null edgesPerFace');
xlabel('WT edgesPerFace');
lim = ylim;
lim(1) = 0;
ylim([0 40]);
xlim([0 40]);
text(0.2,0.3,['m = ' num2str(qh(6).scaleFactor,3)],'Units','normalized');

%% Second Panel
set(h,'Position',[680 667 781 430]);
print(h,'fig6_2','-dpdf');

%% perimeter

h = figure;
subplot(3,3,1);
qh(1) = lamins.plot.qqplot_enhanced(stats(20).perimeter(2).all.data,stats(2).perimeter(1).all.data,quantiles);
hold on;
scatter_prctile(stats(20).perimeter(2).all.data,stats(2).perimeter(1).all.data,pvec_low,'r','Marker','+');
scatter_prctile(stats(20).perimeter(2).all.data,stats(2).perimeter(1).all.data,pvec,'g','Marker','+');
ylabel('LAC null perimeter');
xlabel('WT perimeter');
lim = ylim;
lim(1) = 0;
% ylim([0 3]);
text(0.2,0.3,['m = ' num2str(qh(1).scaleFactor,3)],'Units','normalized');

subplot(3,3,2);
qh(2) = lamins.plot.qqplot_enhanced(stats(20).perimeter(1).all.data,stats(7).perimeter(1).all.data,quantiles);
hold on;
scatter_prctile(stats(20).perimeter(1).all.data,stats(7).perimeter(1).all.data,pvec_low,'r','Marker','+');
scatter_prctile(stats(20).perimeter(1).all.data,stats(7).perimeter(1).all.data,pvec,'g','Marker','+');
ylabel('LB1 null perimeter');
xlabel('WT  perimeter');
lim = ylim;
lim(1) = 0;
% ylim([0 3]);
text(0.2,0.3,['m = ' num2str(qh(2).scaleFactor,3)],'Units','normalized');


subplot(3,3,3);
% movies 3 and 7 did not acquire a good focal plane
goodLB2null = [stats(12).perimeter(1).movie([1 2 4 5 6 8 9 10]).data];
qh(3) = lamins.plot.qqplot_enhanced(stats(20).perimeter(1).all.data,goodLB2null,quantiles);
hold on;
scatter_prctile(stats(20).perimeter(1).all.data,goodLB2null,pvec_low,'r','Marker','+');
scatter_prctile(stats(20).perimeter(1).all.data,goodLB2null,pvec,'g','Marker','+');
ylabel('LB2 null perimeter');
xlabel('WT perimeter');
lim = ylim;
lim(1) = 0;
% ylim([0 3]);
text(0.2,0.3,['m = ' num2str(qh(3).scaleFactor,3)],'Units','normalized');

%% average edge length

for i=1:length(stats)
    for j=1:length(stats(i).perimeter)
        stats(i).avgEdgeLength(j).all.data = stats(i).perimeter(j).all.data ./ stats(i).edgesPerFace(j).all.data;
        stats(i).circularity(j).all.data = 4*pi*stats(i).area(j).all.data ./ stats(i).perimeter(j).all.data.^2;
        for k=1:length(stats(i).perimeter(j).movie)
            stats(i).circularity(j).movie(k).data = 4*pi*stats(i).area(j).movie(k).data ./ stats(i).perimeter(j).movie(k).data.^2;
        end
    end
end

subplot(3,3,4);
qh(4) = lamins.plot.qqplot_enhanced(stats(20).avgEdgeLength(2).all.data, stats(2).avgEdgeLength(1).all.data,quantiles);
hold on;
scatter_prctile(stats(20).avgEdgeLength(2).all.data,stats(2).avgEdgeLength(1).all.data,pvec_low,'r','Marker','+');
scatter_prctile(stats(20).avgEdgeLength(2).all.data,stats(2).avgEdgeLength(1).all.data,pvec,'g','Marker','+');
ylabel('LAC null avgEdgeLength');
xlabel('WT avgEdgeLength');
lim = ylim;
lim(1) = 0;
% ylim([0 3]);
text(0.2,0.3,['m = ' num2str(qh(4).scaleFactor,3)],'Units','normalized');

subplot(3,3,5);
qh(5) = lamins.plot.qqplot_enhanced(stats(20).avgEdgeLength(1).all.data,stats(7).avgEdgeLength(1).all.data,quantiles);
hold on;
scatter_prctile(stats(20).avgEdgeLength(1).all.data,stats(7).avgEdgeLength(1).all.data,pvec_low,'r','Marker','+');
scatter_prctile(stats(20).avgEdgeLength(1).all.data,stats(7).avgEdgeLength(1).all.data,pvec,'g','Marker','+');
ylabel('LB1 null avgEdgeLength');
xlabel('WT  avgEdgeLength');
lim = ylim;
lim(1) = 0;
% ylim([0 3]);
text(0.2,0.3,['m = ' num2str(qh(5).scaleFactor,3)],'Units','normalized');


subplot(3,3,6);
% movies 3 and 7 did not acquire a good focal plane
goodLB2null = [stats(12).avgEdgeLength(1).movie([1 2 4 5 6 8 9 10]).data];
qh(6) = lamins.plot.qqplot_enhanced(stats(20).avgEdgeLength(1).all.data,goodLB2null,quantiles);
hold on;
scatter_prctile(stats(20).avgEdgeLength(1).all.data,goodLB2null,pvec_low,'r','Marker','+');
scatter_prctile(stats(20).avgEdgeLength(1).all.data,goodLB2null,pvec,'g','Marker','+');
ylabel('LB2 null avgEdgeLength');
xlabel('WT avgEdgeLength');
lim = ylim;
lim(1) = 0;
% ylim([0 3]);
text(0.2,0.3,['m = ' num2str(qh(6).scaleFactor,3)],'Units','normalized');

%% Edges Per Face

subplot(3,3,7);
qh(4) = lamins.plot.qqplot_enhanced(stats(20).edgesPerFace(2).all.data,stats(2).edgesPerFace(1).all.data,quantiles)
hold on;
scatter_prctile(stats(20).edgesPerFace(2).all.data,stats(2).edgesPerFace(1).all.data,pvec_low,'r','Marker','+');
scatter_prctile(stats(20).edgesPerFace(2).all.data,stats(2).edgesPerFace(1).all.data,pvec,'g','Marker','+');
ylabel('LAC null edgesPerFace');
xlabel('WT edgesPerFace');
lim = ylim;
lim(1) = 0;
ylim([0 40]);
xlim([0 40]);
text(0.2,0.3,['m = ' num2str(qh(4).scaleFactor,3)],'Units','normalized');

subplot(3,3,8);
qh(5) = lamins.plot.qqplot_enhanced(stats(20).edgesPerFace(1).all.data,stats(7).edgesPerFace(1).all.data,quantiles)
hold on;
scatter_prctile(stats(20).edgesPerFace(1).all.data,stats(7).edgesPerFace(1).all.data,pvec_low,'r','Marker','+');
scatter_prctile(stats(20).edgesPerFace(1).all.data,stats(7).edgesPerFace(1).all.data,pvec,'g','Marker','+');
ylabel('LB1 null edgesPerFace');
xlabel('WT edgesPerFace');
lim = ylim;
lim(1) = 0;
ylim(lim);
ylim([0 40]);
xlim([0 40]);
text(0.2,0.3,['m = ' num2str(qh(5).scaleFactor,3)],'Units','normalized');



subplot(3,3,9);
% movies 3 and 7 did not acquire a good focal plane
goodLB2null = [stats(12).edgesPerFace(1).movie([1 2 4 5 6 8 9 10]).data];
qh(6) = lamins.plot.qqplot_enhanced(stats(20).edgesPerFace(1).all.data,goodLB2null,quantiles)
hold on;
scatter_prctile(stats(20).edgesPerFace(1).all.data,goodLB2null,pvec_low,'r','Marker','+');
scatter_prctile(stats(20).edgesPerFace(1).all.data,goodLB2null,pvec,'g','Marker','+');
ylabel('LB2 null edgesPerFace');
xlabel('WT edgesPerFace');
lim = ylim;
lim(1) = 0;
ylim([0 40]);
xlim([0 40]);
text(0.2,0.3,['m = ' num2str(qh(6).scaleFactor,3)],'Units','normalized');

%% Third Panel
set(h,'Position',[680 667 781 645]);
print(h,'fig6_3','-dpdf');


%% Circularity
figure;
D = dir('/project/biophysics/jaqaman_lab/lamins/2015/20150602/MEF*');
circularity = arrayfun(@(s) nanmean(vertcat(s.nucleusCircularity.maskCircularity)),stats)';
circularityTable = cell2table([{D.name}' num2cell(circularity)],'VariableNames',{'Set','Circularity'});
circularityTable.NinetyNine_FaceArea = arrayfun(@(s) prctile(s.area(1).all.data,99),stats)';
scatter(circularityTable.Circularity,circularityTable.NinetyNine_FaceArea,'k','filled')
text(circularityTable.Circularity,circularityTable.NinetyNine_FaceArea,circularityTable.Set)
ylabel('99th Pernctile of Meshwork Face Area');
xlabel('Nucleus Circularity (Shape Factor)');
corr(circularityTable.Circularity,circularityTable.NinetyNine_FaceArea)