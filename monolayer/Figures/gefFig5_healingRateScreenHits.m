function [] = gefFig5_healingRateScreenHits()
addpath(genpath('/home2/azaritsky/code/applications/monolayer/Figures/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

dataDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ProjectAll20160516/dayGeneControlFollowup/Speed/';
figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Fig5/';


% load([dataDname '/beta-PIX/Speed_stats.mat']);
% betaPixCntl = controlRepHealingRate;
% betaPixKD = geneRepHealingRate;
% NbetaPix = length(controlMeanPC1);
% nbetaPix = length(betaPixCntl);
% 
% load([dataDname '/SOS1/Speed_stats.mat']);
% sos1Cntl = controlRepHealingRate;
% sos1KD = geneRepHealingRate;
% Nsos1 = length(controlMeanPC1);
% nsos1 = length(sos1Cntl);

load([dataDname '/ARHGEF18/Speed_stats.mat']);
arhgef18Cntl = controlRepHealingRate;
arhgef18KD = geneRepHealingRate;
Narhgef18 = length(controlMeanPC1);
narhgef18 = length(arhgef18Cntl);

load([dataDname '/ARHGEF11/Speed_stats.mat']);
arhgef11Cntl = controlRepHealingRate;
arhgef11KD = geneRepHealingRate;
Narhgef11 = length(controlMeanPC1);
narhgef11 = length(arhgef11Cntl);

load([dataDname '/ARHGEF28/Speed_stats.mat']);
arhgef28Cntl = controlRepHealingRate;
arhgef28KD = geneRepHealingRate;
Narhgef28 = length(controlMeanPC1);
narhgef28 = length(arhgef28Cntl);

load([dataDname '/ARHGEF3/Speed_stats.mat']);
arhgef3Cntl = controlRepHealingRate;
arhgef3KD = geneRepHealingRate;
Narhgef3 = length(controlMeanPC1);
narhgef3 = length(arhgef3Cntl);

%% Stats.
% fprintf(sprintf('betaPix pval (N = %d, n = %d) = %.4f\n',NbetaPix,nbetaPix,signrank(betaPixCntl,betaPixKD)));
% fprintf(sprintf('sos1 pval (N = %d, n = %d) = %.4f\n',Nsos1,nsos1,signrank(sos1Cntl,sos1KD)));
fprintf(sprintf('arhgef18 pval (N = %d, n = %d) = %.4f\n',Narhgef18,narhgef18,signrank(arhgef18Cntl,arhgef18KD)));
fprintf(sprintf('arhgef11 pval (N = %d, n = %d) = %.4f\n',Narhgef11,narhgef11,signrank(arhgef11Cntl,arhgef11KD)));
fprintf(sprintf('arhgef28 pval (N = %d, n = %d) = %.4f\n',Narhgef28,narhgef28,signrank(arhgef28Cntl,arhgef28KD)));
fprintf(sprintf('arhgef3 pval (N = %d, n = %d) = %.4f\n',Narhgef3,narhgef3,signrank(arhgef3Cntl,arhgef3KD)));


%% Figure
% permbetaPix = randperm(nbetaPix);
% permsos1 = randperm(nsos1);
permarhgef18 = randperm(narhgef18);
permarhgef11 = randperm(narhgef11);
permarhgef28 = randperm(narhgef28);
permarhgef3 = randperm(narhgef3);

% diffbetaPix = betaPixCntl - betaPixKD;
% diffsos1 = sos1Cntl - sos1KD;
diffarhgef18 = arhgef18Cntl - arhgef18KD;
diffarhgef11 = arhgef11Cntl - arhgef11KD;
diffarhgef28 = arhgef28Cntl - arhgef28KD;
diffarhgef3 = arhgef3Cntl - arhgef3KD;


% maxDiff = max(abs([diffbetaPix,diffsos1,diffarhgef18,diffarhgef11,diffarhgef28,diffarhgef3]));
maxDiff = max(abs([diffarhgef18,diffarhgef11,diffarhgef28,diffarhgef3]))+3;

% ns = [1 cumsum([nbetaPix,nsos1,narhgef18,narhgef11,narhgef28,narhgef3])];
ns = [1 cumsum([narhgef18,narhgef11,narhgef28,narhgef3])];

% cmap = colormap(hsv(6));
cmap = colormap(hsv(4));

% [left bottom width height]
FPosition = [0 0 800 300]; % FPosition = [0 0 1000 300];
APosition = [0.1 0.2 0.7 0.75]; 

fontsize = 10;

h = figure;
xlabel('Experiment','FontSize',fontsize);
ylabel('Control - KD','FontSize',fontsize);
hold on;
% plot(ns(1):ns(2),diffbetaPix(permbetaPix),'o','MarkerEdgeColor',cmap(1,:),'LineWidth',2,'MarkerSize',8);
% plot((ns(2)+1):ns(3),diffsos1(permsos1),'o','MarkerEdgeColor',cmap(2,:),'LineWidth',2,'MarkerSize',8);
% plot((ns(3)+1):ns(4),diffarhgef18(permarhgef18),'o','MarkerEdgeColor',cmap(3,:),'LineWidth',2,'MarkerSize',8);
% plot((ns(4)+1):ns(5),diffarhgef11(permarhgef11),'o','MarkerEdgeColor',cmap(4,:),'LineWidth',2,'MarkerSize',8);
% plot((ns(5)+1):ns(6),diffarhgef28(permarhgef28),'o','MarkerEdgeColor',cmap(5,:),'LineWidth',2,'MarkerSize',8);
% plot((ns(6)+1):ns(7),diffarhgef3(permarhgef3),'o','MarkerEdgeColor',cmap(6,:),'LineWidth',2,'MarkerSize',8);
% legend('betaPix','sos1','arhgef18','arhgef11','arhgef28','arhgef3','Location','eastoutside');
plot(ns(1):ns(2),diffarhgef18(permarhgef18),'o','MarkerEdgeColor',cmap(1,:),'LineWidth',2,'MarkerSize',8);
plot((ns(2)+1):ns(3),diffarhgef11(permarhgef11),'o','MarkerEdgeColor',cmap(2,:),'LineWidth',2,'MarkerSize',8);
plot((ns(3)+1):ns(4),diffarhgef28(permarhgef28),'o','MarkerEdgeColor',cmap(3,:),'LineWidth',2,'MarkerSize',8);
plot((ns(4)+1):ns(5),diffarhgef3(permarhgef3),'o','MarkerEdgeColor',cmap(4,:),'LineWidth',2,'MarkerSize',8);
legend('arhgef18','arhgef11','arhgef28','arhgef3','Location','eastoutside');
plot([1,ns(5)],[0,0],'--k','LineWidth',2);
haxes = get(h,'CurrentAxes');
% set(haxes,'XLim',[-3,(ns(7)+4)]);
set(haxes,'XLim',[-3,(ns(5)+4)]);
set(haxes,'YLim',[-maxDiff,maxDiff]);
set(haxes,'XTick',0:50:150);
set(haxes,'XTickLabel',0:50:150);
set(haxes,'YTick',[-15,0,15]);
set(haxes,'YTickLabel',[-15,0,15]);
set(haxes,'FontSize',fontsize);
set(h,'Color','w');
set(h,'Position',FPosition,'PaperPositionMode','auto');
axisHandle= findobj(h,'type','axes');
set(axisHandle,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
set(get(axisHandle,'XLabel'),'FontSize',fontsize); set(get(axisHandle,'YLabel'),'FontSize',fontsize);
export_fig([figDname 'healingRateScreenHits_diff_legend.eps']);
legend off;

FPosition = [0 0 600 300];
APosition = [0.1 0.2 0.85 0.75]; 
set(h,'Position',FPosition,'PaperPositionMode','auto');
axisHandle= findobj(h,'type','axes');
set(axisHandle,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
set(get(axisHandle,'XLabel'),'FontSize',fontsize); set(get(axisHandle,'YLabel'),'FontSize',fontsize);
export_fig([figDname 'healingRateScreenHits_diff.eps']);
hold off;

end