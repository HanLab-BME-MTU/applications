function [] = gefFig1_healingRateRhoGTPAses()
addpath(genpath('/home2/azaritsky/code/applications/monolayer/Figures/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

dataDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ProjectAll20160516/dayGeneControlFollowup/Speed/';
figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Fig1/';


load([dataDname '/RHOA/Speed_stats.mat']);
rhoaCntl = controlRepHealingRate;
rhoaKD = geneRepHealingRate;
NRHOA = length(controlMeanPC1);
nRHOA = length(rhoaCntl);

load([dataDname '/CDC42/Speed_stats.mat']);
cdc42Cntl = controlRepHealingRate;
cdc42KD = geneRepHealingRate;
NCDC42 = length(controlMeanPC1);
nCDC42 = length(cdc42Cntl);

load([dataDname '/RAC1/Speed_stats.mat']);
rac1Cntl = controlRepHealingRate;
rac1KD = geneRepHealingRate;
NRAC1 = length(controlMeanPC1);
nRAC1 = length(rac1Cntl);

%% Stats.
fprintf(sprintf('CDC42 pval (N = %d, n = %d) = %.4f\n',NCDC42,nCDC42,signrank(cdc42Cntl,cdc42KD)));
fprintf(sprintf('RAC1 pval (N = %d, n = %d) = %.4f\n',NRAC1,nRAC1,signrank(rac1Cntl,rac1KD)));
fprintf(sprintf('RHOA pval (N = %d, n = %d) = %.4f\n',NRHOA,nRHOA,signrank(rhoaCntl,rhoaKD)));

%% Figure
minWH = min([cdc42Cntl,cdc42KD,rac1Cntl,rac1KD,rhoaCntl,rhoaKD]);
maxWH = max([cdc42Cntl,cdc42KD,rac1Cntl,rac1KD,rhoaCntl,rhoaKD]);

cmap = colormap(hsv(3));

fontsize = 16;
h = figure;
xlabel('Control ','FontSize',fontsize);
ylabel('Condition','FontSize',fontsize);
hold on;
plot(cdc42Cntl,cdc42KD,'o','MarkerEdgeColor',cmap(1,:),'LineWidth',2,'MarkerSize',8);
plot(rac1Cntl,rac1KD,'o','MarkerEdgeColor',cmap(2,:),'LineWidth',2,'MarkerSize',8);
plot(rhoaCntl,rhoaKD,'o','MarkerEdgeColor',cmap(3,:),'LineWidth',2,'MarkerSize',8);
legend('CDC42','RAC1','RHOA','Location','northwest');
plot([minWH,maxWH],[minWH,maxWH],'--k','LineWidth',2);
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[minWH,maxWH]);
set(haxes,'YLim',[minWH,maxWH]);
set(haxes,'FontSize',fontsize);
set(h,'Color','w');
axis square;
export_fig([figDname 'healingRateRhoGTPAses_cntlVsKD.eps']);
hold off;

%%

permCDC42 = randperm(nCDC42);
permRAC1 = randperm(nRAC1);
permRHOA = randperm(nRHOA);

diffCDC42 = cdc42KD - cdc42Cntl;
diffRAC1 = rac1KD - rac1Cntl;
diffRHOA = rhoaKD - rhoaCntl;

maxDiff = max(abs([diffCDC42,diffRAC1,diffRHOA]));

% [left bottom width height]
FPosition = [0 0 500 300];
APosition = [0.1 0.2 0.85 0.75]; 

fontsize = 10;

h = figure;
xlabel('Experiment','FontSize',fontsize);
ylabel('KD - Control (\mum hour^{-1})','FontSize',fontsize);
hold on;
plot(1:nCDC42,diffCDC42(permCDC42),'o','MarkerEdgeColor',cmap(1,:),'LineWidth',2,'MarkerSize',8);
plot((nCDC42+1):(nCDC42+nRAC1),diffRAC1(permRAC1),'o','MarkerEdgeColor',cmap(2,:),'LineWidth',2,'MarkerSize',8);
plot((nCDC42+nRAC1+1):(nCDC42+nRAC1+nRHOA),diffRHOA(permRHOA),'o','MarkerEdgeColor',cmap(3,:),'LineWidth',2,'MarkerSize',8);
% legend('CDC42','RAC1','RHOA','Location','southwest');
plot([1,nCDC42+nRAC1+nRHOA],[0,0],'--k','LineWidth',2);
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[-3,(nCDC42+nRAC1+nRHOA+4)]);
set(haxes,'YLim',[-maxDiff,maxDiff]);
set(haxes,'XTick',[0,50,100]);
set(haxes,'XTickLabel',[0,50,100]);
set(haxes,'YTick',[-30,0,30]);
set(haxes,'YTickLabel',[-30,0,30]);
set(haxes,'FontSize',fontsize);
set(h,'Color','w');
set(h,'Position',FPosition,'PaperPositionMode','auto');
axisHandle= findobj(h,'type','axes');
set(axisHandle,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
set(get(axisHandle,'XLabel'),'FontSize',fontsize); set(get(axisHandle,'YLabel'),'FontSize',fontsize);
export_fig([figDname 'healingRateRhoGTPAses_diff.eps']);
hold off;

end