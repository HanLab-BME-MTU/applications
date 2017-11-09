function [] = gefFig2_PCsRhoGTPAses_directionality()
addpath(genpath('/home2/azaritsky/code/applications/monolayer/Figures/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

dataDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ProjectAll20160516/dayGeneControlFollowup/Directionality/';
figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Fig2/RhoGTPAses/directionality/';


load([dataDname '/RHOA/Directionality_stats.mat']);
rhoaCntlPC1 = controlRepPC1;
rhoaKDPC1 = geneRepPC1;
rhoaCntlPC2 = controlRepPC2;
rhoaKDPC2 = geneRepPC2;
rhoaCntlPC3 = controlRepPC3;
rhoaKDPC3 = geneRepPC3;
NRHOA = length(controlMeanPC1);
nRHOA = length(controlRepPC1);

load([dataDname '/CDC42/Directionality_stats.mat']);
cdc42CntlPC1 = controlRepPC1;
cdc42KDPC1 = geneRepPC1;
cdc42CntlPC2 = controlRepPC2;
cdc42KDPC2 = geneRepPC2;
cdc42CntlPC3 = controlRepPC3;
cdc42KDPC3 = geneRepPC3;
NCDC42 = length(controlMeanPC1);
nCDC42 = length(controlRepPC1);

load([dataDname '/RAC1/Directionality_stats.mat']);
rac1CntlPC1 = controlRepPC1;
rac1KDPC1 = geneRepPC1;
rac1CntlPC2 = controlRepPC2;
rac1KDPC2 = geneRepPC2;
rac1CntlPC3 = controlRepPC3;
rac1KDPC3 = geneRepPC3;
NRAC1 = length(controlMeanPC1);
nRAC1 = length(controlRepPC1);

%% Stats.
fprintf(sprintf('CDC42 (N = %d, n = %d)\n',NCDC42,nCDC42));
fprintf(sprintf('PC1 pval = %.4f\n',signrank(cdc42CntlPC1,cdc42KDPC1)));
fprintf(sprintf('PC2 pval = %.4f\n',signrank(cdc42CntlPC2,cdc42KDPC2)));
fprintf(sprintf('PC3 pval = %.4f\n\n',signrank(cdc42CntlPC3,cdc42KDPC3)));

fprintf(sprintf('RAC1 (N = %d, n = %d)\n',NRAC1,nRAC1));
fprintf(sprintf('PC1 pval = %.4f\n',signrank(rac1CntlPC1,rac1KDPC1)));
fprintf(sprintf('PC2 pval = %.4f\n',signrank(rac1CntlPC2,rac1KDPC2)));
fprintf(sprintf('PC3 pval = %.4f\n\n',signrank(rac1CntlPC3,rac1KDPC3)));

fprintf(sprintf('RHOA (N = %d, n = %d)\n',NRHOA,nRHOA));
fprintf(sprintf('PC1 pval = %.4f\n',signrank(rhoaCntlPC1,rhoaKDPC1)));
fprintf(sprintf('PC2 pval = %.4f\n',signrank(rhoaCntlPC2,rhoaKDPC2)));
fprintf(sprintf('PC3 pval = %.4f\n',signrank(rhoaCntlPC3,rhoaKDPC3)));

%% Figure
permCDC42 = randperm(nCDC42);
permRAC1 = randperm(nRAC1);
permRHOA = randperm(nRHOA);

diffCDC42PC1 = -(cdc42CntlPC1 - cdc42KDPC1);
diffRAC1PC1 = -(rac1CntlPC1 - rac1KDPC1);
diffRHOAPC1 = -(rhoaCntlPC1 - rhoaKDPC1);

diffCDC42PC2 = -(cdc42CntlPC2 - cdc42KDPC2);
diffRAC1PC2 = -(rac1CntlPC2 - rac1KDPC2);
diffRHOAPC2 = -(rhoaCntlPC2 - rhoaKDPC2);

diffCDC42PC3 = -(cdc42CntlPC3 - cdc42KDPC3);
diffRAC1PC3 = -(rac1CntlPC3 - rac1KDPC3);
diffRHOAPC3 = -(rhoaCntlPC3 - rhoaKDPC3);

maxDiffPC1 = max(abs([diffCDC42PC1,diffRAC1PC1,diffRHOAPC1]));
maxDiffPC2 = max(abs([diffCDC42PC2,diffRAC1PC2,diffRHOAPC2]));
maxDiffPC3 = max(abs([diffCDC42PC3,diffRAC1PC3,diffRHOAPC3]));
    
plotPCsRhoGTPAses(diffCDC42PC1(permCDC42),diffRAC1PC1(permRAC1),diffRHOAPC1(permRHOA),maxDiffPC1,'PC1',[figDname 'PC1_RhoGTPAses_cntlVsKD.eps']);
plotPCsRhoGTPAses(diffCDC42PC2(permCDC42),diffRAC1PC2(permRAC1),diffRHOAPC2(permRHOA),maxDiffPC2,'PC2',[figDname 'PC2_RhoGTPAses_cntlVsKD.eps']);
plotPCsRhoGTPAses(diffCDC42PC3(permCDC42),diffRAC1PC3(permRAC1),diffRHOAPC3(permRHOA),maxDiffPC3,'PC3',[figDname 'PC3_RhoGTPAses_cntlVsKD.eps']);
end

%%
function [] = plotPCsRhoGTPAses(diffCDC42,diffRAC1,diffRHOA,maxDiff,pcStr,outFname)

nCDC42 = length(diffCDC42);
nRAC1 = length(diffRAC1);
nRHOA = length(diffRHOA);

cmap = colormap(hsv(3));

% [left bottom width height]
FPosition = [0 0 500 300];
APosition = [0.1 0.2 0.85 0.75]; 

fontsize = 10;

h = figure;
xlabel('Experiment','FontSize',fontsize);
ylabel([pcStr ': KD - Control'],'FontSize',fontsize);
hold on;
plot(1:nCDC42,diffCDC42,'o','MarkerEdgeColor',cmap(1,:),'LineWidth',2,'MarkerSize',8);
plot((nCDC42+1):(nCDC42+nRAC1),diffRAC1,'o','MarkerEdgeColor',cmap(2,:),'LineWidth',2,'MarkerSize',8);
plot((nCDC42+nRAC1+1):(nCDC42+nRAC1+nRHOA),diffRHOA,'o','MarkerEdgeColor',cmap(3,:),'LineWidth',2,'MarkerSize',8);
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
export_fig(outFname);
hold off;
end