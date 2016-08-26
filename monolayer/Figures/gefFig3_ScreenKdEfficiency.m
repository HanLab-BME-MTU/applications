function [] = gefFig3_ScreenKdEfficiency()
addpath(genpath('/home2/azaritsky/code/applications/monolayer/Figures/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Fig3/';

load('/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ScreenFinal/kd0time200/hairpinEfficiency/shSeqEfficiencyDistribution.mat');

% [left bottom width height]
FPosition = [0 0 300 200];
APosition = [0.2 0.2 0.75 0.75]; 

fontsize = 10;

h = figure;
hold on;
bar(xcenters,allSeqKnockdownDistribution,'k');
xlabel('GEF KD efficiency (%)','FontSize',fontsize);
ylabel('Frequency','FontSize',fontsize);
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[0,100]);
set(haxes,'XTick',0:25:100);
set(haxes,'XTickLabel',0:25:100);
set(haxes,'YLim',[0,0.17]);
set(haxes,'YTick',0:0.1:0.1);
set(haxes,'YTickLabel',0:0.1:0.1);
set(haxes,'FontSize',fontsize);
set(h,'Color','w');
set(h,'Position',FPosition,'PaperPositionMode','auto');
axisHandle= findobj(h,'type','axes');
set(axisHandle,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
set(get(axisHandle,'XLabel'),'FontSize',fontsize); set(get(axisHandle,'YLabel'),'FontSize',fontsize);
hold off;
outFname = [figDname 'shSeqEfficiencyDistribution.eps'];
export_fig(outFname);
end