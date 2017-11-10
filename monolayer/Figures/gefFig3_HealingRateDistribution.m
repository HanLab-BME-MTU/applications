function [] = gefFig3_HealingRateDistribution()
addpath(genpath('/home2/azaritsky/code/applications/monolayer/Figures/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Fig3/';

load('/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ScreenFinal/kd0time200/healingRate/control/healingRateDistribution.mat');

% [left bottom width height]
FPosition = [0 0 300 200];
APosition = [0.2 0.2 0.75 0.75]; 

fontsize = 10;

h = figure;
hold on;
bar(centers,distribution,'k');
xlabel('Healing rate (\mum hour{-1})','FontSize',22);
ylabel('Frequency','FontSize',22);
haxes = get(h,'CurrentAxes');
set(haxes,'XTick',0:20:80);
set(haxes,'XTickLabel',0:20:80);
set(haxes,'YTick',0:0.1:0.1);
set(haxes,'YTickLabel',0:0.1:0.1);
set(haxes,'FontSize',fontsize);
set(h,'Color','w');
set(h,'Position',FPosition,'PaperPositionMode','auto');
axisHandle= findobj(h,'type','axes');
set(axisHandle,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
set(get(axisHandle,'XLabel'),'FontSize',fontsize); set(get(axisHandle,'YLabel'),'FontSize',fontsize);
export_fig([figDname 'healingRateDistribution.eps']);

fprintf('Wound healing rate (n = %d): mean = %.1f, std = %.1f\n',length(healingRates),mean(healingRates),std(healingRates));

end