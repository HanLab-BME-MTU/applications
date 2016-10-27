function [] = gefFig2_controlFeatures()
addpath(genpath('/home2/azaritsky/code/applications/monolayer/Figures/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Fig2/';

load('/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ScreenFinal/pcaParamsScreen.mat');

n = size(outControls.speed.curFeatsNorm,1);

% [left bottom width height]
FPosition = [0 0 500 200];
APosition = [0.1 0.2 0.75 0.75]; 

fontsize = 10;

h = figure; 
imagesc(outControls.speed.curFeatsNorm(randperm(n),:)'); 
hold on;
colormap('jet'); 
caxis([-1,1]);
colorbar; 
xlabel('Experiments','FontSize',fontsize);
ylabel('Features','FontSize',fontsize);
haxes = get(h,'CurrentAxes');
% set(haxes,'XLim',[1,25]);
set(haxes,'XTick',100:100:300);
set(haxes,'XTickLabel',100:100:300);
% set(haxes,'YLim',[-0.07,0.83]);
% set(haxes,'YTick',0:0.4:0.8);
set(haxes,'YTickLabel',[]);
set(h,'Color','w');
set(h,'Position',FPosition,'PaperPositionMode','auto');
axisHandle= findobj(h,'type','axes');
set(axisHandle,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
set(get(axisHandle,'XLabel'),'FontSize',fontsize); set(get(axisHandle,'YLabel'),'FontSize',fontsize);
hold off;
export_fig([figDname 'controlNormalizedFeatures.eps']);
end