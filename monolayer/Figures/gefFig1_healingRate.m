function [] = gefFig1_healingRate()
addpath(genpath('/home2/azaritsky/code/applications/monolayer/Figures/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;
% '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/POC_Nov21/YunYu1/';
% '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/allData/Angeles_20150402_14hrs_5min_AA01_2/';
dataDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/allData/Angeles_20150402_14hrs_5min_AA01_7/';
figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Fig1/EdgeEvolution/';

roisDname = [dataDname 'ROI/roi/'];

timePerFrame = 5;
pixelSize = 1.25;

nFrames = 99; 
nPropFrames = 25; % 0-120 minutes
% nEdgeEvolution = 100;

load([roisDname '001_roi.mat']); % ROI
area0 = sum(ROI(:));
edgeLength = size(ROI,1);

healingRate = nan(1,nFrames);
edgeEvolution = nan(1,nFrames);

healingRate(1) = 0;
edgeEvolution(1) = 0;

%  params.toMuPerHour * nDiffPixels / size(ROI0,1);
% params.toMuPerHour = params.pixelSize * 60/(params.timePerFrame*params.frameJump);

for it = 2 : nFrames    
    load(sprintf('%s%03d_roi.mat', roisDname,it-1)); % ROI 
    prevArea = sum(ROI(:));
    load(sprintf('%s%03d_roi.mat', roisDname,it)); % ROI 
    curArea = sum(ROI(:));
    healingRate(it) = ((curArea - prevArea) / edgeLength) * pixelSize * (60/timePerFrame);
    edgeEvolution(it) = ((curArea - area0) / edgeLength) * pixelSize;
end


% [left bottom width height]
FPosition = [0 0 350 350];
APosition = [0.2 0.2 0.75 0.75]; 

fontsize = 10;

%% healing rate
h = figure;
hold on;
plot(1:nPropFrames,healingRate(1:nPropFrames),'o','MarkerEdgeColor','k','MarkerSize',6,'LineWidth',2);
xlabel('Time (minutes)','FontSize',fontsize);
ylabel('Healing rate (\mum hour^{-1})','FontSize',fontsize);
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[1,25]);
set(haxes,'XTick',[1,13,25]);
set(haxes,'XTickLabel',[0,12,24]*timePerFrame);
% set(haxes,'YLim',[-0.07,0.83]);
% set(haxes,'YTick',0:0.4:0.8);
% set(haxes,'YTickLabel',0:0.4:0.8);
set(h,'Color','w');
set(h,'Position',FPosition,'PaperPositionMode','auto');
axisHandle= findobj(h,'type','axes');
set(axisHandle,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
set(get(axisHandle,'XLabel'),'FontSize',fontsize); set(get(axisHandle,'YLabel'),'FontSize',fontsize);
hold off;
export_fig([figDname 'healingRate.eps']);

%% edge evolution
h = figure;
hold on;
plot(1:nPropFrames,edgeEvolution(1:nPropFrames),'o','MarkerEdgeColor','k','MarkerSize',6,'LineWidth',2);
xlabel('Time (minutes)','FontSize',fontsize);
ylabel('Edge evolution (\mum)','FontSize',fontsize);
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[1,25]);
set(haxes,'XTick',[1,13,25]);
set(haxes,'XTickLabel',[0,13,25]*timePerFrame);
% set(haxes,'YLim',[-0.07,0.83]);
% set(haxes,'YTick',0:0.4:0.8);
% set(haxes,'YTickLabel',0:0.4:0.8);
set(h,'Color','w');
set(h,'Position',FPosition,'PaperPositionMode','auto');
axisHandle= findobj(h,'type','axes');
set(axisHandle,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
set(get(axisHandle,'XLabel'),'FontSize',fontsize); set(get(axisHandle,'YLabel'),'FontSize',fontsize);
hold off;
export_fig([figDname 'edgeEvolution.eps']);

%% edge evolution all
% FPosition = [0 0 350 300];
% APosition = [0.2 0.2 0.75 0.75]; 

h = figure;
hold on;
plot(1:nFrames,edgeEvolution,'-k','LineWidth',2);%
xlabel('Time (minutes)','FontSize',fontsize);
ylabel('Edge evolution (\mum)','FontSize',fontsize);
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[1,99]);
set(haxes,'XTick',[1,21,41,61,81]);
set(haxes,'XTickLabel',[0,20,40,60,80]*timePerFrame);
% set(haxes,'YLim',[-0.07,0.83]);
% set(haxes,'YTick',0:0.4:0.8);
% set(haxes,'YTickLabel',0:0.4:0.8);
set(h,'Color','w');
set(h,'Position',FPosition,'PaperPositionMode','auto');
axisHandle= findobj(h,'type','axes');
set(axisHandle,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
set(get(axisHandle,'XLabel'),'FontSize',fontsize); set(get(axisHandle,'YLabel'),'FontSize',fontsize);
hold off;
export_fig([figDname 'edgeEvolutionAll.eps']);


end