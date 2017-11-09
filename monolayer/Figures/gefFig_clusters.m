function [] = gefFig_clusters()

addpath(genpath('/home2/azaritsky/code/applications/monolayer/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

inDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/allData/Angeles_20140308_16hr_5min_0001_0002_AB01_03/coordination/';
kymoFname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/allData/kymographs/coordination/Angeles_20140308_16hr_5min_0001_0002_AB01_03_coordinationKymograph.mat';
figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/FigCoordination/';

timePerFrame = 5;
pixelSize = 1.25;

nFrames = floor(200/timePerFrame);
nTime = nFrames; % same name, just for conviniency

% [left bottom width height]
FPosition = [0 0 350 350];
APosition = [0.2 0.2 0.75 0.75];

fontsize = 10;

%% norm area
coordinationOverTime = [];
for i = 1 : nFrames
    load([inDname sprintf('%03d',i) '_coordination.mat']);
    load([inDname '/../ROI/roi/' sprintf('%03d',i) '_roi.mat']);
    coordinationOverTime = [coordinationOverTime sum(ROIclusters(:))./sum(ROI(:))];
end

h = figure;
hold on;
plot(1:nFrames,coordinationOverTime,'o','MarkerEdgeColor','k','MarkerSize',6,'LineWidth',2);
xlabel('Time (minutes)','FontSize',fontsize);
ylabel('Coordintaion (pcnt)','FontSize',fontsize);
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[1,nFrames]);
set(haxes,'XTick',[1,20,40]);
set(haxes,'XTickLabel',[1,20,40]*timePerFrame);
% set(haxes,'YLim',[-0.07,0.83]);
% set(haxes,'YTick',0:0.4:0.8);
% set(haxes,'YTickLabel',0:0.4:0.8);
set(h,'Color','w');
set(h,'Position',FPosition,'PaperPositionMode','auto');
axisHandle = findobj(h,'type','axes');
set(axisHandle,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
set(get(axisHandle,'XLabel'),'FontSize',fontsize); set(get(axisHandle,'YLabel'),'FontSize',fontsize);
hold off;
export_fig([figDname 'coordinationOverTime.eps']);

[rCoord,pCoord] = corr((1:nFrames)',coordinationOverTime');

%% kymograph
load(kymoFname);
nDist = 12; % 12 x 15 = 180um
kymograph = coordinationKymograph(1:nDist,1:nTime);
params.fname = [figDname 'coordinationKymograph.eps'];
params.caxis = [0 1];
params.timePerFrame = timePerFrame;
params.patchSize = ceil(15.0/pixelSize); % 15 um in pixels
plotKymograph(kymograph,params);

extraColor = 100;

%% Plot coordination
for i = 1 : nTime
    imgFname = [inDname '/../images/' sprintf('%03d',i) '.tif'];
    
    I = imread(imgFname);
    [sizeY,sizeX] = size(I);
    %     load([inDname '/../ROI/roi/' sprintf('%03d',i) '_roi.mat']); % ROI
    load([inDname sprintf('%03d',i) '_coordination.mat']);% ROIclusters
    
    % motionClustersMap = imresize(ROIclusters,[sizeY,sizeX]);
    % motionClustersContours = bwboundaries(motionClustersMap);
    
    I8 = imadjust(I);
    Imotion = uint8(zeros(sizeY,sizeX,3));%16
    I8red = I8;
    %     I8red(motionClustersMap) = min(I8(motionClustersMap) + 100,255);
    I8red(ROIclusters) = min(I8(ROIclusters) + extraColor, 255);% 65535
    Imotion(:,:,1) = I8red;
    Imotion(:,:,2) = I8;
    Imotion(:,:,3) = I8;
    hMotion = figure;
    imagesc(Imotion);
    hold on;
    haxes = get(hMotion,'CurrentAxes');
    set(haxes,'XTick',[]);
    set(haxes,'YTick',[]);
    set(haxes,'XTickLabel',[]);
    set(haxes,'YTickLabel',[]);
    set(hMotion,'Color','none');
    axis equal; % no image stretching
    axis off;
    axis tight;
    hold off;
    export_fig([figDname filesep 'visPerTime' filesep sprintf('%.3d',i) '_clusters.eps']);
    % savefig(hMotion,[dirs.visMotionStressClustersDir sprintf('%.3d',t) '_velocityClusters.fig']);
    
    close all;
end

end