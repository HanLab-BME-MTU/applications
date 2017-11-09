cd('P:\Basic_Sciences\CMB\GoldmanLab\Takeshi\N-SIM\041015\MEFLB1-LAmAb414-006_Reconstructed');
MD = MovieData.load('MEFLB1-LAmAb414-006_Reconstructed.mat');

import lamins.classes.*;

p.x = [571.310319085489 650.321267541018];
p.y = [381.488029332372 460.498977787901];

figpos = [100 1214 117.333333333333 117.333333333333];

LD = LaminsData(MD);
LD.params.steerable.sigma = 2;
images = LD.getImages;
LI = images(1,1,6);

%% Setup Folder

mkdir('fig1');

%% Configure Positions

pos(1).x = [581.8423  591.9164];
pos(1).y = [408.8796  418.9537];
pos(1).nrect = [0.3604 0.2666 0.1240 0.1240];

pos(2).x = [599.4928  609.6502];
pos(2).y = [429.3608  439.5182];
pos(2).nrect = [0.7129 0.7910 0.1240 0.1240];

pos(3).x = [627.5217  637.3877];
pos(3).y = [388.0377  397.9037];
pos(3).nrect = [0.1357 0.5273 0.1240 0.1240];

%% Make figures

hfig = figure; imshow(LI,[]);
hfig.Position = figpos + [figpos(3) 0 0 0];
xlim(p.x); ylim(p.y);
annotation('rectangle',[0.3604 0.2666 0.1240 0.1240],'Color',[1 0 1],'LineWidth',1)
annotation('rectangle',[0.7129 0.7910 0.1240 0.1240],'Color',[1 0 1],'LineWidth',1)
annotation('rectangle',[0.1357 0.5273 0.1240 0.1240],'Color',[1 0 1],'LineWidth',1);
saveas(hfig,['fig1' filesep '01.png']);

hfig = figure; imshow(LI.steerable.nms,[]);
hfig.Position = figpos + [figpos(3) 0 0 0]*2;
xlim(p.x); ylim(p.y);
annotation('rectangle',[0.3604 0.2666 0.1240 0.1240],'Color',[1 0 1],'LineWidth',1)
annotation('rectangle',[0.7129 0.7910 0.1240 0.1240],'Color',[1 0 1],'LineWidth',1)

annotation('rectangle',[0.1357 0.5273 0.1240 0.1240],'Color',[1 0 1],'LineWidth',1);
saveas(hfig,['fig1' filesep '02.png']);

hfig = figure; imshow(LI.steerable.theta,[]);
hfig.Position = figpos + [figpos(3) 0 0 0]*3;
xlim(p.x); ylim(p.y);
colormap(gca,hsv);
annotation('rectangle',[0.3604 0.2666 0.1240 0.1240],'Color',[0 0 0],'LineWidth',1)
annotation('rectangle',[0.7129 0.7910 0.1240 0.1240],'Color',[0 0 0],'LineWidth',1)
annotation('rectangle',[0.1357 0.5273 0.1240 0.1240],'Color',[0 0 0],'LineWidth',1);
saveas(hfig,['fig1' filesep '03.png']);

hfig = figure; imshow(LI.steerable.res,[]);
hfig.Position = figpos + [figpos(3) 0 0 0]*4;
xlim(p.x); ylim(p.y);
annotation('rectangle',[0.3604 0.2666 0.1240 0.1240],'Color',[1 0 1],'LineWidth',1)
annotation('rectangle',[0.7129 0.7910 0.1240 0.1240],'Color',[1 0 1],'LineWidth',1)
annotation('rectangle',[0.1357 0.5273 0.1240 0.1240],'Color',[1 0 1],'LineWidth',1);
saveas(hfig,['fig1' filesep '04.png']);

%% Make Zoom Figures

[X,Y] = meshgrid(1:size(LI.image,2),1:size(LI.image,1));

hfig = figure; imshow(LI.steerable.nms,[]);
hfig.Position = figpos + [figpos(3) 0 0 0]*5;
xlim(pos(1).x);
ylim(pos(1).y);
hold on; quiver( ...
    X-0.5*cos(LI.steerable.theta+pi/2), ...
    Y-0.5*sin(LI.steerable.theta+pi/2), ...
    cos(LI.steerable.theta+pi/2), ...
    sin(LI.steerable.theta+pi/2), ...
    0,'y','ShowArrowhead','off');
saveas(hfig,['fig1' filesep '05.png']);

hfig = figure; imshow(LI.steerable.nms,[]);
hfig.Position = figpos + [figpos(3) 0 0 0]*6;
xlim(pos(2).x);
ylim(pos(2).y);
hold on; quiver( ...
    X-0.5*cos(LI.steerable.theta+pi/2), ...
    Y-0.5*sin(LI.steerable.theta+pi/2), ...
    cos(LI.steerable.theta+pi/2), ...
    sin(LI.steerable.theta+pi/2), ...
    0,'y','ShowArrowhead','off');
saveas(hfig,['fig1' filesep '06.png']);

hfig = figure; imshow(LI.steerable.nms,[]);
hfig.Position = figpos + [figpos(3) 0 0 0]*7;
xlim(pos(3).x);
ylim(pos(3).y);
hold on; quiver( ...
    X-0.5*cos(LI.steerable.theta+pi/2), ...
    Y-0.5*sin(LI.steerable.theta+pi/2), ...
    cos(LI.steerable.theta+pi/2), ...
    sin(LI.steerable.theta+pi/2), ...
    0,'y','ShowArrowhead','off');
saveas(hfig,['fig1' filesep '07.png']);