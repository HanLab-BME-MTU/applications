%% ---- Settings --- %%
%Figures made for paper under matlab r2011b with ubuntu 12.04 LTS
%1/2015 - Updated to matlab r2014a

saveFigs = true;



%% ------- Common Parameters ------ %%



% ---- Fig Output ---- %

%Parent directory for all panels
%figParentDir = '/home/he19/.gvfs/idac on nucleus/Hunter/orhcestra_home_june30_and_backup_merged/home/Papers/windowing methods paper/Figures/Panels';%from Ubuntu desktop, but has export problems.
figParentDir = 'W:\Hunter\orhcestra_home_june30_and_backup_merged\home\Papers\windowing methods paper\Figures\Panels';%From windows PC
fdName = 'Figure';

%Printing options
pOptEPS = {'-depsc','-r600'};
pOptTIFF = {'-dtiff','-r600'};
%options for the export_fig function
expFigOps = {'-eps','-tiff','-r300'};



%Display options
satPct = 1;

% ---- Data ---- %%

iErkChan = 1;
iWaveChan = 2;
iActChan = 3;
iDAPIChan = 4;

%Parameters for plots/legends/axes labels
plotPars = {'LineWidth',3,'MarkerSize',10};%Parameters for line plots
axLabPars = {'FontSize',16,'FontName','Arial'};%Parameters for axes labels
axPars = {'FontSize',14,'FontName','Arial'};%Parameters for axes 
cBarPars = {'FontSize',14,'FontName','Arial'};
scBarPars = {'FontName','Arial','FontSize',14};
legPars = {'FontSize',16};%Param for legends

cBarFigPars  = {'PaperPositionMode','auto'};%colorbar figure props. Prevents font resizing on EPS export


fpOff = .001;


%Function for getting roi commands from current axis limits
roiStrFunX = @(x)(clipboard('copy',['roi' num2str(x) 'X = [' num2str(xlim) '];']));
roiStrFunY = @(x)(clipboard('copy',['roi' num2str(x) 'Y = [' num2str(ylim) '];']));

%% %%%%%%%% ================ FIGURE 1 ======================== %%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Protrusion calculation illustration stuff

%Use example Arp2/3 cell per G's request
%arpList = '/home/he19/files/LCCB/gtpases/Hunter/methods_paper_data/Kwonmoo_Arp3/movieListHaloOnly.mat';%From ubuntu desktop
arpList = 'P:\gtpases\Hunter\methods_paper_data\Kwonmoo_Arp3\movieListHaloOnly.mat';%From windows desktop
if ~exist('MLarp','var')
    MLarp = MovieList.load(arpList,0);
end

scBarSz = 5e3; %Scale bar size in nm

%% ------Arp  Example image and Edge Motion

panelName = 'Arp23 example image and edge motion';
panelFile = [figParentDir filesep panelName];
panelFig = figure;
axis
panelAxes = get(panelFig,'CurrentAxes');

iArpEx = 1;
iFrame = 1;
pixSize = MLarp.movies_{iArpEx}.pixelSize_;

%MLarp.movies_{iArpEx}.channels_(1).draw(iFrame,'hAxes',panelAxes);
%Transpose so it fits in the figure better
im = MLarp.movies_{iArpEx}.channels_(1).loadImage(iFrame);
imshow(im',[]);
hold on
saturateImageColormap(panelAxes,satPct*2);
plotScaleBar(scBarSz / pixSize,'Handle',panelAxes,'Location','SouthWest');
set(panelAxes,'XDir','reverse');


prot = MLarp.movies_{iArpEx}.processes_{MLarp.movies_{iArpEx}.getProcessIndex('ProtrusionProcess',1,0)}.loadChannelOutput;
nFrames = numel(prot.smoothedEdge);
frameCols = jet(nFrames);
for j =  1:nFrames    
    plot(prot.smoothedEdge{j}(:,2),prot.smoothedEdge{j}(:,1),'Color',frameCols(j,:))    
end

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end

%% --- Colorbar for above figure

panelName = 'Arp23 example image and edge motion colorbar';
panelFile = [figParentDir filesep panelName];
panelFig = figure(cBarFigPars{:});

cBarAx = colorbar(cBarPars{:});
axis off
caxis([0 MLarp.movies_{iArpEx}.timeInterval_ * MLarp.movies_{iArpEx}.nFrames_]) 
set(get(cBarAx,'YLabel'),'String','Time, Seconds',cBarPars{:})

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end

%% ----- Edge ROI w/ vectors example
%This part still uses old example movie because it doesn't matter and who
%gives a shit

%exampMovProt = '/home/he19/files/LCCB/gtpases/Hunter/methods_paper_data/Kwonmoo_Arp3/Arp3_GFP_w_Shutter_stack2_timeCropped';%Ubuntu desktop
exampMovProt = 'P:\gtpases\Hunter\methods_paper_data\Kwonmoo_Arp3\Arp3_GFP_w_Shutter_stack2_timeCropped';%Windows desktop

if ~exist('MDprot','var')
    MDprot = MovieData.load([exampMovProt filesep 'movieData.mat'],0);
end
iFrame = 1;
tInt = MDprot.timeInterval_;
pSize = MDprot.pixelSize_;
nFrames = MDprot.nFrames_;
protVecs = MDprot.processes_{MDprot.getProcessIndex('ProtrusionProcess',1,0)}.loadChannelOutput;


scBarSz = 1e3;%scale bar size for closeup

roiX =  [401.5674  444.4368];
roiY =  [213.0369  244.0495];




panelName = 'edge displacement calc example';
panelFile = [figParentDir filesep panelName];

panelFig = fsFigure(.5);
panelAxes = gca;
%Show image
im1 = mat2gray(MDprot.channels_(1).loadImage(iFrame));
im2 = mat2gray(MDprot.channels_(1).loadImage(iFrame+1));
%im1 = im1(round(roi1Y(1)):round(roi1Y(2)),round(roi1X(1)):round(roi1X(2)));
%im2 = im2(round(roi1Y(1)):round(roi1Y(2)),round(roi1X(1)):round(roi1X(2)));

gam = 1;
f = 10;
sVal = [satPct/100 1-(satPct/100/f)];
%image(cat(3,imadjust(im1,[],stretchlim(im1,sVal),gam),imadjust(im2,[],stretchlim(im2,sVal),gam),zeros(size(im1))));
imshow(im1,[]),colormap gray
axis image
axis xy
axis off
%set(panelAxes,'YDir','reverse')
saturateImageColormap(panelAxes,satPct*3);
hold on


[ix,iy] = intersectionsHLE(protVecs.smoothedEdge{iFrame}(:,1),protVecs.smoothedEdge{iFrame}(:,2),...
                           protVecs.smoothedEdge{iFrame+1}(:,1),protVecs.smoothedEdge{iFrame+1}(:,2));


plot(protVecs.smoothedEdge{iFrame}(:,1),protVecs.smoothedEdge{iFrame}(:,2),'r',plotPars{:})
plot(protVecs.smoothedEdge{iFrame+1}(:,1),protVecs.smoothedEdge{iFrame+1}(:,2),'g',plotPars{:})
quiver(protVecs.smoothedEdge{iFrame}(:,1),protVecs.smoothedEdge{iFrame}(:,2),...
        protVecs.protrusion{iFrame}(:,1),protVecs.protrusion{iFrame}(:,2),0,plotPars{:});
plot(ix,iy,'mo',plotPars{:},'MarkerSize',15);
plotScaleBar(scBarSz / MDprot.pixelSize_,panelAxes);
xlim(roiX);
ylim(roiY);


lh = legend('cell edge \itt','cell edge \itt+1','edge mapping vectors','edge intersections');
set(lh,legPars{:},'FontName','Arial')

set(panelAxes,axPars{:})

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end

%% ----- Temporal Correspondence Illustration ---- %%
%Shows how this lets us track an edge section over time

showFrames = 10:25;
nShow = numel(showFrames);

panelName = 'edge temporal correspondence example';
panelFile = [figParentDir filesep panelName];
panelFig = fsFigure(.5);
hold on

edgeCols = gray(nShow);
edgeCols = edgeCols(end:-1:1,:);
%edgeCols(:,1) = 1;
vecCols = gray(nShow);
vecCols = vecCols(end:-1:1,:);
vecCols(:,3) = 1;

roi1X = [363.8776  437.8305];
roi1Y = [253.9915  309.1217];

for j = 1:nShow-1
    
    plot(protVecs.smoothedEdge{showFrames(j)}(:,1),protVecs.smoothedEdge{showFrames(j)}(:,2),'color',edgeCols(j,:),plotPars{:})    
    quiver(protVecs.smoothedEdge{showFrames(j)}(:,1),protVecs.smoothedEdge{showFrames(j)}(:,2),...
           protVecs.protrusion{showFrames(j)}(:,1),protVecs.protrusion{showFrames(j)}(:,2),0,'color',vecCols(j,:),plotPars{:});
end

plot(protVecs.smoothedEdge{showFrames(end)}(:,1),protVecs.smoothedEdge{showFrames(end)}(:,2),'color',edgeCols(end,:),plotPars{:})

axis image
set(gca,'XTick',[])
set(gca,'YTick',[]);
xlim(roi1X)
ylim(roi1Y)
box on

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end

%% --- Colorbar for above figure

panelName = 'edge temporal correspondence example colorbar1';
panelFile = [figParentDir filesep panelName];
panelFig = figure(cBarFigPars{:});

colormap(edgeCols)
cBarAx = colorbar(cBarPars{:});
axis off
caxis([0 MDprot.timeInterval_ * MDprot.nFrames_]) 
set(get(cBarAx,'YLabel'),'String','Time, Seconds',cBarPars{:})

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end

panelName = 'edge temporal correspondence example colorbar2';
panelFile = [figParentDir filesep panelName];
panelFig = figure;

colormap(vecCols)
cBarAx = colorbar(cBarPars{:});
axis off
caxis([0 MDprot.timeInterval_ * MDprot.nFrames_])
set(get(cBarAx,'YLabel'),'String','Time, Seconds',cBarPars{:})



if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end


%% %%%%%%%% ================ FIGURE 2 ======================== %%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Windowing Methods Figure - construction and properties

% ---- Init ---- %

%We use the same cell from the protrusion example, as per G's suggestion.
MDwin = MLarp.movies_{iArpEx};
exampMovWin =  MLarp.movies_{iArpEx}.outputDirectory_;
winSize = 2e3;
winType = 'constant_number';
pixSize = MDwin.pixelSize_;

winDir = [exampMovWin filesep 'windows_' num2str(winSize) 'nm_' winType ];
winFiles = dir([winDir filesep '*.mat']);
windows = load([winDir filesep winFiles(iFrame).name]);
windows = windows.windows;

scBarSz = 5e3; %Scale bar size in nm


%% ---- Windowing Example Cell ----- %%

panelName = 'Windowing example';
panelFile = [figParentDir filesep panelName];
panelFig = figure;
panelAxes = gca;
%Show image. Transpose to fit in figure and match with fig 1.
imWinEx = MLarp.movies_{iArpEx}.channels_(1).loadImage(iFrame);
imshow(imWinEx',[]);
hold on
saturateImageColormap(panelAxes,satPct*2);
plotScaleBar(scBarSz / pixSize,'Handle',panelAxes,'Location','SouthWest');

%Transpose the windows to match the image.
windowsTP = cellfun(@(z)(cellfun(@(y)(cellfun(@(x)(x([2 1],:)),y,'Unif',0)),z,'Unif',0)),windows,'Unif',0);

plotWindows(windowsTP,[{'r','FaceAlpha',0,'FaceColor','none','EdgeColor','r'} plotPars{:}])


%Strips and bands to highlight
iBandHi = 2;
iStripHi = 37;
xl = xlim;yl=ylim;
winAlph = 1;


plotWindows(windowsTP,[{'b','FaceAlpha',winAlph} plotPars{:}],'bandMin',iBandHi,'bandMax',iBandHi);
plotWindows(windowsTP{iStripHi},[{'g','FaceAlpha',winAlph} plotPars{:}]);
plotWindows(windowsTP{iStripHi},[{'y','FaceAlpha',1} plotPars{:}],'bandMin',iBandHi,'bandMax',iBandHi);

roi1X = [535.0331      785.0331];
roi1Y = [420.6424      640.5651];

roi2X = [9.725717      201.5701];
roi2Y = [40.20587      208.9696];

patch(roi1X([1 2 2 1]),roi1Y([1 1 2 2 ]),'w','FaceColor','none','EdgeColor',ones(1,3)-1e-3,plotPars{:})
patch(roi2X([1 2 2 1]),roi2Y([1 1 2 2 ]),'w','FaceColor','none','EdgeColor','m',plotPars{:})
%patch(roi3X([1 2 2 1]),roi3Y([1 1 2 2 ]),'w','FaceColor','none','EdgeColor','m',plotPars{:})
set(gcf, 'InvertHardcopy', 'off');
set(gcf,'DefaultLineLineSmoothing','on');
set(gcf,'DefaultPatchLineSmoothing','on');
axis off,axis tight
axis ij
set(panelAxes,'XDir','reverse')
xlim(xl);ylim(yl);


if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end

%% --- Windowing Close-Up 1 ----- %

panelName = 'Windowing close up example 1';
panelFile = [figParentDir filesep panelName];
panelFig = figure;
hold on

imHan = imshow(imWinEx',[]);
panelAxes = get(panelFig,'CurrentAxes');

hold on
saturateImageColormap(panelAxes,satPct*2);



winHan = plotWindows(windowsTP,[{'r','FaceAlpha',0,'EdgeColor','r','FaceColor','none'},'LineWidth',9]);
axis ij
set(panelAxes,'XDir','reverse');
set(gca,'XTick',[])
set(gca,'YTick',[]);
%axis on


%For some reason this is showing up behind the windwos, so get handle and
%move to top

xlim(roi1X);
ylim(roi1Y);

scBarHan = plotScaleBar(scBarSz / pixSize,'Handle',panelAxes,'Location','NorthWest');

uistack(imHan,'top')
uistack(winHan,'top')
uistack(scBarHan,'top')



if saveFigs
    %print(panelFig,panelFile,pOptTIFF{:});
    %print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
    
    export_fig(panelFile,expFigOps{:})
end

%% --- Windowing Close-Up 2 ----- %

panelName = 'Windowing close up example 2';
panelFile = [figParentDir filesep panelName];
panelFig = figure;
hold on

imHan = imshow(imWinEx',[]);
panelAxes = get(panelFig,'CurrentAxes');

hold on
saturateImageColormap(panelAxes,satPct*2);



winHan = plotWindows(windowsTP,[{'r','FaceAlpha',0,'EdgeColor','r','FaceColor','none'} ,'LineWidth',9]);
axis ij
set(panelAxes,'XDir','reverse');
set(gca,'XTick',[])
set(gca,'YTick',[]);
%axis on


%For some reason this is showing up behind the windwos, so get handle and
%move to top

xlim(roi2X);
ylim(roi2Y);

scBarHan = plotScaleBar(scBarSz / pixSize,'Handle',panelAxes,'Location','NorthWest');

uistack(imHan,'top')
uistack(winHan,'top')
uistack(scBarHan,'top')



if saveFigs
    %print(panelFig,panelFile,pOptTIFF{:});
    %print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
    
    export_fig(panelFile,expFigOps{:})
end

%% ----------- "Master Band" illustration ------- %%
% Illustrates effect of changing startcontour

winSz = 40;


mSize = [500 500];

x = 0:mSize(2);
y = sin( (x ./ mSize(2) * 2 * pi) - 0)*mSize(1)/4;

mPolyx = [0 0:mSize(2) mSize(2)];
yLoc = round(mSize(1)*2/3);
mPolyy = [0 yLoc + y 0];

exMask = poly2mask(mPolyx,mPolyy,mSize(1),mSize(2));


iContShow = [1 3 5];
nContShow = numel(iContShow);

for j = 1:nContShow

    panelName = ['Windowing master band example contour ' num2str(iContShow(j))];
    panelFile = [figParentDir filesep panelName];
    panelFig = fsFigure(.5);    
    
    imshow(exMask,[])
    hold on

    win = getMaskWindows(exMask,winSz,winSz,'StartContour',iContShow(j));
    isoCont = contourc(double(bwdist(~exMask)),0:winSz:1e3);
    isoCont = separateContours(isoCont);
    plotWindows(win,[{'r','FaceColor','none','EdgeColor','r'} plotPars{:}]);
    plot(isoCont{iContShow(j)}(1,:),isoCont{iContShow(j)}(2,:),'g',plotPars{:});

    if saveFigs
        print(panelFig,panelFile,pOptTIFF{:});
        print(panelFig,panelFile,pOptEPS{:});
        hgsave(panelFig,panelFile);
    end
    
end

%% %%%%%%%% ================ FIGURE 3 ======================== %%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Windowing propagation figure - method and activity matrix illustration

%Use same cell as 1 & 2 used to use pre-G's suggested change. Doesn't
%matter for illustration
exampMovWin = '/home/he19/files/LCCB/gtpases/Hunter/methods_paper_data/Kwonmoo_Arp3/Arp3_GFP_w_Shutter_stack2_timeCropped';
mdWin = MovieData.load([exampMovWin filesep 'movieData.mat']);

winDir = [exampMovWin filesep 'windows_' num2str(winSize) 'nm_' winType ];
winFiles = dir([winDir filesep '*.mat']);
windows = load([winDir filesep winFiles(iFrame).name]);
windows = windows.windows;

scBarSz = 5e3; %Scale bar size in nm


%% ----- Windowing Propagation Illusration ----- %%
%Closeup showing window in two frames and vectors

panelName = 'Window propagation illustration';
panelFile = [figParentDir filesep panelName];
panelFig = figure('Position',[628   326   891   669]);
hold on

winSize = 500;
winType = 'protrusion_based';

iFrame = 3;

winDir = [exampMovWin filesep 'windows_' num2str(winSize) 'nm_' winType ];
winFiles = dir([winDir filesep '*.mat']);
windows = load([winDir filesep winFiles(iFrame).name]);
windows = windows.windows;
windowst2 = load([winDir filesep winFiles(iFrame+1).name]);
windowst2 = windowst2.windows;

%iWin = [205 1];
iWin = [177 1];

plotWindows(windows{iWin(1)}{iWin(2)},[{'k','EdgeColor','none','FaceColor',.8*ones(1,3)} plotPars{:}]);    
plotWindows(windows,[{'k','FaceAlpha',0,'FaceColor','none','EdgeColor',.6*ones(1,3)} plotPars{:}])
quiver(protVecs.smoothedEdge{iFrame}(:,1),protVecs.smoothedEdge{iFrame}(:,2),...
        protVecs.protrusion{iFrame}(:,1),protVecs.protrusion{iFrame}(:,2),0,plotPars{:});
plotWindows(windowst2,[{'k','FaceAlpha',0,'FaceColor','none','EdgeColor','r'} plotPars{:}])
plotWindows(windowst2{iWin(1)}{iWin(2)},[{'r','EdgeColor','none','FaceColor','r','FaceAlpha',.3} plotPars{:}]);    

axis equal

roi1X = [544.5202      564.7019];
roi1Y = [296.5981      312.5335];

xlim(roi1X);
ylim(roi1Y);

set(gca,'Units', 'normalized', 'Position', [0+fpOff 0+fpOff 1-fpOff*2 1-fpOff*2])
set(gca,'XTick',[]);
set(gca,'YTick',[]);

set(gcf, 'InvertHardcopy', 'off');
set(gcf,'DefaultLineLineSmoothing','on');
set(gcf,'DefaultPatchLineSmoothing','on');
box on


if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end

%% ----- Windowing Propagation Example ----- %%
%Shows highlighted window and edge motion

panelName = 'Window propagation example';
panelFile = [figParentDir filesep panelName];
panelFig = figure('Position',[628   326   891   669]);
hold on

iWin = [185 1];

showFrames = 1:1:40;%Focuses on retraction event so windows can be seen easily
nShow = numel(showFrames);
nBuf = 5;%Make colormap initially larger than needed to avoid having frames be completely white (and invisible)

frameCols = gray(nShow+nBuf);
frameCols = frameCols(end-nBuf:-1:1,:);
winCols = gray(nShow+nBuf);
winCols = winCols(end-nBuf:-1:1,:);
winCols(:,1) = 1;
winBordCols = gray(nShow+2*nBuf);
winBordCols= winBordCols(end-2*nBuf:-1:1,:);
winBordCols(:,1) = 1;

vecFade = 2.5;
nVecCol = round(nShow*vecFade);
vecCols = gray(nVecCol);
vecCols = vecCols(end-nBuf:-1:end-nShow-nBuf+1,:);
vecCols(:,3) = 1;



%Do the edges first so they are behind the windows
for j = 1:nShow
    plot(protVecs.smoothedEdge{showFrames(j)}(:,1),protVecs.smoothedEdge{showFrames(j)}(:,2),'Color',frameCols(j,:),plotPars{:})
%Get's too confusing with prot vecs shown....
    quiver(protVecs.smoothedEdge{showFrames(j)}(:,1),protVecs.smoothedEdge{showFrames(j)}(:,2),...
           protVecs.protrusion{showFrames(j)}(:,1),protVecs.protrusion{showFrames(j)}(:,2),0,'color',vecCols(j,:),plotPars{:});
% end
% for j = 1:nShow
    windows = load([winDir filesep winFiles(showFrames(j)).name]);
    windows = windows.windows;
    if ~isempty(windows{iWin(1)}) && ~isempty(windows{iWin(1)}{iWin(2)})
        plotWindows(windows{iWin(1)}{iWin(2)},[{'r','FaceAlpha',1,'EdgeColor',winBordCols(j,:),'FaceColor',winCols(j,:)} plotPars{:},'LineWidth',4]);
        %plotWindows(windows{iWin(1)}{iWin(2)},[{'r','FaceAlpha',1,'EdgeColor','none','FaceColor',winCols(j,:)} plotPars{:},'LineWidth',4]);        
    end
end


roi2X = [466.9339      538.8578];
roi2Y = [241.8571      295.2817];

xlim(roi2X);
ylim(roi2Y);

set(gca,'Units', 'normalized', 'Position', [0+fpOff 0+fpOff 1-fpOff*2 1-fpOff*2])
set(gcf, 'InvertHardcopy', 'off');
set(gcf,'DefaultLineLineSmoothing','on');
set(gcf,'DefaultPatchLineSmoothing','on');
set(gca,'XTick',[]);
set(gca,'YTick',[]);
box on

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end

%% --- Colorbars for above figure


panelName = [panelName ' colorbar1'];
panelFile = [figParentDir filesep panelName];
panelFig = figure;

colormap(frameCols);
cBarAx = colorbar(cBarPars{:});
axis off
caxis([0 MDwin.timeInterval_ * MDwin.nFrames_]) 
set(get(cBarAx,'YLabel'),'String','Time, Seconds',cBarPars{:})

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end

panelName = [panelName ' colorbar2'];
panelFile = [figParentDir filesep panelName];
panelFig = figure;

colormap(winCols);
cBarAx = colorbar(cBarPars{:});
axis off
caxis([0 MDwin.timeInterval_ * MDwin.nFrames_]) 
set(get(cBarAx,'YLabel'),'String','Time, Seconds',cBarPars{:})

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end

panelName = [panelName ' colorbar2'];
panelFile = [figParentDir filesep panelName];
panelFig = figure;

colormap(vecCols);
cBarAx = colorbar(cBarPars{:});
axis off
caxis([0 MDwin.timeInterval_ * MDwin.nFrames_]) 
set(get(cBarAx,'YLabel'),'String','Time, Seconds',cBarPars{:})

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end



%% ------ Activity Map Illustration ------ %%
%Shows assembly of samples into matrix
%(Mostly made in illustrator, but some from here)

winSize = 500;
winType = 'constant_number';

winDir = [exampMovWin filesep 'windows_' num2str(winSize) 'nm_' winType ];
winFiles = dir([winDir filesep '*.mat']);
windows = load([winDir filesep winFiles(iFrame).name]);
windows = windows.windows;

winSampDir = [exampMovWin filesep 'window_samples_' num2str(winSize) 'nm_' winType ];
winSampFiles = dir([winSampDir filesep '*.mat']);
winSamps = load([winSampDir filesep winSampFiles(1).name]);
sampFields = fieldnames(winSamps);
winSamps = winSamps.(sampFields{1});

tData = 0:MDwin.timeInterval_:MDwin.timeInterval_*MDwin.nFrames_;

nBands = 4;%Number of bands to show activity maps from
nSnap = 3;%Number of frames to show windows/images from

iSnap = round(linspace(1,MDwin.nFrames_,nSnap));


%% ---- Activity Maps for Illustration ---- %%

cStp = 6;
nColTot = nBands*cStp;%Use more colors for bands than for maps so we can see differences in both
bandMasterMap = hsv(nColTot*5);%Use beginning of HSV map so you can see diff from jet used on prot map
bandMasterMap = bandMasterMap(1:nColTot,:);
bandCols = bandMasterMap(1:cStp:nBands*cStp,:);
nCols = 256;
fadeTo = .33;
bandMap = zeros(nCols,3,nBands);

for j = 1:nBands
        
    panelName = ['activity map for assembly illustration band ' num2str(j)];
    panelFile = [figParentDir filesep panelName];
    panelFig = figure('Position',[397         359        1278         502]);
    set(gcf,'PaperPositionMode', 'auto');
    hold on
    
    bandMap(:,:,j) = cat(2,linspace(bandCols(j,1),fadeTo,nCols)',linspace(bandCols(j,2),fadeTo,nCols)',linspace(bandCols(j,3),fadeTo,nCols)');            
    
    smMap = smoothActivityMap(squeeze(winSamps.avg(:,j,:)));
    imHan = imagesc(tData,1:size(winSamps.avg,1),smMap);   
    colormap(bandMap(:,:,j))
    saturateImageColormap(imHan,satPct*2);    
    if j == 1
        %Only show on first band
        xlabel('Time, Seconds',axLabPars{:})
        ylabel({'Along Cell Edge','(Slice #)'},axLabPars{:});
    else
        set(imHan,'AlphaData',~isnan(smMap))
        set(gca,'color','none')
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        set(gca,'Units', 'normalized', 'Position', [0+fpOff 0+fpOff 1-fpOff*2 1-fpOff*2])
    end
    
    xlim([min(tData) max(tData)])
    ylim([1 size(winSamps.avg,1)])    
    box on
    
    
    if j == 1
        currCa = caxis;
    else
        caxis(currCa);
    end
    
    set(gca,axPars{:})
%     set(gca,'XColor',bandCols(j,:))
%     set(gca,'YColor',bandCols(j,:))
%     
    if saveFigs
        print(panelFig,panelFile,pOptTIFF{:});
        print(panelFig,panelFile,pOptEPS{:});
        hgsave(panelFig,panelFile);
    end
    
    % ---- Colorbar for each map
    panelName = [panelName ' colorbar'];
    panelFile = [figParentDir filesep panelName];
    panelFig = figure;

    colormap(bandMap(:,:,j))
    cBarAx = colorbar(cBarPars{:});
    axis off
    caxis(currCa) 
    set(get(cBarAx,'YLabel'),'String',['Activity, A.U.'],cBarPars{:})

    if saveFigs
        print(panelFig,panelFile,pOptTIFF{:});
        print(panelFig,panelFile,pOptEPS{:});
        hgsave(panelFig,panelFile);
    end            

end



%% ----- Window Images for Illustration ----- %%
%Use bigger windows so we can see them,
%switch back to the arp example so the images are familiar
exampMovWin = MLarp.movies_{iArpEx}.outputDirectory_;



winSize = 2e3;
winType = 'constant_number';

winDir = [exampMovWin filesep 'windows_' num2str(winSize) 'nm_' winType ];
winFiles = dir([winDir filesep '*.mat']);



bandCols = bandMasterMap;
for j = 1:nSnap
    
    
    windows = load([winDir filesep winFiles(iSnap(j)).name]);
    windows = windows.windows;
    
    panelName = ['windows for activity map assembly illustration snapshot ' num2str(j)];
    panelFile = [figParentDir filesep panelName];
    panelFig = fsFigure(.5);    
    imHan = MDwin.channels_(1).draw(iSnap(j));
    saturateImageColormap(imHan,satPct*5)    
    hold on
    
    
    plotWindows(windows,[{'k','EdgeColor',ones(1,3) * .6,'FaceColor','none'} plotPars{:}])
    
    for k = 1:max(cellfun(@numel,windows))

        plotWindows(windows,[{'k','EdgeColor','none','FaceColor',bandCols(k,:),'FaceAlpha',.6} plotPars{:}],'bandMax',k,'bandMin',k)

    end
    
    if saveFigs
        print(panelFig,panelFile,pOptTIFF{:});
        print(panelFig,panelFile,pOptEPS{:});
        hgsave(panelFig,panelFile);
    end
           
    

end

%% --- Colorbar for above figure - band color coding (dist from edge)

panelName = [panelName ' colorbar'];
panelFile = [figParentDir filesep panelName];
panelFig = figure;

colormap(bandCols);
cBarAx = colorbar(cBarPars{:});
axis off
caxis([1 size(bandCols,1)]);
set(get(cBarAx,'YLabel'),'String','Distance from Cell Edge, Band #',cBarPars{:})

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end






%% ----- Protrusion Map for Matrix creation Illustration

protMapColMap = jet;
panelName = 'Protrusion map example';
panelFile = [figParentDir filesep panelName];
panelFig = figure('Position',[397         359        1278         502]);

tInt = MDprot.timeInterval_;
nFrames = MDprot.nFrames_;
uSamp = 5;%Upsampling for prot map smoothing

protSamples = load([exampMovProt filesep 'protrusion_samples_' num2str(winSize) 'nm_' winType filesep 'protrusion_samples.mat']);
protSamples = protSamples.protSamples;
tmp = protSamples.avgNormal(:,1:end-1);
tmp = tmp .* MDprot.pixelSize_ / MDprot.timeInterval_;%Convert to physical units
protErrs = detectOutliers(tmp,4);
tmp(protErrs) = NaN;
imagesc(0:tInt/uSamp:tInt*nFrames,1:size(tmp,1),smoothActivityMap(tmp,'FillNaN',1,'UpSample',uSamp));
colormap(protMapColMap);
saturateImageColormap(panelFig,satPct);
%Make sure the color axis is symmetric. Just use the max vel after saturation to keep it
%simple
caxis(max(caxis) * [-1 1])
currCa = caxis;
xlabel('Time, Seconds',axLabPars{:})
ylabel('Along Edge, Strip #',axLabPars{:})
set(gca,axPars{:});
set(gca,'YDir','normal');
set(gcf,'PaperPositionMode', 'auto');


if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end

%% --- Colorbar for above figure


panelName = [panelName ' colorbar'];
panelFile = [figParentDir filesep panelName];
panelFig = figure;

colormap(protMapColMap);
cBarAx = colorbar(cBarPars{:});
axis off
caxis(currCa);
set(get(cBarAx,'YLabel'),'String','Edge Velocity, nm/s',cBarPars{:})

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end


%% %%%%%%%% ================ FIGURE 4 ======================== %%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Whole-cell windowing and propagation comparison

%Use highly motile rac1 cell from marco
wcExampMov = '/home/he19/files/LCCB/gtpases/Hunter/methods_paper_data/rac1_cell_actuallyMoves_from_marco';
if ~exist('MDwc','var')
    MDwc = MovieData.load([wcExampMov filesep 'movieData.mat']);
end
protVecs = MDwc.processes_{MDwc.getProcessIndex('ProtrusionProcess',1,0)}.loadChannelOutput;

showFrames = 1:1:100;

winSize = 1.6e3;%Arbitratily selected for visualization. Gives 5x5 pixel windows.

%Axis lmits so Gaudenz doesn't cry about the extra white space
wcFigXlim = [50 260];
wcFigYlim = [53 351];


%% --------- Example Images Figures ------- %%
%Shows some snapshots from the movie
scBarSz = 2.5e4;%Scale bar size in nm
nSnap = 3;
iSnap = round(linspace(min(showFrames),max(showFrames),nSnap));

wcCmap = jet;
for j = 1:nSnap

    panelName = ['whole cell windowing example image ' num2str(j)];
    panelFile = [figParentDir filesep panelName];        
    panelFig = figure;
    
    im = MDwc.channels_(1).loadImage(iSnap(j));
    im = im ./ mean(im(:));
    imHan = imshow(im,[]);    
    panelAxes = get(imHan,'Parent');
    %panelFig = get(panelAxes,'Parent');
    colormap(wcCmap)        
    saturateImageColormap(imHan,satPct)    
    wcCaxis = caxis;        
    xlim(wcFigXlim)
    ylim(wcFigYlim)
    plotScaleBar(scBarSz / MDwc.pixelSize_,'Location','Northeast','Color',zeros(1,3),'Label',['t=' num2str((iSnap(j)-1) * MDwc.timeInterval_) 's'],scBarPars{:});
    im = get(imHan,'CData');    
    tmpMask = imfill(im > 0,'holes');
    tmpMask = bwareaopen(tmpMask,100);
    set(imHan,'AlphaData',tmpMask)    
    axis on
    set(panelAxes,'Color','w');
    set(panelAxes,'XTick',[]);
    set(panelAxes,'YTick',[]);

    
    if saveFigs
        print(panelFig,panelFile,pOptTIFF{:});
        print(panelFig,panelFile,pOptEPS{:});
        hgsave(panelFig,panelFile);
    end
    
end

%% --- Colorbar for above figure


panelName = [panelName ' colorbar'];
panelFile = [figParentDir filesep panelName];
panelFig = figure;

colormap(wcCmap);
cBarAx = colorbar(cBarPars{:});
axis off
caxis([1.5 2]);
set(get(cBarAx,'YLabel'),'String','Normalized FRET Ratio',cBarPars{:})

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end





%% ---------- Windowing Propagation Comparison -------- %%
%Compares two propagation methods
scBarSz = 1e4;%Scale bar size in nm

% ----- WHole-Cell Figure ---- %

panelName = 'Window propagation comparison';
panelFile = [figParentDir filesep panelName];
panelFig = fsFigure(.5);
hold on


winType1 = 'protrusion_based';
winDir1 = [wcExampMov filesep 'windows_' num2str(winSize) 'nm_' winType1 ];
winFiles1 = dir([winDir1 filesep '*.mat']);

winType2 = 'constant_number';
winDir2 = [wcExampMov filesep 'windows_' num2str(winSize) 'nm_' winType2 ];
winFiles2 = dir([winDir2 filesep '*.mat']);


nShow = numel(showFrames);

iWin = [29 1 ;
        70 1];    
    
nWinShow = size(iWin,1);    
nBuf = 10;

frameCols = gray(nShow+nBuf);
frameCols = frameCols(end-nBuf:-1:1,:);
winCols1 = gray(nShow+nBuf);
winCols1 = winCols1(end-nBuf:-1:1,:);
winCols1(:,1) = 1;
winBordCols1 = gray(nShow+2*nBuf);
winBordCols1= winBordCols1(end-2*nBuf:-1:1,:);
winBordCols1(:,1) = 1;
winCols2 = gray(nShow+nBuf);
winCols2 = winCols2(end-nBuf:-1:1,:);
winCols2(:,3) = 1;
winBordCols2 = gray(nShow+2*nBuf);
winBordCols2= winBordCols2(end-2*nBuf:-1:1,:);
winBordCols2(:,3) = 1;


%Do the edges first so they are behind the windows
for j = 1:nShow
    plot(protVecs.smoothedEdge{showFrames(j)}(:,1),protVecs.smoothedEdge{showFrames(j)}(:,2),'Color',frameCols(j,:))
end

for j = 1:nShow
    
    windows = load([winDir1 filesep winFiles1(showFrames(j)).name]);
    windows = windows.windows;
    
    for k = 1:nWinShow
        if ~isempty(windows{iWin(k,1)}) && ~isempty(windows{iWin(k,1)}{iWin(k,2)})
            plotWindows(windows{iWin(k,1)}{iWin(k,2)},[{'r','FaceAlpha',1,'EdgeColor',winBordCols1(j,:),'FaceColor',winCols1(j,:)} plotPars{:}]);        
        else
            disp(['no win 1 frame ' num2str(showFrames(j))])
        end
        
        %plotWindows(windows,[{'r','FaceAlpha',1,'EdgeColor',winCols1(j,:),'FaceColor','none'} ]);        
    end
    windows = load([winDir2 filesep winFiles2(showFrames(j)).name]);
    windows = windows.windows;
    
    for k = 1:nWinShow
        if ~isempty(windows{iWin(k,1)}) && ~isempty(windows{iWin(k,1)}{iWin(k,2)})
            plotWindows(windows{iWin(k,1)}{iWin(k,2)},[{'r','FaceAlpha',1,'EdgeColor',winBordCols2(j,:),'FaceColor',winCols2(j,:)} plotPars{:}]);
        else
            disp(['no win 2 frame ' num2str(showFrames(j))])
        end
        %plotWindows(windows,[{'r','FaceAlpha',1,'EdgeColor',winCols2(j,:),'FaceColor','none'} ]);        
    end
end


xlim(wcFigXlim)
ylim(wcFigYlim)
%ylim auto
%xlim auto

set(gca,'XTick',[]);
set(gca,'YTick',[]);
box on

roi1X = [123.679      216.6944];
roi1Y = [223.1792      353.4006];


roiCol = [0 .6 .2];
%patch(roi1X([1 2 2 1]),roi1Y([1 1 2 2 ]),roiCol,'FaceColor','none','EdgeColor',roiCol,plotPars{:})

plotScaleBar(scBarSz / MDwc.pixelSize_,'Location','Southwest','Color',zeros(1,3))

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end

%% ----- ROI of Whole-Cell Figure ---- %

scBarSz = 5e3;

panelName = 'Window propagation comparison closeup';
panelFile = [figParentDir filesep panelName];

xlim(roi1X);
ylim(roi1Y);

plotScaleBar(scBarSz / MDwc.pixelSize_,'Location','Southeast','Color',zeros(1,3));

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end

%% --- Colorbars for above figure


panelName = [panelName ' colorbar1'];
panelFile = [figParentDir filesep panelName];
panelFig = figure;

colormap(frameCols);
cBarAx = colorbar(cBarPars{:});
axis off
caxis([0 MDwc.timeInterval_ * MDwc.nFrames_]) 
set(get(cBarAx,'YLabel'),'String','Time, Seconds',cBarPars{:})

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end

panelName = [panelName ' colorbar2'];
panelFile = [figParentDir filesep panelName];
panelFig = figure;

colormap(winCols);
cBarAx = colorbar(cBarPars{:});
axis off
caxis([0 MDwc.timeInterval_ * MDwc.nFrames_]) 
set(get(cBarAx,'YLabel'),'String','Time, Seconds',cBarPars{:})

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end

panelName = [panelName ' colorbar2'];
panelFile = [figParentDir filesep panelName];
panelFig = figure;

colormap(vecCols);
cBarAx = colorbar(cBarPars{:});
axis off
caxis([0 MDwc.timeInterval_ * MDwc.nFrames_]) 
set(get(cBarAx,'YLabel'),'String','Time, Seconds',cBarPars{:})

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end




%% ----- Whole-Cell activity map for prop comparison - constant number

panelName = 'whole cell activity map for prop comparison protrusion based';
panelFile = [figParentDir filesep panelName];
panelFig = figure('Position',[397         359        1278         502]);
set(gcf,'PaperPositionMode', 'auto');


winType1 = 'protrusion_based_StartContour1_StartPoint_132_193';
sampDir1 = [wcExampMov filesep 'WindowingPackage_' winType1 filesep 'window_sampling'];
sampFiles1 = dir([sampDir1 filesep '*.mat']);
samps1 = load([sampDir1 filesep sampFiles1.name]);
samps1 = samps1.samples;

tData = 0:MDwin.timeInterval_:MDwin.timeInterval_*MDwin.nFrames_;
hold on

smMap = smoothActivityMap(squeeze(samps1.avg(:,1,:)));
imHan = imagesc(tData,1:size(samps1.avg,1),smMap);   
colormap jet
saturateImageColormap(imHan,satPct);    

xlabel('Time, Seconds',axLabPars{:})
ylabel('Along Cell Edge',axLabPars{:});

xlim([min(tData) max(tData)])
ylim([1 size(samps1.avg,1)])    
box on
set(gca,axPars{:})


%%

panelName = 'whole cell activity map for prop comparison constant number';
panelFile = [figParentDir filesep panelName];
panelFig = figure('Position',[397         359        1278         502]);
set(gcf,'PaperPositionMode', 'auto');

winType2 = 'constant_number_StartContour1_StartPoint_132_193';
sampDir2 = [wcExampMov filesep 'WindowingPackage_' winType2 filesep 'window_sampling'];
sampFiles2 = dir([sampDir2 filesep '*.mat']);
samps2 = load([sampDir2 filesep sampFiles2.name]);
samps2 = samps2.samples;

tData = 0:MDwin.timeInterval_:MDwin.timeInterval_*MDwin.nFrames_;
hold on

smMap = smoothActivityMap(squeeze(samps2.avg(:,1,:)));
imHan = imagesc(tData,1:size(samps2.avg,1),smMap);   
colormap jet
saturateImageColormap(imHan,satPct);    

xlabel('Time, Seconds',axLabPars{:})
ylabel('Along Cell Edge',axLabPars{:});

xlim([min(tData) max(tData)])
ylim([1 size(samps2.avg,1)])    
box on
set(gca,axPars{:})


%% %%%%%%%% ================ FIGURE 5 ======================== %%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Whole-cell pattern extraction 

wcExampMov = '/home/he19/files/LCCB/gtpases/Hunter/methods_paper_data/rac1_cell_actuallyMoves_from_marco';
%Use highly motile rac1 cell from marco for time illustration
if ~exist('MDwc','var')
    MDwc = MovieData.load([wcExampMov filesep 'movieData.mat']);
end
%
% protVecs = MDwc.processes_{MDwc.getProcessIndex('ProtrusionProcess',1,0)}.loadChannelOutput;
winSize = 1.6e3;%Arbitratily selected for visualization. Gives 5x5 pixel windows.

scBarSz = 5e3; %Scale bar size in nm



%% ------ Rac1 Cell Spatial pattern circle mapping Timepoint Comparison ----


showFrames = 1:1:100;
nSnap = 3;
iSnap = round(linspace(min(showFrames),max(showFrames),nSnap));


winType2 = 'constant_number';
sampDir2 = [wcExampMov filesep 'window_samples_' num2str(winSize) 'nm_' winType2 ];
sampFiles2 = dir([sampDir2 filesep '*.mat']);
samps2 = load([sampDir2 filesep sampFiles2.name]);
samps2 = samps2.samples;

nBandShow = size(samps2.avg,2);
sampRng = zeros(nSnap,2);
for j = 1:nSnap

    panelName = ['whole cell spatial pattern over time snapshot ' num2str(j)];
    panelFile = [figParentDir filesep panelName];
       
    currSamp = squeeze(samps2.avg(:,1:nBandShow,iSnap(j)));
    smMap = smoothActivityMap(currSamp,'FillNaN',true);        
    
    sampRng(j,:) = prctile(currSamp(:),[satPct (100-satPct)]);        
    panelFig = figure('Position',[ 984   368   876   806]);
    plotActivityMapPolar(smMap);           
    caxis(sampRng(j,:))
    axis off
    axis image
        
    if saveFigs
        print(panelFig,panelFile,pOptTIFF{:});
        print(panelFig,panelFile,pOptEPS{:});
        hgsave(panelFig,panelFile);
    end

end

%% ----- Activity matrix for circle-mapping method illustration

panelName = 'whole cell activity matrix for circle mapping illustration';
panelFile = [figParentDir filesep panelName];

panelFig = figure('Position',[397         490        1878         588]);

winType = 'constant_number';

iSnapUse = 1;%Number of snapshot to use activity matrix from

sampDir = [wcExampMov filesep 'window_samples_' num2str(winSize) 'nm_' winType ];
sampFiles = dir([sampDir filesep '*.mat']);
samps = load([sampDir filesep sampFiles.name]);
samps = samps.samples;

sampShow = samps.avg(:,:,iSnapUse(1))';
tmpMask = ~isnan(sampShow);
imHan = imagesc(sampShow);
set(imHan,'AlphaData',tmpMask)
panelAxes = get(panelFig,'CurrentAxes');
set(panelAxes,axPars{:});
axis image
set(panelAxes,'YDir','normal')
xlabel('Along cell edge, Slice #',axLabPars{:})
ylabel('Away from cell edge, Band #',axLabPars{:})


 if saveFigs
    
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});

     %export_fig(panelFile,expFigOps{:});
    hgsave(panelFig,panelFile);
end



%% ---- Erk Whole-cell Image, circle map and windowing comparison  -----

exList = '/home/he19/files/LCCB/gtpases/Hunter/methods_paper_data/Michelle_IF_Erk_Wave_converted/movieListUseHMEC.mat';
if ~exist('MLerk','var')
    MLerk = MovieList.load(exList,0);
end

scBarSz = 1e4;
scBarLoc = 'Northeast';
scBarCol = zeros(1,3);


iEx = [1 4 7 6 8 3];%Weird order to preserve filenames and linking in illustrator
isEGF = [1 1 0 0 0 1]>0;%if EGF stimulated
nEx = numel(iEx);

%subplot(3,nEx,1);
nShow = 20;
intSize = [200,80];
intSamps = nan([intSize, nEx]);

for j = 1:nEx
    
    panelName = 'Erk circle mapping';
                
    iWinProc = MLerk.movies_{iEx(j)}.getProcessIndex('WindowingProcess',1,0);
    windows = MLerk.movies_{iEx(j)}.processes_{iWinProc}.loadChannelOutput(1);
    iSampProc = MLerk.movies_{iEx(j)}.getProcessIndex('WindowSamplingProcess',1,0);
    samps(j) = MLerk.movies_{iEx(j)}.processes_{iSampProc}.loadChannelOutput(iErkChan);
    sampSize(j,:) = size(samps(j).avg);
    [Xi,Yi] = meshgrid(linspace(1,sampSize(j,2),intSize(2)),linspace(1,sampSize(j,1),intSize(1)));
    intSamps(:,:,j) = interp2(samps(j).avg,Xi,Yi);

    iMaskProc = MLerk.movies_{iEx(j)}.getProcessIndex('MaskRefinementProcess',1,0);
    currMask = MLerk.movies_{iEx(j)}.processes_{iMaskProc}.loadChannelOutput(iActChan,1);
    
    sampRng = prctile(samps(j).avg(:),[satPct/2 (100-satPct/2)]);
    satSamps = samps(j).avg;
    satSamps(satSamps>sampRng(2)) = sampRng(2);
    satSamps(satSamps<sampRng(1)) = sampRng(1);        
    
    cMap = hot;
    
    panelFile = [figParentDir filesep panelName ' image cell ' num2str(j)];
    panelFig = figure('Position',[810         262        1021         861]);
    hold on            
    %imageViewer(MLerk.movies_{iEx(j)},'ChannelIndex',iErkChan,'Saturate',satPct/100,'AxesHandle',gca);colormap(cMap)    
    imHan = MLerk.movies_{iEx(j)}.channels_(iErkChan).draw(1);
    set(imHan,'AlphaData',currMask);
    colormap(cMap)
    saturateImageColormap(imHan,satPct);
    plotScaleBar(scBarSz / MLerk.movies_{iEx(j)}.pixelSize_,'Location',scBarLoc,'Color',scBarCol);
    axis xy
    axis on
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    set(gca,'Color','w')
    
    if saveFigs
        print(panelFig,panelFile,pOptTIFF{:});
        print(panelFig,panelFile,pOptEPS{:});
        hgsave(panelFig,panelFile);
    end

    
    panelFile = [figParentDir filesep panelName ' windows cell ' num2str(j)];
    panelFig = figure('Position',[810         262        1021         861]);
    
    plotWindowsColormapped(windows,satSamps,cMap);    
    xlim([1 MLerk.movies_{iEx(j)}.imSize_(2)]);
    ylim([1 MLerk.movies_{iEx(j)}.imSize_(1)]);
    plotScaleBar(scBarSz / MLerk.movies_{iEx(j)}.pixelSize_,'Location',scBarLoc,'Color',scBarCol);
    axis on
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    set(gca,'Color','w')
    box off
    
    if saveFigs
        print(panelFig,panelFile,pOptTIFF{:});
        print(panelFig,panelFile,pOptEPS{:});
        hgsave(panelFig,panelFile);
    end
    
    panelFile = [figParentDir filesep panelName ' circle map cell ' num2str(j)];
    panelFig = figure('Position',[810         262        1021         861]);
    
    %nShowCurr = min(size(satSamps,2),nShowMax);
    %smMap = smoothActivityMap(intSamps(:,1:nShowCurr,j),'FillNaN',true);
    smMap = smoothActivityMap(intSamps(:,1:nShow,j),'FillNaN',true);
    plotActivityMapPolar(smMap);
    colormap(cMap)
    caxis(sampRng);
    axis off,axis equal 

    if saveFigs
        print(panelFig,panelFile,pOptTIFF{:});
        print(panelFig,panelFile,pOptEPS{:});
        hgsave(panelFig,panelFile);
    end    
    
    panelName = [panelName ' colorbar'];
    panelFile = [figParentDir filesep panelName];
    panelFig = figure;

    colormap(cMap);
    cBarAx = colorbar(cBarPars{:});
    axis off
    caxis(sampRng);
    set(get(cBarAx,'YLabel'),'String','p-Erk Fluorescence, a.u.',cBarPars{:})

    
    if saveFigs
        print(panelFig,panelFile,pOptTIFF{:});
        print(panelFig,panelFile,pOptEPS{:});
        hgsave(panelFig,panelFile);
    end
    
    
    
end

%% %%%%%%%% ================ FIGURE 6 ======================== %%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Erk data spatial pattern averaging/comparison/quant


%% ------ pErk and Wave Localization ------ %%

erkExample = '/home/he19/files/LCCB/gtpases/Hunter/methods_paper_data/Michelle_IF_Erk_Wave_converted/IFdata-new_HMEC 2.5 min_cell 10';
if ~exist('MDerk','var')
    MDerk = MovieData.load([erkExample filesep 'movieData.mat']);
end

scBarSz = 1e4;

%% --- Erk Raw Images Panel ---- %%


panelName = 'Erk raw images';
panelFile = [figParentDir filesep panelName];
panelFig = figure;

imDirs = MDerk.getChannelPaths;
imNames = MDerk.getImageFileNames;

erkImage = imread([imDirs{iErkChan} filesep imNames{iErkChan}{1}]);
waveImage = imread([imDirs{iWaveChan} filesep imNames{iWaveChan}{1}]);

erkLim = stretchlim(erkImage,satPct*2/100);
waveLim = stretchlim(waveImage,satPct*2/100);

imshow(cat(3,mat2gray(imadjust(erkImage,erkLim)),mat2gray(imadjust(waveImage,waveLim)),zeros(size(erkImage)) ))
plotScaleBar(scBarSz / MDerk.pixelSize_,'Location','southeast')
if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    hgsave(panelFig,panelFile);
end


%% ---- Raw Erk Ratio Panel ---- %
panelName = 'Erk raw ratio';
panelFile = [figParentDir filesep panelName];
panelFig = figure;



imshow(double(erkImage) ./ double(waveImage))
rawErkCmap = jet;
colormap(rawErkCmap);
saturateImageColormap(panelFig,satPct);
plotScaleBar(scBarSz / MDerk.pixelSize_,'Location','southeast')
rawErkCaxis = caxis;
erkRatCmap = colormap;
erkRatPos = get(panelFig,'Position');

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    hgsave(panelFig,panelFile);
end

%% --- Colorbar for above figure


panelName = [panelName ' colorbar'];
panelFile = [figParentDir filesep panelName];
panelFig = figure;

colormap(rawErkCmap);
cBarAx = colorbar(cBarPars{:});
axis off
caxis(rawErkCaxis);
set(get(cBarAx,'YLabel'),'String','p-Erk:Wave2 Ratio',cBarPars{:})

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end

%% ---- Windowed Erk/Wave RGB  Panel ----- %

panelName = 'Erk raw RGB windowed';
panelFile = [figParentDir filesep panelName];
panelFig = figure('Position',erkRatPos);
scBarCol = [0 0 0 ];

nCol = 1024;

iSampProc = MDerk.getProcessIndex('WindowSamplingProcess',1,0);
sampsErk = MDerk.processes_{iSampProc}.loadChannelOutput(iErkChan);
sampsErk = sampsErk.avg;
sampsWave = MDerk.processes_{iSampProc}.loadChannelOutput(iWaveChan);
sampsWave = sampsWave.avg;

iWinProc = MDerk.getProcessIndex('WindowingProcess',1,0);
windows = MDerk.processes_{iWinProc}.loadChannelOutput(1);


%Do the fucking saturation manually so I don't go nuts with the RGB
%conversion etc. This way it matches with above images.
erkMax = erkLim(2) * 65535+1;%Undo the stupid stretchlim conversion
erkMin = erkLim(1) * 65535+1;
sampsErk(sampsErk < erkMin | isnan(sampsErk)) = erkMin;
sampsErk(sampsErk > erkMax) = erkMax;
waveMax = waveLim(2) * 65535+1;%Undo the stupid stretchlim conversion
waveMin = waveLim(1) * 65535+1;
sampsWave(sampsWave < waveMin | isnan(sampsWave)) = waveMin;
sampsWave(sampsWave > waveMax) = waveMax;

rgbIm = cat(3,mat2gray(sampsErk),...
              mat2gray(sampsWave),zeros(size(sampsErk)));
%[sampRGB,erkWinRawMap] = rgb2ind(rgbIm,nCol,'nodither');

plotWindowsColormapped(windows,rgbIm);

panelAxes = get(panelFig,'CurrentAxes');

xlim([1 MDerk.imSize_(2)]);
ylim([1 MDerk.imSize_(1)]);
axis ij
plotScaleBar(scBarSz / MDerk.pixelSize_,'Location','Southeast','Color',scBarCol);
set(panelAxes,'XTick',[]);
set(panelAxes,'YTick',[]);
set(panelAxes,'color','w');

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end



%% ---- Windowed Erk Ratio Panel ---- %
panelName = 'Erk windowed ratio';
panelFile = [figParentDir filesep panelName];
panelFig = figure('Position',erkRatPos);

erkExample = '/home/he19/files/LCCB/gtpases/Hunter/methods_paper_data/Michelle_IF_Erk_Wave_converted/IFdata-new_HMEC 2.5 min_cell 10';
if ~exist('MDerk','var')
    MDerk = MovieData.load([erkExample filesep 'movieData.mat']);
end

imDirs = MDerk.getChannelPaths;
imNames = MDerk.getImageFileNames;

erkImage = imread([imDirs{iErkChan} filesep imNames{iErkChan}{1}]);
waveImage = imread([imDirs{iWaveChan} filesep imNames{iWaveChan}{1}]);

iWinProc = MDerk.getProcessIndex('WindowingProcess',1,0);
windows = MDerk.processes_{iWinProc}.loadChannelOutput(1);
iSampProc = MDerk.getProcessIndex('WindowSamplingProcess',1,0);
erkSamp = MDerk.processes_{iSampProc}.loadChannelOutput(iErkChan);
wavSamp = MDerk.processes_{iSampProc}.loadChannelOutput(iWaveChan);
sampRat = erkSamp.avg ./ wavSamp.avg;

%Make sure we have the same range
sampRat(sampRat>rawErkCaxis(2)) = rawErkCaxis(2);
sampRat(sampRat<rawErkCaxis(1)) = rawErkCaxis(1);

plotWindowsColormapped(windows,sampRat,erkRatCmap);
panelAxes = get(panelFig,'CurrentAxes');
axis ij,axis tight
set(panelAxes,'XTick',[]);
set(panelAxes,'YTick',[]);
set(panelAxes,'Color','w');
box on

xlim([0 size(erkImage,2)]),ylim([0 size(erkImage,1)])
plotScaleBar(scBarSz / MDerk.pixelSize_,'Location','Southeast','Color',scBarCol);

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end



%% ----- Erk Cell-Averaged Local Stuff Fig ---- %%
%Uses cells which are loaded and processed in prev figure - that cell must be run first.

%cellAvgFigDir = '/home/he19/files/LCCB/gtpases/Hunter/methods_paper_data/Michelle_IF_Erk_Wave_converted/HMEC use subset prominent lamellae';
cellAvgFigDir = '/home/he19/files/LCCB/gtpases/Hunter/methods_paper_data/Michelle_IF_Erk_Wave_converted/HMEC use';
cellAvgFile = 'combined sample analysis.mat';
%Reload the movieList because it doesnt load from saved file right


load([cellAvgFigDir filesep cellAvgFile]);

szLam = size(sampLamIntBothAvg);
szOth = size(sampOthIntBothAvg);
nChan = size(sampAllIntBothAvg,4);
nCell = size(sampAllIntBothAvg,3);

%Re-arrange samples so they agree with indices created above in stupid
%order
combSamp = cat(1,sampLamIntBothAvg,sampOthIntBothAvg);
combSamp = combSamp(:,:,iEx,:);
combSampEGF = combSamp(:,:,isEGF,:);
combSampNoEGF = combSamp(:,:,~isEGF,:);
combRegAvgEGF = squeeze(nanmean(combSamp(:,:,isEGF,:),3));
combRegAvgNoEGF = squeeze(nanmean(combSamp(:,:,~isEGF,:),3));
nShow = 18;

%% ---- Area labelling figure
%Make a stupid figure labelling lamellar and non-lamellar in combined
%circle map

panelName = 'Erk circle mapping region label';
panelFile = [figParentDir filesep panelName];
panelFig = figure;

areaLabel = zeros(size(combSamp(:,:,1,1)));
areaLabel(1:szLam(1),:) = 1;
areaLabel(szLam(1)+1:szOth(1)+szLam(1),:) = 2;

plotActivityMapPolar(areaLabel)
axis equal,axis off

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end




%% ----- Circle-Mapped EGF Cell Avg ----- %%

panelNameBase = 'Erk circle mapped regionally averaged EGF ';

useSat = satPct * 3;

subplot(1,nChan,1)
%nShow = 25;

cMap = jet;
mapRngEGF = zeros(nChan,2);

nCol = 512;
for j = 1:nChan
    
    
    panelName = [panelNameBase chanNames{j}];
    panelFile = [figParentDir filesep panelName];
    panelFig = figure;
    hold on
        
    tmp = combRegAvgEGF(:,:,j);
    
    tmp = smoothActivityMap(tmp(:,1:nShow),'FillNaN',true);
    
    plotActivityMapPolar(tmp)    
    
    if j == iWaveChan || j == iActChan%For some reason the actin looks like shit, so stretch range
        mapRngEGF(j,:) = prctile(tmp(:),[useSat*3 (100-useSat*4)]); 
    else
        mapRngEGF(j,:) = prctile(tmp(:),[useSat/2 (100-useSat/2)]);        
    end
    colormap(cMap);
    caxis(mapRngEGF(j,:))
    axis equal,axis off
    
    if saveFigs
        print(panelFig,panelFile,pOptTIFF{:});
        print(panelFig,panelFile,pOptEPS{:});
        hgsave(panelFig,panelFile);
    end    
        
    
    % -- Second version with fluorescence-style colormap
    if j <= 2
        %Matches other panel for ERK/WAVE
        cMapFluor = zeros(nCol,3);
        cMapFluor(:,j) = linspace(0,1,nCol);                
    elseif j == 4
        %DAPI is blue
        cMapFluor = zeros(nCol,3);
        cMapFluor(:,3) = linspace(0,1,nCol);
    elseif j == 3
        %Phalloidin as magenta....
        cMapFluor = zeros(nCol,3);
        cMapFluor(:,1) = linspace(0,1,nCol);
        cMapFluor(:,3) = linspace(0,1,nCol);
    end
        
        
    
    panelName = [panelNameBase chanNames{j} ' fluorescence colormap'];
    panelFile = [figParentDir filesep panelName];
    colormap(cMapFluor)
    if saveFigs
        print(panelFig,panelFile,pOptTIFF{:});
        print(panelFig,panelFile,pOptEPS{:});
        hgsave(panelFig,panelFile);
    end
    
    % --- Colorbar for fluorescence map


    panelName = [panelNameBase chanNames{j} ' fluorescence colormap colorbar'];
    panelFile = [figParentDir filesep panelName];
    panelFig = figure;

    colormap(cMapFluor);
    cBarAx = colorbar(cBarPars{:});
    axis off
    caxis(mapRngEGF(j,:));
    set(get(cBarAx,'YLabel'),'String',[chanNames{j} ' Fluorescence, a.u.'],cBarPars{:})

    if saveFigs
        print(panelFig,panelFile,pOptTIFF{:});
        print(panelFig,panelFile,pOptEPS{:});
        hgsave(panelFig,panelFile);
    end
    
    % --- Colorbar for jet map


    panelName = [panelNameBase chanNames{j} ' colorbar'];
    panelFile = [figParentDir filesep panelName];
    panelFig = figure;

    colormap(cMap);
    cBarAx = colorbar(cBarPars{:});
    axis off
    caxis(mapRngEGF(j,:));
    set(get(cBarAx,'YLabel'),'String',[chanNames{j} ' Fluorescence, a.u.'],cBarPars{:})

    if saveFigs
        print(panelFig,panelFile,pOptTIFF{:});
        print(panelFig,panelFile,pOptEPS{:});
        hgsave(panelFig,panelFile);
    end


    
end




%% ----- Circle-Mapped Unstimulated Cell Avg ----- %%
 
 panelNameBase = 'Erk circle mapped regionally averaged no EGF ';
useSat = satPct * 3;

subplot(1,nChan,1)
%nShow = 25;

cMap = jet;

for j = 1:nChan
    
    
    panelName = [panelNameBase chanNames{j}];
    panelFile = [figParentDir filesep panelName];
    panelFig = figure;
    hold on
        
    tmp = combRegAvgNoEGF(:,:,j);
    
    tmp = smoothActivityMap(tmp(:,1:nShow),'FillNaN',true);
    
    plotActivityMapPolar(tmp)    
    
    if j == iActChan %For some reason the actin looks like shit, so stretch range
        mapRng = prctile(tmp(:),[useSat*3 (100-useSat*4)]); 
    else
        mapRng = prctile(tmp(:),[useSat/2 (100-useSat/2)]);        
    end
    colormap(cMap);
    caxis(mapRng)
    axis equal,axis off
    
    if saveFigs
        print(panelFig,panelFile,pOptTIFF{:});
        print(panelFig,panelFile,pOptEPS{:});
        hgsave(panelFig,panelFile);
    end    
    
    
    %--- Second version with same range as EGF image
    
    caxis(mapRngEGF(j,:))
    
    panelName = [panelNameBase chanNames{j} ' using EGF colorscale'];
    panelFile = [figParentDir filesep panelName];
    
    if saveFigs
        print(panelFig,panelFile,pOptTIFF{:});
        print(panelFig,panelFile,pOptEPS{:});
        hgsave(panelFig,panelFile);
    end    
        
    % -- Second version with fluorescence-style colormap
    if j <= 2
        %Matches other panel for ERK/WAVE
        cMapFluor = zeros(nCol,3);
        cMapFluor(:,j) = linspace(0,1,nCol);                
    elseif j == 4
        %DAPI is blue
        cMapFluor = zeros(nCol,3);
        cMapFluor(:,3) = linspace(0,1,nCol);
    elseif j == 3
        %Phalloidin as magenta....
        cMapFluor = zeros(nCol,3);
        cMapFluor(:,1) = linspace(0,1,nCol);
        cMapFluor(:,3) = linspace(0,1,nCol);
    end
        
        
    
    panelName = [panelNameBase chanNames{j} ' fluorescence colormap'];
    panelFile = [figParentDir filesep panelName];
    colormap(cMapFluor)
    caxis(mapRng)
    if saveFigs
        print(panelFig,panelFile,pOptTIFF{:});
        print(panelFig,panelFile,pOptEPS{:});
        hgsave(panelFig,panelFile);
    end
    
    panelName = [panelNameBase chanNames{j} ' fluorescence colormap using EGF colorscale'];
    panelFile = [figParentDir filesep panelName];
    colormap(cMapFluor)
    caxis(mapRngEGF(j,:))
    if saveFigs
        print(panelFig,panelFile,pOptTIFF{:});
        print(panelFig,panelFile,pOptEPS{:});
        hgsave(panelFig,panelFile);
    end
    
    
    
    % --- Colorbar for above figure (with own color range)


    panelName = [panelNameBase chanNames{j} ' colorbar'];
    panelFile = [figParentDir filesep panelName];
    panelFig = figure;

    colormap(cMap);
    cBarAx = colorbar(cBarPars{:});
    axis off
    caxis(mapRng);
    set(get(cBarAx,'YLabel'),'String',[chanNames{j} ' Fluorescence, a.u.'],cBarPars{:})

    if saveFigs
        print(panelFig,panelFile,pOptTIFF{:});
        print(panelFig,panelFile,pOptEPS{:});
        hgsave(panelFig,panelFile);
    end
    
    % --- Colorbar for fluorescence map (own color range)


    panelName = [panelNameBase chanNames{j} ' fluorescence colormap colorbar'];
    panelFile = [figParentDir filesep panelName];
    panelFig = figure;

    colormap(cMapFluor);
    cBarAx = colorbar(cBarPars{:});
    axis off
    caxis(mapRng);
    set(get(cBarAx,'YLabel'),'String',[chanNames{j} ' Fluorescence, a.u.'],cBarPars{:})

    if saveFigs
        print(panelFig,panelFile,pOptTIFF{:});
        print(panelFig,panelFile,pOptEPS{:});
        hgsave(panelFig,panelFile);
    end

        

end




%% ----- Erk Circle Mapped Comparison Fig ---- %%


panelNameBase = 'Erk circle mapped comparison ';



for j = 1:nChan

    
    panelName = [panelNameBase chanNames{j}];
    panelFile = [figParentDir filesep panelName];
    panelFig = figure;
    hold on


    %Intepolate first, since this is just for visualization. Otherwise we get
    %too many NaN's in resulting diff map
    diffMap = smoothActivityMap(combRegAvgEGF(:,1:nShow,j),'FillNaN',true) ./ smoothActivityMap(combRegAvgNoEGF(:,1:nShow,j),'FillNaN',true);


    mapRng = prctile(diffMap(:),[satPct/2 (100-satPct/2)]);
    %mapRng = [-mapRng(2) mapRng(2)];%Use symmetric cmap so pos/neg is
    %clear (for difference, for ratio use regular saturation)
    
    
    %smMap = smoothActivityMap(diffMap,'FillNaN',true);

    plotActivityMapPolar(diffMap);
    caxis(mapRng);axis off,axis equal
    
    if j == 1
        erkDiffRng = mapRng;
    end

    if saveFigs
        print(panelFig,panelFile,pOptTIFF{:});
        print(panelFig,panelFile,pOptEPS{:});
        hgsave(panelFig,panelFile);
    end
    
    if j > 1
        %Second version with same map as Erk for simplicity        
        panelName = [panelNameBase chanNames{j} ' with Erk cAxis'];
        panelFile = [figParentDir filesep panelName];
        
        caxis(erkDiffRng)
        if saveFigs
            print(panelFig,panelFile,pOptTIFF{:});
            print(panelFig,panelFile,pOptEPS{:});
            hgsave(panelFig,panelFile);
        end
    end
        

    % --- Colorbar for above figure
    panelName = [panelNameBase chanNames{j} ' colorbar'];
    panelFile = [figParentDir filesep panelName];
    panelFig = figure;

    colormap(cMap);
    cBarAx = colorbar(cBarPars{:});
    axis off
    caxis(mapRng);
    set(get(cBarAx,'YLabel'),'String',['Fluorescence Ratio, +:- EGF'],cBarPars{:})

    if saveFigs
        print(panelFig,panelFile,pOptTIFF{:});
        print(panelFig,panelFile,pOptEPS{:});
        hgsave(panelFig,panelFile);
    end

        
end






%% ----- Line-Scans of Cell-Averaged ---- %%


%Output for EGF-only
cellAvgFigDir = '/home/he19/files/LCCB/gtpases/Hunter/methods_paper_data/Michelle_IF_Erk_Wave_converted/HMEC EGF use for figs';


%% ---- Lamellar

panelName = 'Erk line scan lamellar';
panelFile = [figParentDir filesep panelName];
%panelFig = fsFigure(.75);


panelFig = open([cellAvgFigDir filesep 'away from edge profile lamellar average interp both directions all channels.fig']);

% chanCols = jet(nChan);
% 
% sampLamMean = nanmean(sampLamIntBothAvg,3);
% for j = 1:nChan
%     plot(squeeze(mean(sampLamMean(:,:,:,j),1)),'Color',chanCols(j,:))    
% end
% legend(chanNames)
% 

%Set the line styles
panelAxes = get(panelFig,'CurrentAxes');
lineHans = get(panelAxes,'Children');
arrayfun(@(x)(set(x,plotPars{:})),lineHans);

xlabel('Distance From Cell Edge, Band #',axLabPars{:})
ylabel('Mean Fluorescence Intensity, a.u.',axLabPars{:})
title('');
set(panelAxes,axPars{:});
legend('hide')%turn legends off on these, we'll have one legend for all panels


if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end

%% ---- Non-Lamellar


panelName = 'Erk line scan non-lamellar';
panelFile = [figParentDir filesep panelName];
%panelFig = fsFigure(.75);

panelFig = open([cellAvgFigDir filesep 'away from edge profile oth average interp both directions all channels.fig']);
%Set the line styles
panelAxes = get(panelFig,'CurrentAxes');
lineHans = get(panelAxes,'Children');
arrayfun(@(x)(set(x,plotPars{:})),lineHans);

xlabel('Distance From Cell Edge, Band #',axLabPars{:})
ylabel('Mean Fluorescence Intensity, a.u.',axLabPars{:})
title('');
set(gca,axPars{:});

legend('hide')%turn legends off on these, we'll have one legend for all panels

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end

%% ---- Across-Lamellar


panelName = 'Erk line scan across lamellar';
panelFile = [figParentDir filesep panelName];
%panelFig = fsFigure(.75);

panelFig = open([cellAvgFigDir filesep 'along edge profile lamellar average interp both directions all channels.fig']);

%Set the line styles
panelAxes = get(panelFig,'CurrentAxes');
lineHans = get(panelAxes,'Children');
arrayfun(@(x)(set(x,plotPars{:})),lineHans);


xlabel('Distance Along Cell Edge, Strip #',axLabPars{:})
ylabel('Mean Fluorescence Intensity, a.u.',axLabPars{:})
title('');
set(gca,axPars{:});

legend('hide')%turn legends off on these, we'll have one legend for all panels

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end

%% ---- Across Non-Lamellar


panelName = 'Erk line scan across non-lamellar';
panelFile = [figParentDir filesep panelName];
%panelFig = fsFigure(.75);

panelFig = open([cellAvgFigDir filesep 'along edge profile oth average interp both directions all channels.fig']);
%Set the line styles
panelAxes = get(panelFig,'CurrentAxes');
lineHans = get(panelAxes,'Children');
arrayfun(@(x)(set(x,plotPars{:})),lineHans);


xlabel('Distance Along Cell Edge, Strip #',axLabPars{:})
ylabel('Mean Fluorescence Intensity, a.u.',axLabPars{:})
title('');
set(panelAxes,axPars{:});

legend('hide')%turn legends off on these, we'll have one legend for all panels

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end

%% ---- Legend for Erk-Linescan panels
%Save one figure separately with just legend

panelName = 'Erk line scan legend';
panelFile = [figParentDir filesep panelName];
%panelFig = fsFigure(.75);

panelFig = open([cellAvgFigDir filesep 'along edge profile oth average interp both directions all channels.fig']);
%Set the line styles
panelAxes = get(panelFig,'CurrentAxes');
lineHans = get(panelAxes,'Children');
set(panelAxes,axPars{:})
arrayfun(@(x)(set(x,plotPars{:})),lineHans);%Set plot params so legend shows up correctly
arrayfun(@(x)(delete(x)),lineHans);%Remove all the plot crap

title('');
axis off


if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end




%% %%%%%%%% ================ FIGURE 6 ======================== %%%%%%%% %%
%%%%%%%%%%%%%%% ARP2/3 & Protrusion FIG%%%%%%%%%%%%%%%%%%%%%%%%%%%


arpList = '/home/he19/files/LCCB/gtpases/Hunter/methods_paper_data/Kwonmoo_Arp3/movieListHaloOnly.mat';
if ~exist('MLarp','var')
    MLarp = MovieList.load(arpList,0);
end

scBarSz = 5e3; %Scale bar size in nm

iArpEx = 1;
szUse = 1000;
typeUse = 'protrusion_based';

iFrame = 1;%For example images/windows
nFrames = MLarp.movies_{iArpEx}.nFrames_;
tIntArpEx = MLarp.movies_{iArpEx}.timeInterval_;

fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];
winDir = [MLarp.movies_{iArpEx}.outputDirectory_ filesep 'windows_' num2str(szUse) 'nm_' typeUse];
winFiles = dir([winDir filesep '*.mat']);
windows = load([winDir filesep winFiles(iFrame).name]);
windows = windows.windows;


ccDataArpEx = load([MLarp.movies_{iArpEx}.outputDirectory_ filesep 'sample_crosscorrelation_' num2str(szUse) 'nm_' typeUse filesep 'temporal crosscorrelation channel 1.mat']);
[nStripCC,nBandCC,nLagTot,~] = size(ccDataArpEx.ccProtActPerWin);
nLags = (nLagTot-1)/2;

tDataArpCC = nLags*tIntArpEx:-tIntArpEx:-nLags*tIntArpEx;


tLagLim = [-100 100];%Use common time-lag limits;
corrLim = [-.5 .5];%Use common correlation value limits

%% ---------Arp2/3 per-window correlation curves ----- %%

panelName = 'Arp23 example per window correlation';
panelFile = [figParentDir filesep panelName];
panelFig = fsFigure(.5);
hold on;
panelAxes = get(panelFig,'CurrentAxes');    

ccPerWinCols = [0 0 1];

%tmpCC = reshape(ccPerWinCols,[nStripCC*nBandCC,nLagTot,2]);

warning('off','PLOTTRANSPARENT:NAN');

for j = 1:nStripCC
    
    for k = 1:3%nBandCC ONly show bands where correlation is meaningful, otherwise the figure just looks like noise
        
        if ~any(isnan(ccDataArpEx.ccProtActPerWin(j,k,:,1)))
            %plotTransparent(tDataArpCC,squeeze(ccDataArpEx.ccProtActPerWin(j,k,:,1)),.01*ones(1,nLagTot),ccPerWinCols,.05,0);
            
            %Works pretty well with only bands 1:3
            plotTransparent(tDataArpCC,squeeze(ccDataArpEx.ccProtActPerWin(j,k,:,1)),squeeze(ccDataArpEx.ccProtActPerWin(j,k,:,2)),ccPerWinCols,.0135,0);                                                            
            
            %SHows all bands, but looks like shit
            %plotTransparent(tDataArpCC,squeeze(ccDataArpEx.ccProtActPerWin(j,k,:,1)),.03,ccPerWinCols,.01,0);            
            
        end
        
    end
    
    
end

xlim(tLagLim)
%ylim(corrLim)
%Use diff ylim for this due to variation / noise
ylim([-.6 .6])

plot(xlim,[0 0],'k',plotPars{:});
plot([0 0],ylim,'k',plotPars{:});
xlabel('Delay, Seconds',axLabPars{:});
ylabel('Correlation',axLabPars{:});
set(panelAxes,axPars{:})
set(panelFig,'DefaultLineLineSmoothing','on');
set(panelFig,'DefaultPatchLineSmoothing','on');
box on



if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end


%% ------ Arp2/3 Spatial Variation in Correlation

panelName = 'Arp23 example correlation variation';
panelFile = [figParentDir filesep panelName];
panelFig = fsFigure(.5);
hold on;



%ccPerWinMat = max(abs(ccDataArpEx.ccProtActPerWin(:,:,:,1)),[],3);
%ccPerWinMat = mean(abs(ccDataArpEx.ccProtActPerWin(:,:,:,1)),3);
%ccPerWinMat = sum(abs(ccDataArpEx.ccProtActPerWin(:,:,:,1)),3);
%ccPerWinMat = sum(ccDataArpEx.ccProtActPerWin(:,:,:,1) .^2,3);
%ccPerWinMat =  mean(abs(ccDataArpEx.ccProtActPerWin(:,:,:,1)),3) ./ mean(abs(ccDataArpEx.ccProtActPerWin(:,:,:,2)),3);
%ccPerWinMat =  max(abs(ccDataArpEx.ccProtActPerWin(:,:,:,1)),[],3) ./ mean(abs(ccDataArpEx.ccProtActPerWin(:,:,:,2)),3);
%ccPerWinMat =  mean(ccDataArpEx.ccProtActPerWin(:,:,:,1) .^2,3) ./ mean(abs(ccDataArpEx.ccProtActPerWin(:,:,:,2)),3);
%ccPerWinMat =  mean(ccDataArpEx.ccProtActPerWin(:,:,:,1) .^2,3);

ccPerWinMat =  1 ./ mean(ccDataArpEx.ccProtActPerWin(:,:,:,2) ./ sqrt(nFrames) ,3);



mapRng = prctile(ccPerWinMat(~isnan(ccPerWinMat(:))),[satPct/2 100-(satPct/2)]);
ccPerWinMat(ccPerWinMat> mapRng(2)) = mapRng(2);
ccPerWinMat(ccPerWinMat< mapRng(1)) = mapRng(1);

%ccPerWinMat(isnan(ccPerWinMat)) = min(ccPerWinMat(:));

ccPerWinCmap = jet;
%ccPerWinCmap = ccPerWinCmap(end:-1:1,:);

%Transpose the windows since that's how we've been showing this image
%throughout the paper
windowsTP = cellfun(@(z)(cellfun(@(y)(cellfun(@(x)(x([2 1],:)),y,'Unif',0)),z,'Unif',0)),windows,'Unif',0);

plotWindows(windowsTP,{'k','EdgeColor','None'})
plotWindowsColormapped(windowsTP,ccPerWinMat,ccPerWinCmap);
box on
set(gca,'XTick',[])
set(gca,'YTick',[])
axis ij
set(gca,'XDir','reverse')
set(gca,'Color','w')
plotScaleBar(scBarSz ./ MLarp.movies_{iArpEx}.pixelSize_,'Color',zeros(1,3));
%title('Mean Variance of Cross-Correlation')

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end

%% ---- Colorbar for above figure

panelName = [panelName ' colorbar'];
panelFile = [figParentDir filesep panelName];
panelFig = figure;

colormap(ccPerWinCmap);
cBarAx = colorbar(cBarPars{:});
axis off
caxis(mapRng);
set(get(cBarAx,'YLabel'),'String','1 / (Mean SEM of Correlation, all Delays)',cBarPars{:})

if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end





%% ------- Arp2/3 Per-Cell and Cobined CC ------ %

panelName = 'arp23 combined cross correlation';
panelFile = [figParentDir filesep panelName];
panelFig = fsFigure(.5);


nMov = numel(MLarp.movies_);
sizeUse = 1e3;%Windowing result size to use.
%winUse = 'constant_number';
winUse = 'protrusion_based';
ccString = ['sample_crosscorrelation_' num2str(sizeUse) 'nm_' winUse];
ccFile = 'temporal crosscorrelation channel 1.mat';

iBand = 1;%Band to use to combine crosscorr

for j =  1:nMov        
    cc(j) = load([MLarp.movies_{j}.outputDirectory_ filesep ccString filesep ccFile]);              
    nFrames(j) = MLarp.movies_{j}.nFrames_;
    nCorr(j,:) = size(cc(j).bandMeanCC);      
    nDelay(j) = (nCorr(j,2)-1)/2;
    tInt(j) = MLarp.movies_{j}.timeInterval_;
    pSize(j) = MLarp.movies_{j}.pixelSize_;
    
    %Bootstrap the per-cell cc    
    [cellMeanCC{j},cellMeanCI{j}] = correlationBootstrap(squeeze(cc(j).ccProtActPerWin(:,iBand,:,1))',1.96 ./ sqrt(cc(j).nObsComb(:,iBand))');
    
end

if numel(unique(tInt)) > 1 %|| numel(unique(pSize)) > 1
    error('Inconsistent spatiotemporal resolution!!! Can''t combine crosscorr!')
else
    tInt = tInt(1);
end

nCorrUse = min(nCorr,[],1);
ccAllCell = cell(nMov,1);
cbAllCell = cell(nMov,1);
nLagsUse = (nCorrUse(2)-1)/2;
tData = nLagsUse*tInt:-tInt:-nLagsUse*tInt;

for j = 1:nMov
    
    lagDiff = nCorr(j,2) - nCorrUse(2);
    if lagDiff > 0
        ccBandMeanAllCell(:,j) = cc(j).bandMeanCC(iBand,lagDiff/2+1:end-(lagDiff/2));
        ccAllCell{j} = squeeze(cc(j).ccProtActPerWin(:,iBand,lagDiff/2+1:end-(lagDiff/2),1));             
    else
        ccAllCell{j} = squeeze(cc(j).ccProtActPerWin(:,iBand,:,1));
    end
    %tmp = cc(j).nObsComb(:,iBand);           
    cbAllCell{j} = 1.96 ./ sqrt(cc(j).nObsComb(:,iBand));
end

ccAllCell = cat(1,ccAllCell{:});
cbAllCell = cat(1,cbAllCell{:});



%Plot the per-cell curves first so they're behind
perCellCol = [0 0 0]+.4;
for j = 1:nMov
    currTData = nDelay(j)*tInt:-tInt:-nDelay(j)*tInt;
    plot(currTData,cellMeanCC{j},'color',perCellCol);
    %hold on
    %plotTransparent(currTData,cellMeanCC{j},.01,perCellCol,.3,0);
    patch([currTData currTData(end:-1:1)],[cellMeanCI{j}(1,:),cellMeanCI{j}(2,end:-1:1)],perCellCol,'EdgeColor','none','FaceAlpha',.3)
    hold on
    
end
    


%Bootstrap the combined CC curves
[bandMeanCC,bandBootCI] = correlationBootstrap(ccAllCell',cbAllCell');

plot(tData,bandMeanCC,plotPars{:})
hold on
%plot(tData,squeeze(bandBootCI(1,:)),'--b',plotPars{:})
%legend('Mean Cross-Correlation','95% Confidence Interval')
%plot(tData,squeeze(bandBootCI(2,:)),'--b',plotPars{:})

patch([tData tData(end:-1:1)],[bandBootCI(1,:),bandBootCI(2,end:-1:1)],[0 0 1],'EdgeColor','none','FaceAlpha',.4)

xlim(tLagLim)
ylim(corrLim)
plot(xlim,[0 0],'k',plotPars{:})
plot([0 0],ylim,'k',plotPars{:})
set(gca,axPars{:})
xlabel('Delay, Seconds',axLabPars{:})
ylabel('Correlation',axLabPars{:})


if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end

%% ------ Arp2/3 Sample size variation and CC ------- %%

%Do processing in separate cell 

sizesUse = [16 8 4 2 1 .5 .25]*1e3;
nSizes = numel(sizesUse);

for k = 1:nSizes

    winUse = 'protrusion_based';
    ccString = ['sample_crosscorrelation_' num2str(sizesUse(k)) 'nm_' winUse];
    ccFile = 'temporal crosscorrelation channel 1.mat';

    iBand = 1;


    for j =  1:nMov        
        cc(j) = load([MLarp.movies_{j}.outputDirectory_ filesep ccString filesep ccFile]);              
        nFrames(j) = MLarp.movies_{j}.nFrames_;
        nCorr(j,:) = size(cc(j).bandMeanCC);      
        nDelay(j) = (nCorr(j,2)-1)/2;
        tInt(j) = MLarp.movies_{j}.timeInterval_;
        pSize(j) = MLarp.movies_{j}.pixelSize_;

        %Bootstrap the per-cell cc    
        [cellMeanCC{j},cellMeanCI{j}] = correlationBootstrap(squeeze(cc(j).ccProtActPerWin(:,iBand,:,1))',1.96 ./ sqrt(cc(j).nObsComb(:,iBand))');

    end

    if numel(unique(tInt)) > 1 %|| numel(unique(pSize)) > 1
        error('Inconsistent spatiotemporal resolution!!! Can''t combine crosscorr!')
    else
        tInt = tInt(1);
    end

    nCorrUse = min(nCorr,[],1);
    ccAllCell = cell(nMov,1);
    cbAllCell = cell(nMov,1);
    nLagsUse = (nCorrUse(2)-1)/2;
    tData = nLagsUse*tInt:-tInt:-nLagsUse*tInt;

    for j = 1:nMov

        lagDiff = nCorr(j,2) - nCorrUse(2);
        if lagDiff > 0
            ccBandMeanAllCell(:,j) = cc(j).bandMeanCC(iBand,lagDiff/2+1:end-(lagDiff/2));
            ccAllCell{j} = squeeze(cc(j).ccProtActPerWin(:,iBand,lagDiff/2+1:end-(lagDiff/2),1));             
        else
            ccAllCell{j} = squeeze(cc(j).ccProtActPerWin(:,iBand,:,1));
        end
        %tmp = cc(j).nObsComb(:,iBand);           
        cbAllCell{j} = 1.96 ./ sqrt(cc(j).nObsComb(:,iBand));
    end

    ccAllCell = cat(1,ccAllCell{:});
    cbAllCell = cat(1,cbAllCell{:});


    %Bootstrap the combined CC curves
    [allCellBandMeanCC{k},allCellBandMeanCI{k}] = correlationBootstrap(ccAllCell',cbAllCell');
end

%% 


panelName = 'arp23 oversized sample cross correlation';
panelFile = [figParentDir filesep panelName];
panelFig = fsFigure(.5);
hold on
sizeCols = isomorphicColormap('b',nSizes);
%sizeCols = jet(nSizes);

legStr = cell(nSizes,1);
% for k = 1:nSizes
% 
%     plot(tData,allCellBandMeanCC{k},'color',sizeCols(k,:),plotPars{:})
%     hold on
%     %plot(tData,squeeze(allCellBandMeanCI{k}(1,:)),'--b',plotPars{:})
%     %legend('Mean Cross-Correlation','95% Confidence Interval')
%     %plot(tData,squeeze(allCellBandMeanCI{k}(2,:)),'--b',plotPars{:})
% 
%     patch([tData tData(end:-1:1)],[allCellBandMeanCI{k}(1,:),allCellBandMeanCI{k}(2,end:-1:1)],sizeCols(k,:),'EdgeColor','none','FaceAlpha',.2)
%  
%     legStr{(k-1)*2+1} = [num2str(sizesUse(k)/1e3) 'um'];
%     legStr{(k-1)*2+2} = [num2str(sizesUse(k)/1e3) 'um C.I.'];
% 
% end

%Plot the CIs first so they're in the back and its less confusing
for k = 1:nSizes

    
    patch([tData tData(end:-1:1)],[allCellBandMeanCI{k}(1,:),allCellBandMeanCI{k}(2,end:-1:1)],sizeCols(k,:),'EdgeColor','none','FaceAlpha',1)
    
    %plot(tData,allCellBandMeanCC{k},'color',sizeCols(k,:)/1.2)
 
    legStr{k} = [num2str(sizesUse(k)) 'nm'];        
end

legend(legStr)

% for k = 1:nSizes
%     
%     plot(tData,allCellBandMeanCC{k},'color',sizeCols(k,:),plotPars{:})
%     
%     
%     %plot(tData,squeeze(allCellBandMeanCI{k}(1,:)),'--b',plotPars{:})
%     %legend('Mean Cross-Correlation','95% Confidence Interval')
%     %plot(tData,squeeze(allCellBandMeanCI{k}(2,:)),'--b',plotPars{:})
% 
% end



xlim(tLagLim)
ylim(corrLim)
plot(xlim,[0 0],'k',plotPars{:})
plot([0 0],ylim,'k',plotPars{:})
set(gca,axPars{:})
box on
xlabel('Delay, Seconds',axLabPars{:})
ylabel('Correlation',axLabPars{:})


if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end



%% --------------- Control Halo & Protrusion --------- %%

haloList = '/home/he19/files/LCCB/gtpases/Hunter/methods_paper_data/Control_Halo_TMR_12262011/Cropped/movieList.mat';

if ~exist('MLhalo','var')
    MLhalo = MovieList.load(haloList,0);
end

%% ---- Halo combined CC

panelName = 'halo only combined cross correlation';
panelFile = [figParentDir filesep panelName];
panelFig = fsFigure(.5);


nMov = numel(MLhalo.movies_);
sizeUse = 1e3;%Windowing result size to use.
%winUse = 'constant_number';
winUse = 'protrusion_based';
ccString = ['sample_crosscorrelation_' num2str(sizeUse) 'nm_' winUse];
ccFile = 'temporal crosscorrelation channel 1.mat';

iBand = 1;%Band to use to combine crosscorr

for j =  1:nMov        
    cc(j) = load([MLhalo.movies_{j}.outputDirectory_ filesep ccString filesep ccFile]);              
    nFrames(j) = MLhalo.movies_{j}.nFrames_;
    nCorr(j,:) = size(cc(j).bandMeanCC);      
    nDelay(j) = (nCorr(j,2)-1)/2;
    tInt(j) = MLhalo.movies_{j}.timeInterval_;
    pSize(j) = MLhalo.movies_{j}.pixelSize_;
    
    %Bootstrap the per-cell cc    
    [cellMeanCC{j},cellMeanCI{j}] = correlationBootstrap(squeeze(cc(j).ccProtActPerWin(:,iBand,:,1))',1.96 ./ sqrt(cc(j).nObsComb(:,iBand))');
    
end

if numel(unique(tInt)) > 1 %|| numel(unique(pSize)) > 1
    error('Inconsistent spatiotemporal resolution!!! Can''t combine crosscorr!')
else
    tInt = tInt(1);
end

nCorrUse = min(nCorr,[],1);
ccAllCell = cell(nMov,1);
cbAllCell = cell(nMov,1);
nLagsUse = (nCorrUse(2)-1)/2;
tData = nLagsUse*tInt:-tInt:-nLagsUse*tInt;

for j = 1:nMov
    
    lagDiff = nCorr(j,2) - nCorrUse(2);
    if lagDiff > 0
        ccBandMeanAllCell(:,j) = cc(j).bandMeanCC(iBand,lagDiff/2+1:end-(lagDiff/2));
        ccAllCell{j} = squeeze(cc(j).ccProtActPerWin(:,iBand,lagDiff/2+1:end-(lagDiff/2),1));             
    else
        ccAllCell{j} = squeeze(cc(j).ccProtActPerWin(:,iBand,:,1));
    end
    %tmp = cc(j).nObsComb(:,iBand);           
    cbAllCell{j} = 1.96 ./ sqrt(cc(j).nObsComb(:,iBand));
end

ccAllCell = cat(1,ccAllCell{:});
cbAllCell = cat(1,cbAllCell{:});



%Plot the per-cell curves first so they're behind
perCellCol = [0 0 0]+.4;
for j = 1:nMov
    currTData = nDelay(j)*tInt:-tInt:-nDelay(j)*tInt;
    plot(currTData,cellMeanCC{j},'color',perCellCol);
    %hold on
    %plotTransparent(currTData,cellMeanCC{j},.01,perCellCol,.3,0);
    patch([currTData currTData(end:-1:1)],[cellMeanCI{j}(1,:),cellMeanCI{j}(2,end:-1:1)],perCellCol,'EdgeColor','none','FaceAlpha',.3)
    hold on
    
end
    


%Bootstrap the combined CC curves
[bandMeanCC,bandBootCI] = correlationBootstrap(ccAllCell',cbAllCell');

plot(tData,bandMeanCC,plotPars{:})
hold on
%plot(tData,squeeze(bandBootCI(1,:)),'--b',plotPars{:})
%legend('Mean Cross-Correlation','95% Confidence Interval')
%plot(tData,squeeze(bandBootCI(2,:)),'--b',plotPars{:})

patch([tData tData(end:-1:1)],[bandBootCI(1,:),bandBootCI(2,end:-1:1)],[0 0 1],'EdgeColor','none','FaceAlpha',.4)

xlim(tLagLim)
ylim(corrLim)
plot(xlim,[0 0],'k',plotPars{:})
plot([0 0],ylim,'k',plotPars{:})
set(gca,axPars{:})
xlabel('Delay, Seconds',axLabPars{:})
ylabel('Correlation',axLabPars{:})


if saveFigs
    print(panelFig,panelFile,pOptTIFF{:});
    print(panelFig,panelFile,pOptEPS{:});
    hgsave(panelFig,panelFile);
end





