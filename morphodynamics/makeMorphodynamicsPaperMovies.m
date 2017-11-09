%% -------- Common Parameters ----------- %%


%movieDir = '/home/he19/orchestra/home/Papers/windowing methods paper/Movies';
movieDir = '/home/he19/Desktop/TEMP/temp_morpho_mov_writing';%Workaround problems writing to SMB mounted nucleus drive...

makeAvi = true;
makeMov = true;
movQual = .99;

satPct = 1;

iArpEx = 1;

%% %%%%%%%% ================ Movie 1 ======================== %%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Arp3 cell with edge segmentation overlay




%Use example Arp2/3 cell per G's request
arpList = '/home/he19/files/LCCB/gtpases/Hunter/methods_paper_data/Kwonmoo_Arp3/movieListHaloOnly.mat';
if ~exist('MLarp','var')
    MLarp = MovieList.load(arpList,0);
end

scBarSz = 5e3; %Scale bar size in nm

tInt = MLarp.movies_{iArpEx}.timeInterval_;

%% ------Arp  Example w/ Edge Seg ----- %%

movieName = 'Movie 1 - Arp3 example with edge segmentation';
movieFile = [movieDir filesep movieName];
movieFig = figure;
axis;
movieAxes = get(movieFig,'CurrentAxes');


pixSize = MLarp.movies_{iArpEx}.pixelSize_;
nFrames = MLarp.movies_{iArpEx}.nFrames_;

prot = MLarp.movies_{iArpEx}.processes_{MLarp.movies_{iArpEx}.getProcessIndex('ProtrusionProcess',1,0)}.loadChannelOutput;

%I just re-wrote the movie-making explicitly here so it would match the
%figure....


for iFrame = 1:141%nFrames

    %MLarp.movies_{iArpEx}.channels_(1).draw(iFrame,'hAxes',panelAxes);
    %Transpose so it fits in the figure better
    im = MLarp.movies_{iArpEx}.channels_(1).loadImage(iFrame);
    cla(movieAxes);
    imshow(im',[]);
    hold on
    saturateImageColormap(movieAxes,satPct*2);
    plotScaleBar(scBarSz / pixSize,'Handle',movieAxes,'Location','SouthWest','Label',['t=' num2str( (iFrame-1)*tInt ) 's']);
    set(movieAxes,'XDir','reverse');
          
    plot(prot.smoothedEdge{iFrame}(:,2),prot.smoothedEdge{iFrame}(:,1),'r','LineWidth',3);
                
    
    if makeAvi
        movieFrames(iFrame) = getframe(movieFig);
    end
    if makeMov
        if iFrame == 1
            MakeQTMovie('start',[movieFile '.mov']);
            MakeQTMovie('quality',movQual);
        end
        MakeQTMovie('addfigure');
    end
            
        
end

if makeMov
    MakeQTMovie('finish');
end
if makeAvi            
    movie2avi(movieFrames,[movieFile '.avi']);            
end

%% %%%%%%%% ================ Movie 2 ======================== %%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Arp3 cell with windowing overlay


MDwin = MLarp.movies_{iArpEx};
exampMovWin =  MLarp.movies_{iArpEx}.outputDirectory_;
winSize = 2e3;
%winType = 'constant_number';
winType = 'protrusion_based';
pixSize = MDwin.pixelSize_;

winDir = [exampMovWin filesep 'windows_' num2str(winSize) 'nm_' winType ];
winFiles = dir([winDir filesep '*.mat']);

scBarSz = 5e3; %Scale bar size in nm



movieName = 'Movie 2 - Arp3 example with protrusion based windowing';
movieFile = [movieDir filesep movieName];
movieFig = figure;
axis;
movieAxes = get(movieFig,'CurrentAxes');


pixSize = MLarp.movies_{iArpEx}.pixelSize_;
nFrames = MLarp.movies_{iArpEx}.nFrames_;

prot = MLarp.movies_{iArpEx}.processes_{MLarp.movies_{iArpEx}.getProcessIndex('ProtrusionProcess',1,0)}.loadChannelOutput;

%I just re-wrote the movie-making explicitly here so it would match the
%figure....


movieFrames(1:nFrames) = struct('cdata',[],'colormap',[]);

for iFrame = 1:141%nFrames

    %MLarp.movies_{iArpEx}.channels_(1).draw(iFrame,'hAxes',panelAxes);
    %Transpose so it fits in the figure better
    im = MLarp.movies_{iArpEx}.channels_(1).loadImage(iFrame);
    cla(movieAxes);
    imshow(im',[]);
    hold on
    saturateImageColormap(movieAxes,satPct*2);
    plotScaleBar(scBarSz / pixSize,'Handle',movieAxes,'Location','SouthWest','Label',['t=' num2str( (iFrame-1)*tInt ) 's']);
    
          
    %plot(prot.smoothedEdge{iFrame}(:,2),prot.smoothedEdge{iFrame}(:,1),'r');
    
    windows = load([winDir filesep winFiles(iFrame).name]);
    windows = windows.windows;

    %Transpose the windows to match the image.
    windowsTP = cellfun(@(z)(cellfun(@(y)(cellfun(@(x)(x([2 1],:)),y,'Unif',0)),z,'Unif',0)),windows,'Unif',0);
           
    plotWindows(windowsTP,[{'r','FaceAlpha',0,'FaceColor','none','EdgeColor','r'} ,'LineWidth',3])
    
    set(movieAxes,'XDir','reverse');
    
    %---get the image
    if makeAvi
        movieFrames(iFrame) = getframe(movieFig);
    end
    if makeMov
        if iFrame == 1
            MakeQTMovie('start',[movieFile '.mov']);
            MakeQTMovie('quality',movQual);
        end
        MakeQTMovie('addfigure');
    end
            
        
end

if makeMov
    MakeQTMovie('finish');
end
if makeAvi            
    movie2avi(movieFrames,[movieFile '.avi']);            
end

%% %%%%%%%% ================ Movie 3 ======================== %%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rac1 whole cell with highlighted window - prot based

scBarSz = 1e4;%Scale bar size in nm
winSize = 1.6e3;%Arbitratily selected for visualization. Gives 5x5 pixel windows.

%Use highly motile rac1 cell from marco
wcExampMov = '/home/he19/files/LCCB/gtpases/Hunter/methods_paper_data/rac1_cell_actuallyMoves_from_marco';
if ~exist('MDwc','var')
    MDwc = MovieData.load([wcExampMov filesep 'movieData.mat']);
end
winType1 = 'protrusion_based';
winDir1 = [wcExampMov filesep 'windows_' num2str(winSize) 'nm_' winType1 ];
winFiles1 = dir([winDir1 filesep '*.mat']);

protVecs = MDwc.processes_{MDwc.getProcessIndex('ProtrusionProcess',1,0)}.loadChannelOutput;

showFrames = 1:1:100;

tInt = MDwc.timeInterval_;

%Axis lmits so Gaudenz doesn't cry about the extra white space
wcFigXlim = [50 260];
wcFigYlim = [53 351];

movieName = 'Movie 3 - Rac1 example with protrusion based window highlighted';
movieFile = [movieDir filesep movieName];
movieFig = figure('Position',[ 892  236  1032  1035]);
axis;

movieAxes = get(movieFig,'CurrentAxes');



nShow = numel(showFrames);

iWin = [29 1 ;
        70 1];    
    
nWinShow = size(iWin,1);    
clear movieFrames

edgeCol = [1 0 0 ];
winCol = [1 0 0 ];

for iFrame = 1:nShow
    
    
    cla(movieAxes);
    
    plot(protVecs.smoothedEdge{showFrames(iFrame)}(:,1),protVecs.smoothedEdge{showFrames(iFrame)}(:,2),'Color','k','LineWidth',2)
    
    windows = load([winDir1 filesep winFiles1(showFrames(iFrame)).name]);
    windows = windows.windows;
    for k = 1:nWinShow
        if ~isempty(windows{iWin(k,1)}) && ~isempty(windows{iWin(k,1)}{iWin(k,2)})
            plotWindows(windows{iWin(k,1)}{iWin(k,2)},[{'r','FaceAlpha',.6,'EdgeColor',edgeCol,'FaceColor',winCol,'LineWidth',3}]);        
        else
            disp(['no win 1 frame ' num2str(showFrames(iFrame))])
        end
        
        %plotWindows(windows,[{'r','FaceAlpha',1,'EdgeColor',winCols1(j,:),'FaceColor','none'} ]);        
    end
    
    xlim(wcFigXlim)
    ylim(wcFigYlim)

    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    box on
    
    plotScaleBar(scBarSz / MDwc.pixelSize_,'Location','Southwest','Color',zeros(1,3),'Label',['t=' num2str( (iFrame-1)*tInt ) 's']);
    
    
    
    %---get the image
    if makeAvi
        movieFrames(iFrame) = getframe(movieFig);
    end
    if makeMov
        if iFrame == 1
            MakeQTMovie('start',[movieFile '.mov']);
            MakeQTMovie('quality',movQual);
        end
        MakeQTMovie('addfigure');
    end
    
    
    
end


if makeMov
    MakeQTMovie('finish');
end
if makeAvi            
    movie2avi(movieFrames,[movieFile '.avi']);            
end

%% %%%%%%%% ================ Movie 4 ======================== %%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rac1 whole cell with highlighted window - constant number

scBarSz = 1e4;%Scale bar size in nm
winSize = 1.6e3;%Arbitratily selected for visualization. Gives 5x5 pixel windows.

%Use highly motile rac1 cell from marco
wcExampMov = '/home/he19/files/LCCB/gtpases/Hunter/methods_paper_data/rac1_cell_actuallyMoves_from_marco';
if ~exist('MDwc','var')
    MDwc = MovieData.load([wcExampMov filesep 'movieData.mat']);
end
winType1 = 'constant_number';
winDir1 = [wcExampMov filesep 'windows_' num2str(winSize) 'nm_' winType1 ];
winFiles1 = dir([winDir1 filesep '*.mat']);

protVecs = MDwc.processes_{MDwc.getProcessIndex('ProtrusionProcess',1,0)}.loadChannelOutput;

showFrames = 1:1:100;



%Axis lmits so Gaudenz doesn't cry about the extra white space
wcFigXlim = [50 260];
wcFigYlim = [53 351];

movieName = 'Movie 3 - Rac1 example with constant number window highlighted';
movieFile = [movieDir filesep movieName];
movieFig = figure('Position',[ 892  236  1032  1035]);
axis;

movieAxes = get(movieFig,'CurrentAxes');



nShow = numel(showFrames);

iWin = [29 1 ;
        70 1];    
    
nWinShow = size(iWin,1);    
clear movieFrames

edgeCol = [0 0 1 ];
winCol = [0 0 1 ];

for iFrame = 1:nShow
    
    
    cla(movieAxes);
    
    plot(protVecs.smoothedEdge{showFrames(iFrame)}(:,1),protVecs.smoothedEdge{showFrames(iFrame)}(:,2),'Color','k','LineWidth',2)
    
    windows = load([winDir1 filesep winFiles1(showFrames(iFrame)).name]);
    windows = windows.windows;
    for k = 1:nWinShow
        if ~isempty(windows{iWin(k,1)}) && ~isempty(windows{iWin(k,1)}{iWin(k,2)})
            plotWindows(windows{iWin(k,1)}{iWin(k,2)},[{'r','FaceAlpha',.6,'EdgeColor',edgeCol,'FaceColor',winCol,'LineWidth',3}]);        
        else
            disp(['no win 1 frame ' num2str(showFrames(iFrame))])
        end
        
        %plotWindows(windows,[{'r','FaceAlpha',1,'EdgeColor',winCols1(j,:),'FaceColor','none'} ]);        
    end
    
    xlim(wcFigXlim)
    ylim(wcFigYlim)

    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    box on
    
    plotScaleBar(scBarSz / MDwc.pixelSize_,'Location','Southwest','Color',zeros(1,3),'Label',['t=' num2str( (showFrames(iFrame)-1)*tInt ) 's']);
    
    
    
    %---get the image
    if makeAvi
        movieFrames(iFrame) = getframe(movieFig);
    end
    if makeMov
        if iFrame == 1
            MakeQTMovie('start',[movieFile '.mov']);
            MakeQTMovie('quality',movQual);
        end
        MakeQTMovie('addfigure');
    end
    
    
    
end


if makeMov
    MakeQTMovie('finish');
end
if makeAvi            
    movie2avi(movieFrames,[movieFile '.avi']);            
end


%% %%%%%%%% ================ Movie 5 ======================== %%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rac1 whole cell circle mapping side-by-side


scBarSz = 2.5e4;%Scale bar size in nm



movieName = 'Movie 4 - Rac1 example circle mapping';
movieFile = [movieDir filesep movieName];
movieFig = figure('Position',[560         219        1346         977]);


wcCmap = jet;


winType2 = 'constant_number';
sampDir2 = [wcExampMov filesep 'window_samples_' num2str(winSize) 'nm_' winType2 ];
sampFiles2 = dir([sampDir2 filesep '*.mat']);
samps2 = load([sampDir2 filesep sampFiles2.name]);
samps2 = samps2.samples;

sampSiz = size(samps2.avg);

nBandShow = size(samps2.avg,2);

for iFrame = 1:nShow
    
    
    % ---- Draw the ratio image
    
    subplot(1,2,1);
    panelAxes = get(movieFig,'CurrentAxes');
    cla(panelAxes);
    set(panelAxes,'Units','normalized')
    set(panelAxes,'Position',[.05 .05 .43 .90])    
    hold on    
   
    im = MDwc.channels_(1).loadImage(showFrames(iFrame));
    im = im ./ mean(im(:));
    imHan = imshow(im,[]);    
    panelAxes = get(imHan,'Parent');
    %panelFig = get(panelAxes,'Parent');
    colormap(wcCmap)        
    saturateImageColormap(imHan,satPct)    
    wcCaxis = caxis;        
    xlim(wcFigXlim)
    ylim(wcFigYlim)
    plotScaleBar(scBarSz / MDwc.pixelSize_,'Location','Northeast','Color',zeros(1,3),'Label',['t=' num2str((showFrames(iFrame)-1) * MDwc.timeInterval_) 's']);
    im = get(imHan,'CData');    
    tmpMask = imfill(im > 0,'holes');
    tmpMask = bwareaopen(tmpMask,100);
    set(imHan,'AlphaData',tmpMask)    
    axis on
    box on
    set(panelAxes,'Color','w');
    set(panelAxes,'XTick',[]);
    set(panelAxes,'YTick',[]);

    
    % ---- Draw the circle map
    
    subplot(1,2,2);
    panelAxes = get(movieFig,'CurrentAxes');
    cla(panelAxes);
    set(panelAxes,'Units','normalized')
    set(panelAxes,'Position',[.5 .05 .43 .90])
    hold on
    
    currSamp = squeeze(samps2.avg(:,1:nBandShow,showFrames(iFrame)));
    smMap = smoothActivityMap(currSamp,'FillNaN',true);        
        
    sampRng = prctile(currSamp(:),[satPct (100-satPct)]);            
    plotActivityMapPolar(smMap);           
    caxis(sampRng)
    axis image
    axis off

    %So it matches with image
    set(panelAxes,'XDir','reverse')
    set(panelAxes,'YDir','reverse')
    
    
    %---get the image
    if makeAvi
        movieFrames(iFrame) = getframe(movieFig);
    end
    if makeMov
        if iFrame == 1
            MakeQTMovie('start',[movieFile '.mov']);
            MakeQTMovie('quality',movQual);
        end
        MakeQTMovie('addfigure');
    end
    
    
end


if makeMov
    MakeQTMovie('finish');
end
if makeAvi            
    movie2avi(movieFrames,[movieFile '.avi']);            
end





