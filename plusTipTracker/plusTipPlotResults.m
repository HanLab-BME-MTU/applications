function plusTipPlotResults(projData,remBegEnd,timeRange,speedLim,lifeLim,dispLim,saveDir,useFirstImage)
% plusTipPlotResults make 18 summary overlay plots and histograms0
%
% SYNOPSIS: plusTipPlotResults(projData,remBegEnd,timeRange,speedLim,lifeLim,dispLim,saveDir)
%
% INPUT: projData: output of post-processing
%        remBegEnd: 1 to remove data from beginning/end of movie, 0 to
%                   retain it
%        timeRange: [startFrame endFrame] for plots
%        speedLim:  [min max] speed range (microns/min)
%        lifeLim:   [min max] lifetime range (seconds)
%        dispLim:   [min max] displacement range (microns)
%        saveDir:   path to output directory
% 
% OUTPUT: 
%       speed/lifetime/displacement overlays for growth, fgap, and bgap
%             sub-tracks (9)
%       speed/lifetime/displacement histograms showing stacked bars for
%             growth, fgap, and bgap sub-tracks (3)
%       initition/termination locations for fgaps and bgaps, separate and
%             merged (6)

fileExt='.tif';

if nargin<1 || isempty(projData)
    dirName=uigetdir(pwd,'Please select project directory');
    temp=load([dirName filesep 'meta' filesep 'projData.mat']);
    projData=temp.projData;
end

% assume we should plot all data
if nargin<2 || isempty(remBegEnd)
    remBegEnd=0;
end

% get number of time points
nFrames=projData.nFrames;
% check whether a time range for plotting was input
if nargin<3 || isempty(timeRange)
    timeRange=[1 nFrames];
else
    if timeRange(1)<1
        timeRange(1)=1;
    end
    if timeRange(2)>nFrames
        timeRange(2)=nFrames;
    end
end

% format image/analysis directory paths
imDir=formatPath(projData.imDir);
anDir=formatPath(projData.anDir);

if nargin<7 || isempty(saveDir)
    saveDir=uigetdir(anDir,'Please choose output directory');
end

if feature('ShowFigureWindows'), 
    wtBar = waitbar(0,'Initializing...');
else
    wtBar=-1;
end

% extract track info (block matrix)
allData=abs(projData.mergedDataMatAllSubTracksConverted);


if nargin<4 || isempty(speedLim) || strcmpi(speedLim,'max')
    speedLim=prctile(allData(:,4),95);    
end

if nargin<5 || isempty(lifeLim) || strcmpi(lifeLim,'max')
    lifeLim=prctile(allData(:,6),95);
end

if nargin<6 || isempty(dispLim) || strcmpi(dispLim,'max')
    dispLim=prctile(allData(:,7),95);
end

if nargin<8, useFirstImage=0; end

% get first image from imDir
if useFirstImage
    [listOfImages] = searchFiles('.tif',[],imDir,0);
    img = double(imread([listOfImages{1,2} filesep listOfImages{1,1}]));
else
    [filename pathname]  = uigetfile('*','Select an Image File',imDir);
    img = double(imread([pathname filesep filename]));
end

% get track type and matrices containing xy-coordinates for all subtracks
trackType=allData(:,5);
trackStarts=allData(:,2);
trackEnds=allData(:,3);
% convert to each subtrack
[xMat,yMat]=plusTipGetSubtrackCoords(projData,[],1);


if remBegEnd==1
    % any track not entirely contained within the frame range will be
    % excluded
    gIdx=find(trackType==1 & trackStarts>timeRange(1) & trackEnds<timeRange(2)); % just growth
    fIdx=find(trackType==2 & trackStarts>timeRange(1) & trackEnds<timeRange(2)); % just fgap
    bIdx=find(trackType==3 & trackStarts>timeRange(1) & trackEnds<timeRange(2)); % just bgap
else
    % any track which ends before the frame range begins or begins after
    % the frame range ends will be excluded
    gIdx=find(trackType==1 & trackStarts<timeRange(2) & trackEnds>timeRange(1)); % just growth
    fIdx=find(trackType==2 & trackStarts<timeRange(2) & trackEnds>timeRange(1)); % just fgap
    bIdx=find(trackType==2 & trackStarts<timeRange(2) & trackEnds>timeRange(1)); % just bgap
end

% get size of the images
minY=1; maxY=size(img,1);
minX=1; maxX=size(img,2);
scrsz = get(0,'ScreenSize');
screenW=scrsz(3);
screenL=scrsz(4);
magCoef=inf;
maxMagCoefW = (0.8*screenW)/(maxX-minX+1);
maxMagCoefL = (0.8*screenL)/(maxY-minY+1);
if magCoef > min([maxMagCoefW; maxMagCoefL])
    calcMagCoef = min([magCoef; maxMagCoefW; maxMagCoefL]);
else
    calcMagCoef = magCoef;
end
movieL = (calcMagCoef*(maxY-minY+1));
movieW = (calcMagCoef*(maxX-minX+1));
figPos=[round(screenW*(1-movieW/screenW)/2) round(screenL*(1-movieL/screenL)/2) movieW movieL];


overlayFig = figure;
imshow(img,[]);

set(gca,'Position',[0 0 1 .95]);
hold on

cMapLength=128; cMap=jet(cMapLength);

% Create MT dynamics maps of speed, lifetime and displacement for growth,
% fgaps and bgaps events
type= {'speed','lifetime','displacement'};
typeIndx= {4,6,7};
typeMatIndx = {1:3,4:6,7:9};
typeLim = {speedLim;lifeLim;dispLim};
typeUnits = {'microns/min';'s';'microns'};
event= {'growth tracks','fgaps','bgaps'};
eventIndx= {gIdx;fIdx;bIdx};


for i=1:3
    x=xMat(eventIndx{i},1:end-1)';
    y=yMat(eventIndx{i},1:end-1)';
    
    for j=1:3
        % SPEED
        data=allData(eventIndx{i},typeIndx{j});
        dataRange=ceil(typeLim{j});
        data(data>dataRange)=dataRange;
        mapper=linspace(0,dataRange,cMapLength)';
        
        % get closest colormap index for each feature
        D=createDistanceMatrix(data,mapper);
        [sD,idx]=sort(abs(D),2);
        
        % make the plot
        figure(overlayFig);
        delete(findobj(overlayFig,'Type','line'));
        for k=1:cMapLength
            plot(x(:,idx(:,1)==k),y(:,idx(:,1)==k),'color',cMap(k,:),'lineWidth',1);
        end
        title([type{j} ' map of ' event{i} ', max ' type{j}...
            ' ' num2str(dataRange) ' ' typeUnits{j}]);
        
        saveas(gcf,[saveDir filesep 'overlay_' event{i} '_' type{j} '.fig'])
        %saveas(gcf,[saveDir filesep 'overlay_' event{i} '_' type{j} filext])
        saveas(gcf,[saveDir filesep 'overlay_' event{i} '_' type{j} '.eps'],'psc2'); 
        if ishandle(wtBar)
            waitbar(((i-1)*3+ j)/18,wtBar,'Creating overlay maps');
        end
        
    end
end
delete(findobj(overlayFig,'Type','line'));
%%


%% Plot only Coord in Region 



if isfield(projData, 'xCoordInAllTracks')
    dataIn(:,1) = projData.speedInMicPerMinAllTracks; 
dataIn(:,2) = projData.insideSecAllTracks; 
roiMask = double(imread([projData.anDir filesep 'roiMask.tif'])); 
maskImg = img.*roiMask; 
%[y1,x1] = ind2sub(size(roiMask),find(roiMask,1)); 

allBCoordsYX = bwboundaries(roiMask); 

% always use all for now 
    maskFig = figure; 
    imshow(maskImg,[]); 
    hold on; 
    
    set(gca,'Position',[0 0 1 .95]);
hold on
    for j = 1:2 
        xMat = projData.xCoordInAllTracks; 
yMat = projData.yCoordInAllTracks; 

 x=xMat(:,1:end-1)';
    y=yMat(:,1:end-1)';
      data = dataIn(:,j); 
      dataRange = ceil(typeLim{j}); 
      
data(data>dataRange)=dataRange;
        mapper=linspace(0,dataRange,cMapLength)';
        
        % get closest colormap index for each feature
        D=createDistanceMatrix(data,mapper);
        [sD,idx]=sort(abs(D),2);
        
        % make the plot
        figure(maskFig);
        delete(findobj(maskFig,'Type','line'));
        for k=1:cMapLength
            plot(x(:,idx(:,1)==k),y(:,idx(:,1)==k),'color',cMap(k,:),'lineWidth',1);
        end
        
        if j == 1
            hold on
            dataOut = projData.speedOutMicPerMinAllTracks; 
           % dataOut = dataOut(~isnan(dataOut)); 
            dataOut(dataOut>dataRange) = dataRange; 
            
             mapper=linspace(0,dataRange,cMapLength)';
        
        % get closest colormap index for each feature
        D=createDistanceMatrix(dataOut,mapper);
        [sD,idx]=sort(abs(D),2);
            
            
            xMat = projData.xCoordOutAllTracks; 
            yMat = projData.yCoordOutAllTracks; 
            x = xMat(:,1:end-1)'; 
            y = yMat(:,1:end-1)'; 
            
            
            for k=1:cMapLength
            plot(x(:,idx(:,1)==k),y(:,idx(:,1)==k),'color',cMap(k,:),'lineWidth',1);
            end
        else 
            hold on 
            x = projData.xCoordOutAllTracks; 
            y = projData.yCoordOutAllTracks; 
            
            plot(x',y','w'); 
            
            
            
            
            
        end 
    for i = 1:numel(allBCoordsYX) 
        roiYX = allBCoordsYX{i};
    plot(roiYX(:,2),roiYX(:,1),'color','w','linewidth', 1);
    end 
    title([type{j} ' of growth calculated specifically in region: Max ' type{j}...
            ' ' num2str(dataRange) ' ' typeUnits{j}]);
        
        saveas(maskFig,[saveDir filesep 'overlay_growth_tracks_' type{j} 'calculated_specifically_in_ROI.fig'])     
          saveas(maskFig,[saveDir filesep 'overlay_growth_tracks_' type{j} 'calculated_specifically_in_ROI.eps'],'psc2'); 
    end 
    
    
close(maskFig); 
end












%% STACKED HISTOGRAMS
[xMat,yMat]=plusTipGetSubtrackCoords(projData,[],1);

[dummy,speedLifeDispMat]=plusTipDynamParam(allData,[],1,0);

for j=1:3
    data=speedLifeDispMat(:,typeMatIndx{j});
    % create x-axis bins spanning all values
    n=linspace(0,typeLim{j},25);

    % bin the samples
    x1 = histc(data(:,1),n)/numel(data); % growth
    x2 = histc(data(:,2),n)/numel(data); % fgap
    x3 = histc(data(:,3),n)/numel(data); % bgap

    % put the binned values into a matrix for the stacked plot
    M=nan(max([length(x1) length(x2) length(x3)]),3);
    M(1:length(x1),1)=x1;
    M(1:length(x2),2)=x2;
    M(1:length(x3),3)=x3;

    % make the plot
    figure
    bar(n,M,'stack')
    colormap([1 0 0; 0 0 1; 0 1 0])
    legend('growth','fgap','bgap','Location','best')
    xlabel([type{j} ' (' typeUnits{j} ')']);
    ylabel('Frequency');
    
    saveas(gcf,[saveDir filesep 'histogram_' type{j} '.fig'])
    saveas(gcf,[saveDir filesep 'histogram_' type{j} fileExt])
    close(gcf)
    if ishandle(wtBar)
        waitbar(.5+j/12,wtBar,'Creating stacked histograms');
    end
end


        
    



%% RESCUE FROM PAUSE AND SHRINKAGE PLOTS
overlayFig = figure; 
imshow(img,[]); 
set(gca,'Position',[0 0 1 .95]);
hold on

% pause info
nFgaps=length(fIdx);
xCatPause=zeros(nFgaps,1); yCatPause=zeros(nFgaps,1);
xResPause=zeros(nFgaps,1); yResPause=zeros(nFgaps,1);
for iFgap=1:nFgaps
    c=fIdx(iFgap);
    % coordinates for pause initiation (growth to pause)
    idx=find(~isnan(xMat(c,:)),1,'first');
    xCatPause(iFgap)=xMat(c,idx);
    yCatPause(iFgap)=yMat(c,idx);
    % coordinates for pause termination (pause to growth)    
    idx=find(~isnan(xMat(c,:)),1,'last');
    xResPause(iFgap)=xMat(c,idx);
    yResPause(iFgap)=yMat(c,idx);
end

% shrink info
nBgaps=length(bIdx);
xCatShrink=zeros(nBgaps,1); yCatShrink=zeros(nBgaps,1);
xResShrink=zeros(nBgaps,1); yResShrink=zeros(nBgaps,1);
for iBgap=1:nBgaps
    c=bIdx(iBgap);
    % coordinates for shrinkage initiation (growth to shrinkage)
    idx=find(~isnan(xMat(c,:)),1,'first');
    xCatShrink(iBgap)=xMat(c,idx);
    yCatShrink(iBgap)=yMat(c,idx);
    % coordinates for shrinkage termination (shrinkage to growth)    
    idx=find(~isnan(xMat(c,:)),1,'last');
    xResShrink(iBgap)=xMat(c,idx);
    yResShrink(iBgap)=yMat(c,idx);
end

if ishandle(wtBar),  waitbar(.75,wtBar,'Creating sites plots'); end

% make the pause/shrink initiation plots
figure(overlayFig);
  delete(findobj(overlayFig,'Type','line'));
h=scatter(xCatPause,yCatPause,'y','filled');
title('pause initiation sites (growth to pause)')
saveas(gcf,[saveDir filesep 'pause initiation sites (growth to pause)' '.fig'])
saveas(gcf,[saveDir filesep 'pause initiation sites (growth to pause)' fileExt])

delete(h);
h=scatter(xCatShrink,yCatShrink,'r','filled');
title('shrinkage initiation sites (growth to shrinkage)')
saveas(gcf,[saveDir filesep 'shrinkage initiation sites (growth to shrinkage)' '.fig'])
saveas(gcf,[saveDir filesep 'shrinkage initiation sites (growth to shrinkage)' fileExt])

h(2)=scatter(xCatPause,yCatPause,'y','filled');
title('pause (yellow) and shrinkage (red) initiation sites')
saveas(gcf,[saveDir filesep 'pause and shrinkage initiation sites' '.fig'])
saveas(gcf,[saveDir filesep 'pause and shrinkage initiation sites' fileExt])

if ishandle(wtBar),  waitbar(.875,wtBar,'Creating sites plots'); end

% make the pause/shrinkage termination plots
delete(h);
h=scatter(xResPause,yResPause,'y','filled');
title('pause termination sites (pause to growth)')
saveas(gcf,[saveDir filesep 'pause termination sites (pause to growth)' '.fig'])
saveas(gcf,[saveDir filesep 'pause termination sites (pause to growth)' fileExt])

delete(h);
h=scatter(xResShrink,yResShrink,'r','filled');
title('shrinkage termination sites (shrinkage to growth)')
saveas(gcf,[saveDir filesep 'shrinkage termination sites(shrinkage to growth)' '.fig'])
saveas(gcf,[saveDir filesep 'shrinkage termination sites(shrinkage to growth)' fileExt])

h(2) = scatter(xResPause,yResPause,'y','filled');
title('pause (yellow) and shrinkage (red) termination sites')    
saveas(gcf,[saveDir filesep 'pause and shrinkage termination sites' '.fig'])
saveas(gcf,[saveDir filesep 'pause and shrinkage termination sites' fileExt])

delete(h);
p = [1 0 1]; 
% h = scatter(xResShrink,yResShrink,'b','filled'); 
% h(2) = scatter(xCatPause,yCatPause,'c','filled'); 
% h(3) = scatter(xCatShrink,yCatShrink,'y','filled'); 
%  plotScaleBar(100,2,'Label','10 um','Location','NorthEast','FontSize',14),
% saveas(gcf,[saveDir filesep 'for_paper.eps'],'psc2'); 

close(gcf)

function [img2show]=addMaskInColor(img,roiMask,c)
%subfunction to add new polygon outline to composite image - this is needed
%because you can't pass vector graphics info to roipoly function, and we
%want to be able to visualize the regions that have already been selected.
temp=double(bwmorph(roiMask,'remove'));
borderIdx=find(temp);
nPix=numel(roiMask);

img2show=img;
img2show(borderIdx)=c(1);
img2show(borderIdx+nPix)=c(2);
img2show(borderIdx+2*nPix)=c(3);

end

close all
if ishandle(wtBar), close(wtBar); end
end
