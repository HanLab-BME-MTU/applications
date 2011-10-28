function plusTipPlotResults(projData,remBegEnd,timeRange,speedLim,lifeLim,dispLim,saveDir)
% plusTipPlotResults make 18 summary overlay plots and histograms
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

% get first image from imDir
[listOfImages] = searchFiles('.tif',[],imDir,0);
[filename pathname]  = uigetfile('*','Select an Image File',imDir); 
%img = double(imread([listOfImages{1,2} filesep listOfImages{1,1}]));
img = double(imread([pathname filesep filename])); 

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

% get track type and matrices containing xy-coordinates for all subtracks
trackType=allData(:,5);
trackStarts=allData(:,2);
trackEnds=allData(:,3);
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


cMapLength=128; cMap=jet(cMapLength);

% Create MT dynamics maps of speed, lifetime and displacement for growth,
% fgaps and bgaps events
type= {'speed','lifetime','displacement'};
typeIndx= {4,6,7};
typeMatIndx = {1:3,4:6,7:9};
typeLim = {speedLim;lifeLim;dispLim};
typeUnits = {speedLim;lifeLim;dispLim};
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
        figure('Position',figPos)
        imagesc(img); colormap gray;
        axis off
        hold on
        for k=1:cMapLength
            plot(x(:,idx(:,1)==k),y(:,idx(:,1)==k),'color',cMap(k,:));
        end
        xlabel([type{i} 'map of ' event{j} ', max ' type{i}...
            ' ' num2str(dataRange) ' ' typeUnits{i}]);
        
        saveas(gcf,[saveDir filesep 'overlay_' event{j} '_' type{i} '.fig'])
        saveas(gcf,[saveDir filesep 'overlay_' event{j} '_' type{i} fileExt])
        close(gcf)
        
        if ishandle(wtBar)
            waitbar(((i-1)*3+ j)/18,wtBar,'Creating overlay maps');
        end
        
    end
end

%% STACKED HISTOGRAMS

[dummy,speedLifeDispMat]=plusTipDynamParam(allData,[],1,0);

for j=1:3
    data=speedLifeDispMat(:,typeMatIndx{j});
    % create x-axis bins spanning all values
    n=linspace(0,typeLim{j},25);

    % bin the samples
    [x1,dummy] = histc(data(:,1),n); % growth
    [x2,dummy] = histc(data(:,2),n); % fgap
    [x3,dummy] = histc(data(:,3),n); % bgap

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
figure('Position',figPos)
imagesc(img); colormap gray;
hold on
scatter(xCatPause,yCatPause,'y','filled')
xlabel('pause initiation sites (growth to pause)')
saveas(gcf,[saveDir filesep 'pause initiation sites (growth to pause)' '.fig'])
saveas(gcf,[saveDir filesep 'pause initiation sites (growth to pause)' fileExt])
close(gcf)

figure('Position',figPos)
imagesc(img); colormap gray;
hold on
scatter(xCatShrink,yCatShrink,'r','filled')
xlabel('shrinkage initiation sites (growth to shrinkage)')

saveas(gcf,[saveDir filesep 'shrinkage initiation sites (growth to shrinkage)' '.fig'])
saveas(gcf,[saveDir filesep 'shrinkage initiation sites (growth to shrinkage)' fileExt])
close(gcf)

figure('Position',figPos)
imagesc(img); colormap gray;
hold on
scatter(xCatPause,yCatPause,'y','filled')
scatter(xCatShrink,yCatShrink,'r','filled')
xlabel('pause (yellow) and shrinkage (red) initiation sites')

saveas(gcf,[saveDir filesep 'pause and shrinkage initiation sites' '.fig'])
saveas(gcf,[saveDir filesep 'pause and shrinkage initiation sites' fileExt])
close(gcf)

if ishandle(wtBar),  waitbar(.875,wtBar,'Creating sites plots'); end

% make the pause/shrinkage termination plots
figure('Position',figPos)
imagesc(img); colormap gray;
hold on
scatter(xResPause,yResPause,'y','filled')
xlabel('pause termination sites (pause to growth)')

saveas(gcf,[saveDir filesep 'pause termination sites (pause to growth)' '.fig'])
saveas(gcf,[saveDir filesep 'pause termination sites (pause to growth)' fileExt])
close(gcf)


figure('Position',figPos)
imagesc(img); colormap gray;
hold on
scatter(xResShrink,yResShrink,'r','filled')
xlabel('shrinkage termination sites (shrinkage to growth)')

saveas(gcf,[saveDir filesep 'shrinkage termination sites(shrinkage to growth)' '.fig'])
saveas(gcf,[saveDir filesep 'shrinkage termination sites(shrinkage to growth)' fileExt])
close(gcf)


figure('Position',figPos)
imagesc(img); colormap gray;
hold on
scatter(xResPause,yResPause,'y','filled')
scatter(xResShrink,yResShrink,'r','filled')
xlabel('pause (yellow) and shrinkage (red) termination sites')    

saveas(gcf,[saveDir filesep 'pause and shrinkage termination sites' '.fig'])
saveas(gcf,[saveDir filesep 'pause and shrinkage termination sites' fileExt])
close(gcf)


if ishandle(wtBar), close(wtBar); end
