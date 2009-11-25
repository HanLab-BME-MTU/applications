function plusTipPlotResults(projData,remBegEnd,timeRange,speedLim,lifeLim,dispLim,saveDir)
% make plots of speed/lifetime/displacement, histograms, initition/termination pts for fgaps and bgaps

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
numFrames=projData.numFrames;
% check whether a time range for plotting was input
if nargin<3 || isempty(timeRange)
    timeRange=[1 numFrames];
else
    if timeRange(1)<1
        timeRange(1)=1;
    end
    if timeRange(2)>numFrames
        timeRange(2)=numFrames;
    end
end

% format image/analysis directory paths
imDir=formatPath(projData.imDir);
anDir=formatPath(projData.anDir);

if nargin<7 || isempty(saveDir)
    saveDir=uigetdir(anDir,'Please choose output directory');
end


% get first image from imDir
[listOfImages] = searchFiles('.tif',[],imDir,0);
img = double(imread([listOfImages{1,2} filesep listOfImages{1,1}]));

% extract track info (block matrix)
dataMatMerge=plusTipMergeSubtracks(projData); % merged data is first output
allData=abs(dataMatMerge);
allData(:,6)=allData(:,6).* projData.secPerFrame; % convert lifetimes to seconds
allData(:,7)=allData(:,7).*(projData.pixSizeNm/1000); % convert displacements to microns


if nargin<4 || isempty(speedLim)
    speedLim=prctile(allData(:,4),95);
end

if nargin<5 || isempty(lifeLim)
    lifeLim=prctile(allData(:,6),95);
end

if nargin<6 || isempty(dispLim)
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

%% GROWTH SUBTRACKS

x=xMat(gIdx,1:end-1)';
y=yMat(gIdx,1:end-1)';

% SPEED
speed=allData(gIdx,4);
m=ceil(speedLim);
speed(speed>m)=m;
mapper=linspace(0,m,cMapLength)';

% get closest colormap index for each feature
D=createDistanceMatrix(speed,mapper);
[sD,idx]=sort(abs(D),2);

% make the plot
figure('Position',figPos)
imagesc(img); colormap gray;
hold on
for i=1:size(x,2)
    plot(x(:,i),y(:,i),'color',cMap(idx((i),1),:));
end
xlabel(['speed map of growth tracks, max speed ' num2str(m) ' microns/min'])



% LIFETIME
life=allData(gIdx,6);
m=ceil(lifeLim);
life(life>m)=m;
mapper=linspace(0,m,cMapLength)';

% get closest colormap index for each feature
D=createDistanceMatrix(life,mapper);
[sD,idx]=sort(abs(D),2);

% make the plot
figure('Position',figPos)
imagesc(img); colormap gray;
hold on
for i=1:size(x,2)
    plot(x(:,i),y(:,i),'color',cMap(idx((i),1),:));
end
xlabel(['lifetime map of growth tracks, max lifetime ' num2str(m) ' sec'])



% DISPLACEMENT
disp=allData(gIdx,7);
m=ceil(dispLim);
disp(disp>m)=m;
mapper=linspace(0,m,cMapLength)';

% get closest colormap index for each feature
D=createDistanceMatrix(disp,mapper);
[sD,idx]=sort(abs(D),2);

% make the plot
figure('Position',figPos)
imagesc(img); colormap gray;
hold on
for i=1:size(x,2)
    plot(x(:,i),y(:,i),'color',cMap(idx((i),1),:));
end
xlabel(['displacement map of growth tracks, max displacement ' num2str(m) ' microns'])

%% FGAP SUBTRACKS

x=xMat(fIdx,1:end-1)';
y=yMat(fIdx,1:end-1)';

% SPEED
speed=allData(fIdx,4);
m=ceil(speedLim);
speed(speed>m)=m;
mapper=linspace(0,m,cMapLength)';

% get closest colormap index for each feature
D=createDistanceMatrix(speed,mapper);
[sD,idx]=sort(abs(D),2);

% make the plot
figure('Position',figPos)
imagesc(img); colormap gray;
hold on
for i=1:size(x,2)
    plot(x(:,i),y(:,i),'color',cMap(idx((i),1),:));
end
xlabel(['speed map of fgaps, max speed ' num2str(m) ' microns/min'])



% LIFETIME
life=allData(fIdx,6);
m=ceil(lifeLim);
life(life>m)=m;
mapper=linspace(0,m,cMapLength)';

% get closest colormap index for each feature
D=createDistanceMatrix(life,mapper);
[sD,idx]=sort(abs(D),2);

% make the plot
figure('Position',figPos)
imagesc(img); colormap gray;
hold on
for i=1:size(x,2)
    plot(x(:,i),y(:,i),'color',cMap(idx((i),1),:));
end
xlabel(['lifetime map of fgaps, max lifetime ' num2str(m) ' sec'])



% DISPLACEMENT
disp=allData(fIdx,7);
m=ceil(dispLim);
disp(disp>m)=m;
mapper=linspace(0,m,cMapLength)';

% get closest colormap index for each feature
D=createDistanceMatrix(disp,mapper);
[sD,idx]=sort(abs(D),2);

% make the plot
figure('Position',figPos)
imagesc(img); colormap gray;
hold on
for i=1:size(x,2)
    plot(x(:,i),y(:,i),'color',cMap(idx((i),1),:));
end
xlabel(['displacement map of fgaps, max displacement ' num2str(m) ' microns'])


%% BGAP SUBTRACKS

x=xMat(bIdx,1:end-1)';
y=yMat(bIdx,1:end-1)';

% SPEED
speed=allData(bIdx,4);
m=ceil(speedLim);
speed(speed>m)=m;
mapper=linspace(0,m,cMapLength)';

% get closest colormap index for each feature
D=createDistanceMatrix(speed,mapper);
[sD,idx]=sort(abs(D),2);

% make the plot
figure('Position',figPos)
imagesc(img); colormap gray;
hold on
for i=1:size(x,2)
    plot(x(:,i),y(:,i),'color',cMap(idx((i),1),:));
end
xlabel(['speed map of bgaps, max speed ' num2str(m) ' microns/min'])



% LIFETIME
life=allData(bIdx,6);
m=ceil(lifeLim);
life(life>m)=m;
mapper=linspace(0,m,cMapLength)';

% get closest colormap index for each feature
D=createDistanceMatrix(life,mapper);
[sD,idx]=sort(abs(D),2);

% make the plot
figure('Position',figPos)
imagesc(img); colormap gray;
hold on
for i=1:size(x,2)
    plot(x(:,i),y(:,i),'color',cMap(idx((i),1),:));
end
xlabel(['lifetime map of bgaps, max lifetime ' num2str(m) ' sec'])



% DISPLACEMENT
disp=allData(bIdx,7);
m=ceil(dispLim);
disp(disp>m)=m;
mapper=linspace(0,m,cMapLength)';

% get closest colormap index for each feature
D=createDistanceMatrix(disp,mapper);
[sD,idx]=sort(abs(D),2);

% make the plot
figure('Position',figPos)
imagesc(img); colormap gray;
hold on
for i=1:size(x,2)
    plot(x(:,i),y(:,i),'color',cMap(idx((i),1),:));
end
xlabel(['displacement map of bgaps, max displacement ' num2str(m) ' microns'])

%% STACKED HISTOGRAMS

[dummy,speedLifeDispMat]=plusTipDynamParam(allData);

for iParam=1:3
    switch iParam
        case 1
            data=speedLifeDispMat(:,1:3); % speed
            lim=speedLim;
            titleStr='speed';
            xStr='speed (microns/min)';
        case 2
            data=speedLifeDispMat(:,4:6); % lifetime
            lim=lifeLim;
            titleStr='lifetime';
            xStr='lifetime (sec)';
        case 3
            data=speedLifeDispMat(:,7:9); % displacement
            lim=dispLim;
            titleStr='displacement';
            xStr='displacement (microns)';
    end


    % create x-axis bins spanning all values
    n=linspace(0,lim,25);

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
    xlabel(xStr);
    ylabel('frequency');
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


% make the pause/shrink initiation plots
figure('Position',figPos)
imagesc(img); colormap gray;
hold on
scatter(xCatPause,yCatPause,'y','filled')
xlabel('pause initiation sites (growth to pause)')

figure('Position',figPos)
imagesc(img); colormap gray;
hold on
scatter(xCatShrink,yCatShrink,'r','filled')
xlabel('shrinkage initiation sites (growth to shrinkage)')

figure('Position',figPos)
imagesc(img); colormap gray;
hold on
scatter(xCatPause,yCatPause,'y','filled')
scatter(xCatShrink,yCatShrink,'r','filled')
xlabel('pause (yellow) and shrinkage (red) initiation sites')

% make the pause/shrinkage termination plots
figure('Position',figPos)
imagesc(img); colormap gray;
hold on
scatter(xResPause,yResPause,'y','filled')
xlabel('pause termination sites (pause to growth)')

figure('Position',figPos)
imagesc(img); colormap gray;
hold on
scatter(xResShrink,yResShrink,'r','filled')
xlabel('shrinkage termination sites (shrinkage to growth)')

figure('Position',figPos)
imagesc(img); colormap gray;
hold on
scatter(xResPause,yResPause,'y','filled')
scatter(xResShrink,yResShrink,'r','filled')
xlabel('pause (yellow) and shrinkage (red) termination sites')    

%% save all the plots
figNum=1;
for iType=1:3
    switch iType
        case 1
            typeStr='growth';
        case 2
            typeStr='fgap';
        case 3
            typeStr='bgap';
    end

    % save the track overlays
    for iParam=1:3
        switch iParam
            case 1
                paramStr='speed';
            case 2
                paramStr='lifetime';
            case 3
                paramStr='displacement';
        end
        saveas(figNum,[saveDir filesep 'overlay_' typeStr '_' paramStr '.fig'])
        saveas(figNum,[saveDir filesep 'overlay_' typeStr '_' paramStr '.emf'])
        close(figNum)
        figNum=figNum+1;
    end
end
% save the stacked histograms
for iParam=1:3
    switch iParam
        case 1
            paramStr='speed';
        case 2
            paramStr='lifetime';
        case 3
            paramStr='displacement';
    end
    saveas(figNum,[saveDir filesep 'histogram_' paramStr '.fig'])
    saveas(figNum,[saveDir filesep 'histogram_' paramStr '.emf'])
    close(figNum)
    figNum=figNum+1;
end

% save the catastophe/rescue plots
saveas(figNum,[saveDir filesep 'pause initiation sites (growth to pause)' '.fig'])
saveas(figNum,[saveDir filesep 'pause initiation sites (growth to pause)' '.emf'])
close(figNum)
figNum=figNum+1;

saveas(figNum,[saveDir filesep 'shrinkage initiation sites (growth to shrinkage)' '.fig'])
saveas(figNum,[saveDir filesep 'shrinkage initiation sites (growth to shrinkage)' '.emf'])
close(figNum)
figNum=figNum+1;

saveas(figNum,[saveDir filesep 'pause and shrinkage initiation sites' '.fig'])
saveas(figNum,[saveDir filesep 'pause and shrinkage initiation sites' '.emf'])
close(figNum)
figNum=figNum+1;

saveas(figNum,[saveDir filesep 'pause termination sites (pause to growth)' '.fig'])
saveas(figNum,[saveDir filesep 'pause termination sites (pause to growth)' '.emf'])
close(figNum)
figNum=figNum+1;

saveas(figNum,[saveDir filesep 'shrinkage termination sites(shrinkage to growth)' '.fig'])
saveas(figNum,[saveDir filesep 'shrinkage termination sites(shrinkage to growth)' '.emf'])
close(figNum)
figNum=figNum+1;

saveas(figNum,[saveDir filesep 'pause and shrinkage termination sites' '.fig'])
saveas(figNum,[saveDir filesep 'pause and shrinkage termination sites' '.emf'])
close(figNum)
figNum=figNum+1;



