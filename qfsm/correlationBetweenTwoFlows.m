function [] = correlationBetweenTwoFlows(pathMD)
% correlationBetweenTwoFlows(pathMD) takes two flows from the movieData,
% get one or more ROIs from user input, get average flow over space, do
% time-crosscorrelation with flow speeds, and get DCS and VMCS with
% averaging measured velocities over some time period. 
% input:    pathMD:    path to the movieData file
% output:     map of ROIs on each channel
%                   time series of flow mag of each channel
%                   DCS and VMCS of the two flows per some segments of
%                   times (e.g. per 10 frames)
% backing up the original one
% flow_mag is in unit of nm/min.

% Load movieData
movieDataPath = [pathMD '/movieData.mat'];
movieData = MovieData.load(movieDataPath);
nFrames = movieData.nFrames_;
nChannels = length(movieData.channels_);

% Load flow vectors
iFlow = movieData.getProcessIndex('FlowAnalysisProcess');
if isempty(iFlow)
    iFlow = movieData.getProcessIndex('FlowTrackingProcess');
    if isempty(iFlow)
        error('Flow tracking has to be run.')
    else
        display('Flow analysis process has not been run, flow tracking process is used instead.')
    end
end
flowProcess = movieData.getProcess(iFlow);
flow1 = flowProcess.loadChannelOutput(1,'output','Md');
% Load segmented masks
dt = movieData.timeInterval_; 
res = movieData.pixelSize_;

%Create string for current directory
[pathstr,~,~] = fileparts(flowProcess.outFilePaths_{1});
outputPath = [pathstr filesep 'Correlation_Analysis'];
imgPath = [outputPath filesep 'imgs'];
dataPath = [outputPath filesep 'data'];
tifPath = [imgPath filesep 'tifs'];
figPath = [imgPath filesep 'figs'];
if ~exist(figPath,'dir')
    mkdir(imgPath);
    mkdir(dataPath);
    mkdir(tifPath);
    mkdir(figPath);
end

flowMagFrame = cell(nChannels,nFrames);
iChan=1;
flow1 = flowProcess.loadChannelOutput(iChan,'output','Md');
flow2 = flowProcess.loadChannelOutput(2,'output','Md');
% Flow interpolation for all points for all frames to prevent lact of data
% for correlation analysis
posAll = [];
% Getting all the positions
for jj=1:nFrames
    curflow1 = flow1{jj};
    curflow2 = flow2{jj};
    posAll = [posAll; curflow1(:,2:-1:1);curflow2(:,2:-1:1)];
end
% Ball Quenry to get reasonably separated points (5 pix)
max_distance = 5;
idx = KDTreeBallQuery(posAll, posAll, max_distance);
valid = true(numel(idx),1);
for i = 1 : numel(idx)
    if ~valid(i), continue; end
    neighbors = idx{i}(idx{i} ~= i);
    valid(neighbors) = false;
end
posFiltered = posAll(valid, :);

%% Interpolation
for jj=1:nFrames
    curflow1 = flow1{jj};
    curflow2 = flow2{jj};

    flow1interp{jj} = vectorFieldSparseInterp(curflow1,posFiltered(:,2:-1:1), 50, 20,[]);
    flow2interp{jj} = vectorFieldSparseInterp(curflow2,posFiltered(:,2:-1:1), 50, 20,[]);    
end

%% Sampling
for jj=1:nFrames
    if iChan == 1 && jj==1
        % Load segmented masks
        I=double(movieData.channels_(iChan).loadImage(jj));
        dI = double(I)/max(I(:));
        I2=double(movieData.channels_(iChan).loadImage(nFrames));
        dI2 = double(I2)/max(I2(:));
        combI(:,:,1) = dI;
        combI(:,:,2) = dI2;
        combI(:,:,3) = 0;
        h=figure;
        imshow(combI,[]); hold on
        mag = 100/dt;
        flow1f = flow1{jj};
        flow1l = flow1{end};

        quiver(flow1f(:,2),flow1f(:,1),mag*(flow1f(:,4)-flow1f(:,2)),mag*(flow1f(:,3)-flow1f(:,1)),0,'k','Linewidth',0.5);
        quiver(flow1l(:,2),flow1l(:,1),mag*(flow1l(:,4)-flow1l(:,2)),mag*(flow1l(:,3)-flow1l(:,1)),0,'b','Linewidth',0.5);
        flow2f = flow2{jj};
        flow2l = flow2{end};
        quiver(flow2f(:,2),flow2f(:,1),mag*(flow2f(:,4)-flow2f(:,2)),mag*(flow2f(:,3)-flow2f(:,1)),0,'r','Linewidth',0.5);
        quiver(flow2l(:,2),flow2l(:,1),mag*(flow2l(:,4)-flow2l(:,2)),mag*(flow2l(:,3)-flow2l(:,1)),0,'m','Linewidth',0.5);
        % Get user input by impoly that will be excluded
        nPoly = input('How many ROIs do you want to make for flow correlation analysis? ');
        hpoly = cell(nPoly,1);
        polyMask = cell(nPoly,1);
        display(['Please draw ' num2str(nPoly) ' polygon(s) where you analyze the vector field. Finish with right click after a final point for each polygon'])
        for k=1:nPoly
            hpoly{k} = impoly;
            polyMask{k} = createMask(hpoly{k});
            display(['Drawing done for polygon # ' num2str(k) ' out of ' num2str(nPoly) '.'])
        end
        hgexport(h,strcat(tifPath,'/ROIs.tif'),hgexport('factorystyle'),'Format','tiff')
        hgsave(h,strcat(figPath,'/ROIs.fig'),'-v7.3')
        close(h)
    end
    % Now let's do the analysis!
    curFlow1 = flow1interp{jj};
    curFlow2 = flow2interp{jj};
    for k=1:nPoly
        % average flow in each ROI over the space
        % store flow in each ROI
        inMaskIdx1 = arrayfun(@(i,j) maskVectors(i,j,polyMask{k}),curFlow1(:,2),curFlow1(:,1));
        inMaskIdx2 = arrayfun(@(i,j) maskVectors(i,j,polyMask{k}),curFlow2(:,2),curFlow2(:,1));
        roiCurFlow1 = curFlow1(inMaskIdx1,:);
        roiCurFlow2 = curFlow2(inMaskIdx2,:);
        % Averaging flows into one flow vector
        avgRoiFlow1_vec = [nanmean(roiCurFlow1(:,4)-roiCurFlow1(:,2)), nanmean(roiCurFlow1(:,3)-roiCurFlow1(:,1))];
        avgRoiFlow2_vec = [nanmean(roiCurFlow2(:,4)-roiCurFlow2(:,2)), nanmean(roiCurFlow2(:,3)-roiCurFlow2(:,1))];
        avgRoiFlow1_vec_std = [nanstd(roiCurFlow1(:,4)-roiCurFlow1(:,2)), nanstd(roiCurFlow1(:,3)-roiCurFlow1(:,1))];
        avgRoiFlow2_vec_std = [nanstd(roiCurFlow2(:,4)-roiCurFlow2(:,2)), nanstd(roiCurFlow2(:,3)-roiCurFlow2(:,1))];
        if jj==1
            avgRoiFlow1_pos{k} = [nanmean(roiCurFlow1(:,2)), nanmean(roiCurFlow1(:,1))];
            avgRoiFlow2_pos{k} = [nanmean(roiCurFlow2(:,2)), nanmean(roiCurFlow2(:,1))];
        end
        roiFlow1{k,jj} = roiCurFlow1;
        roiFlow2{k,jj} = roiCurFlow2;

        avgRoiFlow1{k}.pos(jj,:) = avgRoiFlow1_pos{k};
        avgRoiFlow1{k}.vec(jj,:) = avgRoiFlow1_vec;
        avgRoiFlow1{k}.vec_std(jj,:) = avgRoiFlow1_vec_std;
        avgRoiFlow1{k}.speed(jj) = (avgRoiFlow1_vec(1).^2+avgRoiFlow1_vec(2).^2).^0.5*res*(60/dt); % unit is nm/min
        avgRoiFlow1{k}.speed_std(jj) = (avgRoiFlow1_vec_std(1).^2+avgRoiFlow1_vec_std(2).^2).^0.5*res*(60/dt); % unit is nm/min

        avgRoiFlow2{k}.pos(jj,:) = avgRoiFlow2_pos{k};
        avgRoiFlow2{k}.vec(jj,:) = avgRoiFlow2_vec;
        avgRoiFlow2{k}.vec_std(jj,:) = avgRoiFlow2_vec_std;
        avgRoiFlow2{k}.speed(jj) = (avgRoiFlow2_vec(1).^2+avgRoiFlow2_vec(2).^2).^0.5*res*(60/dt); % unit is nm/min
        avgRoiFlow2{k}.speed_std(jj) = (avgRoiFlow2_vec_std(1).^2+avgRoiFlow2_vec_std(2).^2).^0.5*res*(60/dt); % unit is nm/min

        % save flow magnitudes:
        avgRoiFlow1{k}.t(jj) = (jj-1)*movieData.timeInterval_/60; % in min
        avgRoiFlow2{k}.t(jj) = (jj-1)*movieData.timeInterval_/60; % in min
    end
end

%% Time-crosscorrelation with flow speeds
h=figure; hold on
for k=1:nPoly
    [c(k,:),lags] = xcorr(avgRoiFlow1{k}.speed(1:end-3),avgRoiFlow2{k}.speed(1:end-3),'coeff');
    plot(lags*dt,c(k,:), 'Color',[0.5 0.5 0.5],'Linewidth',1)
end
plot(lags*dt,nanmean(c), 'k','Linewidth',2)
xlabel('Time (s)')
ylabel('Normalized correlation score')
hgexport(h,strcat(tifPath,'/cross_correlation.tif'),hgexport('factorystyle'),'Format','tiff')
hgsave(h,strcat(figPath,'/cross_correlation.fig'),'-v7.3')
save([dataPath '/cross_correlation.mat'],'lags','c','dt')
close(h)

%%  DCS and VMCS with averaging measured velocities over some time period.
% Time-average the flow vectors
nTavgFrames = floor(nFrames/10);
TavgFlow1 = cell(nTavgFrames,1);
TavgFlow2 = cell(nTavgFrames,1);

% Sampling again
for jj=1:nFrames
    if iChan == 1 && jj==1
        % Load segmented masks
        I=double(movieData.channels_(iChan).loadImage(jj));
        dI = double(I)/max(I(:));
        I2=double(movieData.channels_(iChan).loadImage(nFrames));
        dI2 = double(I2)/max(I2(:));
        combI(:,:,1) = dI;
        combI(:,:,2) = dI2;
        combI(:,:,3) = 0;
        h=figure;
        imshow(combI,[]); hold on
        mag = 100/dt;
        flow1f = flow1interp{jj};
        flow1l = flow1interp{end};
        quiver(flow1f(:,2),flow1f(:,1),mag*(flow1f(:,4)-flow1f(:,2)),mag*(flow1f(:,3)-flow1f(:,1)),0,'k','Linewidth',0.5);
        quiver(flow1l(:,2),flow1l(:,1),mag*(flow1l(:,4)-flow1l(:,2)),mag*(flow1l(:,3)-flow1l(:,1)),0,'b','Linewidth',0.5);
        flow2f = flow2interp{jj};
        flow2l = flow2interp{end};
        quiver(flow2f(:,2),flow2f(:,1),mag*(flow2f(:,4)-flow2f(:,2)),mag*(flow2f(:,3)-flow2f(:,1)),0,'r','Linewidth',0.5);
        quiver(flow2l(:,2),flow2l(:,1),mag*(flow2l(:,4)-flow2l(:,2)),mag*(flow2l(:,3)-flow2l(:,1)),0,'m','Linewidth',0.5);
        % Get user input by impoly that will be excluded
        display(['Please draw  a convex polygon where you analyze DCS and VMCS. Finish with right click after a final point for the polygon'])
        hpoly2 = impoly;
        polyOne = createMask(hpoly2);
        hgexport(h,strcat(tifPath,'/oneROI.tif'),hgexport('factorystyle'),'Format','tiff')
        hgsave(h,strcat(figPath,'/oneROI.fig'),'-v7.3')

        close(h)
    end
    % Now let's do the analysis!
    curFlow1 = flow1interp{jj};
    curFlow2 = flow2interp{jj};
    % average flow in each ROI over the space
    % store flow in each ROI
    inMaskIdx1 = arrayfun(@(i,j) maskVectors(i,j,polyOne),curFlow1(:,2),curFlow1(:,1));
    inMaskIdx2 = arrayfun(@(i,j) maskVectors(i,j,polyOne),curFlow2(:,2),curFlow2(:,1));
    flow1interpFiltered{jj} = curFlow1(inMaskIdx1,:);
    flow2interpFiltered{jj} = curFlow2(inMaskIdx2,:);
end
%%  DCS and VMCS calculation
for jj=1:nTavgFrames
    % every 10 frames, average per each position
    for pp=1:10
        Flow1array(:,:,pp) = flow1interpFiltered{(jj-1)*10+pp};
        Flow2array(:,:,pp) = flow2interpFiltered{(jj-1)*10+pp};
    end
    TavgFlow1{jj} = mean(Flow1array,3);
    TavgFlow2{jj} = mean(Flow2array,3);
    % DCS
    curFlow1 = TavgFlow1{jj};
    curFlow2 = TavgFlow2{jj};
    vec1 = [curFlow1(:,4) - curFlow1(:,2), curFlow1(:,3) - curFlow1(:,1)];
    vec2 = [curFlow2(:,4) - curFlow2(:,2), curFlow2(:,3) - curFlow2(:,1)];
    curDCS = dot(vec1,vec2,2)./(sum(vec1.^2,2).^0.5.*sum(vec2.^2,2).^0.5);
    curVMCS = sum(vec2.^2,2).^0.5./sum(vec1.^2,2).^0.5;
    DCS{jj} = [curFlow1(:,2),curFlow1(:,1),curDCS];
    %VMCS
    VMCS{jj} = [curFlow1(:,2),curFlow1(:,1),curVMCS];
end
%% Showing each vectors, DCS and VMCS
% color map
numColors = 128;

%Blue-Green to yellow.
cMap = [linspace(0,1,numColors/2).' ...
   linspace(0.5,1,numColors/2).' ...
   linspace(0.5,0.05,numColors/2).'];
%Yellow to red.
cMap = [cMap; [ones(numColors/2,1) ...
   linspace(1,0,numColors/2).' zeros(numColors/2,1)]];

tifPath_Flow1 = [tifPath '/flow1'];
figPath_Flow1 = [figPath '/flow1'];
tifPath_Flow2 = [tifPath '/flow2'];
figPath_Flow2 = [figPath '/flow2'];
tifPath_FlowCombined = [tifPath '/flowCombined'];
figPath_FlowCombined = [figPath '/flowCombined'];
tifPath_DCS = [tifPath '/DCS'];
figPath_DCS = [figPath '/DCS'];
tifPath_VMCS = [tifPath '/VMCS'];
figPath_VMCS = [figPath '/VMCS'];
if ~exist(figPath_VMCS,'dir')
    mkdir(tifPath_Flow1);
    mkdir(figPath_Flow1);
    mkdir(tifPath_Flow2);
    mkdir(figPath_Flow2);
    mkdir(tifPath_FlowCombined);
    mkdir(figPath_FlowCombined);
    mkdir(tifPath_DCS);
    mkdir(figPath_DCS);
    mkdir(tifPath_VMCS);
    mkdir(figPath_VMCS);
end
iiformat = ['%.' '3' 'd'];
roiDCS = zeros(nPoly,nTavgFrames);
roiVMCS = zeros(nPoly,nTavgFrames);
for jj=1:nTavgFrames
    I1=double(movieData.channels_(1).loadImage((jj-1)*10+5));
    h1=figure;
    imshow(I1,[]); hold on
    curFlow1 = TavgFlow1{jj};
    quiver(curFlow1(:,2),curFlow1(:,1),mag*(curFlow1(:,4)-curFlow1(:,2)),mag*(curFlow1(:,3)-curFlow1(:,1)),0,'y');
    hgexport(h,strcat(tifPath_Flow1,'/flow1',num2str(jj,iiformat),'.tif'),hgexport('factorystyle'),'Format','tiff')
    hgsave(h,strcat(figPath_Flow1,'/flow1',num2str(jj,iiformat),'.fig'),'-v7.3')
    close(h1)

    I2=double(movieData.channels_(2).loadImage((jj-1)*10+5));
    h2=figure;
    imshow(I2,[]); hold on
    curFlow2 = TavgFlow2{jj};
    quiver(curFlow2(:,2),curFlow2(:,1),mag*(curFlow2(:,4)-curFlow2(:,2)),mag*(curFlow2(:,3)-curFlow2(:,1)),0,'r');
    hgexport(h,strcat(tifPath_Flow2,'/flow2',num2str(jj,iiformat),'.tif'),hgexport('factorystyle'),'Format','tiff')
    hgsave(h,strcat(figPath_Flow2,'/flow2',num2str(jj,iiformat),'.fig'),'-v7.3')
    close(h2)

    h3=figure;
    imshow(I1,[]); hold on
    quiver(curFlow1(:,2),curFlow1(:,1),mag*(curFlow1(:,4)-curFlow1(:,2)),mag*(curFlow1(:,3)-curFlow1(:,1)),0,'y');
    quiver(curFlow2(:,2),curFlow2(:,1),mag*(curFlow2(:,4)-curFlow2(:,2)),mag*(curFlow2(:,3)-curFlow2(:,1)),0,'r');
    hgexport(h,strcat(tifPath_FlowCombined,'/flowComb',num2str(jj,iiformat),'.tif'),hgexport('factorystyle'),'Format','tiff')
    hgsave(h,strcat(figPath_FlowCombined,'/flowComb',num2str(jj,iiformat),'.fig'),'-v7.3')
    close(h3)

    % DCS showing
    dI = double(I)/max(I(:));
    imgDCS = NaN(size(dI));
    curDCS = DCS{jj};
    xMin = round(min(curDCS(:,1)));
    xMax = round(max(curDCS(:,1)));
    yMin = round(min(curDCS(:,2)));
    yMax = round(max(curDCS(:,2)));
    [xi,yi] = meshgrid(xMin:xMax,yMin:yMax);
    DCSinterp = griddata(curDCS(:,1),curDCS(:,2),curDCS(:,3),xi,yi);
%     F = scatteredInterpolant(curDCS(:,1),curDCS(:,2),curDCS(:,3));
%     DCSinterp = F(xi,yi);
    imgDCS(yMin:yMax,xMin:xMax) = DCSinterp;
    imgDCS(~polyOne) = NaN;

    dcsImg = imDataMapOverlay(dI,imgDCS,[0.6 1],cMap);
    h4=figure;
    imshow(dcsImg,[]); hold on
    hgexport(h,strcat(tifPath_DCS,'/DCSMap',num2str(jj,iiformat),'time',num2str(((jj-1)*10+5)*dt),'sec.tif'),hgexport('factorystyle'),'Format','tiff')
    hgsave(h,strcat(figPath_DCS,'/DCSMap',num2str(jj,iiformat),'time',num2str(((jj-1)*10+5)*dt),'sec.fig'),'-v7.3')
    close(h4)

    % VMCS showing
    imgVMCS = NaN(size(dI));
    curVMCS = VMCS{jj};
    VMCSinterp = griddata(curVMCS(:,1),curVMCS(:,2),curVMCS(:,3),  xi,yi);
    imgVMCS(yMin:yMax,xMin:xMax) = VMCSinterp;
    imgVMCS(~polyOne) = NaN;

    vmcsImg = imDataMapOverlay(dI,imgVMCS,[0 2],cMap);
    h5=figure;
    imshow(vmcsImg,[]); hold on
    hgexport(h,strcat(tifPath_VMCS,'/VMCSMap',num2str(jj,iiformat),'time',num2str(((jj-1)*10+5)*dt),'sec.tif'),hgexport('factorystyle'),'Format','tiff')
    hgsave(h,strcat(figPath_VMCS,'/VMCSMap',num2str(jj,iiformat),'time',num2str(((jj-1)*10+5)*dt),'sec.fig'),'-v7.3')
    close(h5)

    % saving DCS and VMCS values for roi samples
    for k=1:nPoly
        roiDCS(jj,k) = imgDCS(round(avgRoiFlow1{k}.pos(1,2)),round(avgRoiFlow1{k}.pos(1,1)));
        roiVMCS(jj,k) = imgVMCS(round(avgRoiFlow1{k}.pos(1,2)),round(avgRoiFlow1{k}.pos(1,1)));
    end
end
%% Plot speed with time per each ROI
h=figure; hold on
for k=1:nPoly
    subplot(nPoly,1,k)
    [hAx,hLine1,hLine2] = plotyy([avgRoiFlow1{k}.t(1:end-3)',avgRoiFlow1{k}.t(1:end-3)'],...
            [avgRoiFlow1{k}.speed(1:end-3)',avgRoiFlow2{k}.speed(1:end-3)'],...
            [avgRoiFlow1{k}.t(5:10:nTavgFrames*10-5)',avgRoiFlow1{k}.t(5:10:nTavgFrames*10-5)'],...
            [roiDCS(:,k), roiVMCS(:,k)]);
    
%     plot(avgRoiFlow2{k}.t(1:end-3), avgRoiFlow2{k}.speed(1:end-3), '--','Color',[1/nPoly*(k-1) 1/nPoly*(k-1) 1/nPoly*(k-1)], 'Linewidth',2)
%     plot(avgRoiFlow1{k}.t(5:10:nTavgFrames*10-5),roiDCS(:,k),'bo-','Linewidth',1)
%     plot(avgRoiFlow1{k}.t(5:10:nTavgFrames*10-5),roiVMCS(:,k),'go-','Linewidth',1)
end
%'Color',[1/nPoly*(k-1) 1/nPoly*(k-1) 1/nPoly*(k-1)],'Linewidth',2

hgexport(h,strcat(tifPath,'/flowMag_time.tif'),hgexport('factorystyle'),'Format','tiff')
hgsave(h,strcat(figPath,'/flowMag_time.fig'),'-v7.3')
save([dataPath '/flowMag_time.mat'],'avgRoiFlow1','avgRoiFlow2','roiDCS','roiVMCS','imgDCS','imgVMCS')
close(h)
end
