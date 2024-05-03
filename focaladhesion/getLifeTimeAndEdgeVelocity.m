function [allLF,allEV,PTprot] = getLifeTimeAndEdgeVelocity(MD,askPlotColorBar,ltMax)
%function [allLF,allEV,PTprot] =
%getLifeTimeAndEdgeVelocity(MD,askPlotColorBar) calculate lifetime and edge
%velocity and persistence time.

if nargin<2
    askPlotColorBar = false;
end
%% Load tracks
faPackage=MD.getPackage(MD.getPackageIndex('FocalAdhesionPackage'));
% Load classification process
classProc = faPackage.getProcess(8);
iChan = find(classProc.checkChannelOutput);
finalProc = faPackage.getProcess(11);
iClassObj = load([finalProc.funParams_.OutputDirectory filesep 'data' ...
    filesep 'idGroups.mat'],'idGroups');
iClasses = iClassObj.idGroups;
idGroupLabel = 1*iClasses{1};
for ii=2:9 
    idGroupLabel = idGroupLabel+ ii*iClasses{ii};
end
tracksNA=finalProc.loadChannelOutput(iChan,'output','tracksNA');
%% sampling only adhesions near the edge
distToEdge = arrayfun(@(x) min(x.distToEdge(x.distToEdge>0)),tracksNA);
thresDist = 700/MD.pixelSize_;
tracksNA2 = tracksNA(distToEdge<thresDist);% & idGroupLabel~=6);
%% Get edge velocity
% I'll use mean(edgeAdvanceDist)/lifetime
% eV = arrayfun(@(x) mean(x.edgeAdvanceDist)/x.lifeTime,tracksNA2);
% meanDist = arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNA2);
% Get adhesion life time
lifeTime = arrayfun(@(x) x.lifeTime,tracksNA2);
% Plot between the two
% Scatter plot
% plot(lifeTime,eV,'ko')
% 2D histogram
%% Window-based edge protrusion
% Get the window package
wPack = MD.getPackage(MD.getPackageIndex('WindowingPackage'));
% inspect each process
% pProc = wPack.getProcess(1);
wProc = wPack.getProcess(2);
psProc = wPack.getProcess(3);
psSamples = psProc.loadChannelOutput;
% I got the information. probably in avgNormal would do it.
% Only problem is to find out the location of each window.
% First, find where each adhesion is located
nTracks = numel(tracksNA2);
meanPersTimeProt = zeros(nTracks,1);
meanPersTimeRet = zeros(nTracks,1);
meanEdgeProtVelAll = zeros(nTracks,1);
meanEdgeRetVelAll = zeros(nTracks,1);
% get windows
windows = wProc.loadChannelOutput(1);
nF = size(windows,2);
WindowRange{1} = 4:(nF-3);
nIMFs = 1; % Input # of IMFs to subtract
[cellData, dataSet] =...
    edgeVelocityQuantification_KRC(MD,nIMFs,'includeWin',WindowRange,...
    'scale',true(1,nF),'outLevel',9*ones(1,nF)); %,'fileName',[name, 'Analysis01']); 

for ii=1:nTracks
    % get the location
    sF = tracksNA2(ii).startingFrameExtra;
    eF = tracksNA2(ii).endingFrameExtra;
    % If the adhesion didn't end until the end of the movie, we shouldn't
    % count those into our lifetime quantification except for long-lived adhesions 
    if (eF-sF)<60/MD.timeInterval_ && eF>= MD.nFrames_-2 %2 is a tolerance
        meanPersTimeProt(ii) = NaN;
        meanPersTimeRet(ii) = NaN;
        meanEdgeProtVelAll(ii) = NaN;
        meanEdgeRetVelAll(ii) = NaN;
        continue
    end
        
    curX = tracksNA2(ii).xCoord(sF);
    curY = tracksNA2(ii).yCoord(sF);
    % Now identify which window contains the adhesion location.
    iClosest = findClosestWindow(windows,[curX curY]);
    % Found it! Now let's get the vector
%     curEdgeVel = psSamples.avgNormal(iClosest,sF:eF);
    % Now get the mean velocity and persistent time
%     %little bit of smoothing?
%     tRange = 1:MD.nFrames_;
%     sd_spline= csaps(tRange,curEdgeVel,0.5);
%     curEdgeVel2=ppval(sd_spline,tRange);
%     tRange = tracksNA2(ii).startingFrameExtraExtra:tracksNA2(ii).endingFrameExtraExtra;
%     d = tracksNA2(ii).edgeAdvanceDist(tRange);
%     % if we are collecting the zero edge motion, that would be from the
%     % image boundary. Need to ignore them.
%     if all(d==0)
%         meanPersTimeProt(ii) = NaN;
%         meanPersTimeRet(ii) = NaN;
%         continue
%     end
%         
%     sd_spline= csaps(tRange,d,0.1);
%     sd=ppval(sd_spline,tRange);
%     
%     %Get the velocity
%     curVel = diff(sd);
%     % The same thing. if the edgeAdvanceDist captures the image boundary,
%     % the velocity will come out as zero. Should ignore them.
%     if all(curVel==0)
%         meanPersTimeProt(ii) = NaN;
%         meanPersTimeRet(ii) = NaN;
%         continue
%     end
%     curEdgeVel = psSamples.avgNormal(iClosest,sF:eF);

    if ismember(iClosest,cellData{1}.data.includedWin{1})
        curWinStatProt = cellData{1}.protrusionAnalysis.windows(...
            find(cellData{1}.data.includedWin{1}==iClosest));
        meanPersTimeProt(ii) = mean(curWinStatProt.persTime);
        meanEdgeProtVelAll(ii) = mean(curWinStatProt.Veloc);
        curWinStatRet = cellData{1}.retractionAnalysis.windows(...
            find(cellData{1}.data.includedWin{1}==iClosest));
        meanPersTimeRet(ii) = mean(curWinStatRet.persTime);
        meanEdgeRetVelAll(ii) = mean(curWinStatRet.Veloc);
    else
        meanPersTimeProt(ii) = NaN;
        meanPersTimeRet(ii) = NaN;
        meanEdgeProtVelAll(ii) = NaN;
        meanEdgeRetVelAll(ii) = NaN;
    end

%     [curProt,curRet] = getPersistenceTime(curEdgeVel,MD.timeInterval_); %,'plotYes',true);
% %     [curProt,curRet] = getPersistenceTime(curEdgeVel2,MD.timeInterval_); %,'plotYes',true);
%     % I will report the mean persistent time for protrusion and retraction
%     meanPersTimeProt(ii) = mean(curProt.persTime);
%     meanPersTimeRet(ii) = mean(curRet.persTime);
% %     meanEdgeProtVelAll(ii) = nanmean(curEdgeVel(curEdgeVel>=0));
% %     meanEdgeRetVelAll(ii) = nanmean(curEdgeVel(curEdgeVel<=0));
%     meanEdgeProtVelAll(ii) = nanmean(curProt.Veloc);
%     meanEdgeRetVelAll(ii) = nanmean(curRet.Veloc);
end


%% refinement
% Some filterring is required
% % 1. We are not interested in adhesions in static or retracting edge
% idValid1 = meanDist>abs(0.5*std(meanDist));
% %  2. We are not interested in adhesions that has movie-long lifetime
% idValid2 = lifeTime <0.8*MD.nFrames_;
% idValid = idValid1 & idValid2;
% 1. We are not interested in adhesions in static or retracting edge
% idValid1 = meanEdgeProtVelAll>abs(0.5*nanstd(meanEdgeProtVelAll));
%  2. We are not interested in adhesions that has movie-long lifetime
idValid2 = lifeTime <0.9*MD.nFrames_; % & lifeTime >=0;
% idValid = idValid1 & idValid2;
idValid = idValid2;

%% Vel vs LT
figEVvsLT=figure('Position',[100,100,300,600]); subplot(2,1,1)
% plot(lifeTime(idValid),meanDist(idValid),'ko')
% plot(lifeTime,meanEdgeProtVelAll,'ko')
plot(lifeTime(idValid),meanEdgeProtVelAll(idValid),'ko')

hold on
binSize=20;
[edges,Y]=plotBinnedBarFromArray(lifeTime(idValid),...
    meanEdgeProtVelAll(idValid),binSize);
% plot(mean(lifeTime(idValid)),mean(meanDist(idValid)),'bx','MarkerSize',20)
% Add barplot by binning data by 10 in x-axis. I can use discretize
% function!
mPath = MD.getPath;
title(mPath(end-17:end))
xlabel('Lifetime (frames)')
ylabel('Edge velocity (px/frame)')

subplot(2,1,2)
eVValid = meanEdgeProtVelAll(idValid);
edgeVelBinArray = arrayfun(@(x) eVValid(Y==x),1:numel(edges),'unif',false);
boxPlotCellArray(edgeVelBinArray,cellstr(num2str(edges'))',1,0,0);
xlabel('Lifetime (frames)')
ylabel('Edge velocity (px/frame)')

%% save the figure
savefig(figEVvsLT,[MD.getPath filesep 'edgeVelVsLifetime.fig']);
print(figEVvsLT,[MD.getPath filesep 'edgeVelVsLifetime.eps'], '-depsc')
close(figEVvsLT)

% subplot(3,1,3)
% plot(lifeTime(idValid),eV(idValid),'ko')
% hold on
% plot(mean(lifeTime(idValid)),mean(eV(idValid)),'rx','MarkerSize',20)
% xlabel('Lifetime (frames)')
% ylabel('Edge velocity (px/frame)')

allLF = lifeTime(idValid);
allEV = meanEdgeProtVelAll(idValid);
%% All - eV
% meanEV = mean(eV(idValid));
% It might be necessary to group all high lifetime (80-140 frames) into one
% I will need to convert the edge advance distance into the velocity
% divided by the same (max) life time.
% But I will concentrate on persistence measurement
%% Persistence measurement about the edge fluctuation.
% It went up.
%% plotting persistent time as a function of lifetime
figPTvsLT=figure('Position',[500,100,400,900]); subplot(3,1,1)
LT = lifeTime(idValid);
PTprot = meanPersTimeProt(idValid);
% plot(LT,PTprot,'ko')
% regression
thresLT = 35;
% full series
try
    mdl = fitlm(LT,PTprot,'VarNames',{'Lifetime','Persistent time'});
    mdl.plot
    % Plot for x>35
    idx = LT>thresLT & ~(LT==83 & PTprot>93);
    mdl35 = fitlm(LT(idx),PTprot(idx),'VarNames',{'Lifetime','Persistent time'});
    mdl35.plot
    mdlunder = fitlm(LT(LT<=thresLT),PTprot(LT<=thresLT),'VarNames',{'Lifetime','Persistent time'});
    hold on;
    mdlunder.plot
catch
    plot(LT,PTprot,'ko')
    hold on
    plotBinnedBarFromArray(LT,PTprot,binSize)
    xlabel('Lifetime (frames)')
    ylabel('Persistance time (sec)')
    title('Protrusion persistance')
end
%%
hold on
% plot(lifeTime(idValid),meanPersTimeRet(idValid),'ro')

subplot(3,1,2)
meanPTPValid = meanPersTimeProt(idValid);
meanPTP_BinArray = arrayfun(@(x) meanPTPValid(Y==x),1:numel(edges),'unif',false);
boxPlotCellArray(meanPTP_BinArray,cellstr(num2str(edges'+10))',1,0,0);
xlabel('Lifetime (frames)')
ylabel('Persistance time (sec)')
title('Protrusion persistance')

subplot(3,1,3)
meanPTRValid = meanPersTimeRet(idValid);
meanPTR_BinArray = arrayfun(@(x) meanPTRValid(Y==x),1:numel(edges),'unif',false);
boxPlotCellArray(meanPTR_BinArray,cellstr(num2str(edges'+10))',1,0,0);
xlabel('Lifetime (frames)')
ylabel('Persistance time (sec)')
title('Retraction persistance')
%% save the figure
savefig(figPTvsLT,[MD.getPath filesep 'PersistenceVsLifetime.fig']);
print(figPTvsLT,[MD.getPath filesep 'PersistenceVsLifetime.eps'], '-depsc')
close(figPTvsLT)
%% Need to look at where they are:
%% cell and cell boundary
iFrame = 100;
% imshow(MD.channels_(iChan).loadImage(iFrame),[])
iEdgeFrames = 1:iFrame; colorMap = @winter;
[figHan,ax1] = makeColoredEdgeOverlayFigure(MD,iFrame,iEdgeFrames,iChan,colorMap);
%colorbar for protrusion
nFramesSel = numel(iEdgeFrames);
frameCols = colorMap(nFramesSel);
if askPlotColorBar
    disp('Draw a rectangle to show a colorbar for cell edges')
    roi = drawrectangle;
    fWidth = floor(ax1.XLim(2));
    fHeight = floor(ax1.YLim(2));
    myPos = [roi.Position(1)/fWidth, (fHeight-sum(roi.Position([2 4])))/fHeight, ...
        0.05 roi.Position(4)/fHeight];
    roi.Visible = 'off';
    
    ax2  = axes('Position',myPos);
else
    ax2  = axes('Position',[0.80  0.2 0.05  0.6]);
end
cBarRGD = ones(length(frameCols),5,3);
cBarRGD(:,:,1) = repmat(frameCols(:,1),1,5);
cBarRGD(:,:,2) = repmat(frameCols(:,2),1,5);
cBarRGD(:,:,3) = repmat(frameCols(:,3),1,5);
imshow(flip(cBarRGD))
% tick
hold on
text(5,iEdgeFrames(end),num2str((iEdgeFrames(1)-1)*MD.timeInterval_/60),'Color','w','FontSize',7);
text(5,iEdgeFrames(1),num2str((iEdgeFrames(end)-1)*MD.timeInterval_/60),'Color','w','FontSize',7);
text(10,50,'Edges (min)','Color','w',...
    'HorizontalAlignment','center','Rotation',90,'FontSize',7);

%% Now let's look at where adhesions in each group are.
axes(ax1)

idxG1 = idValid & lifeTime>thresLT & ~(lifeTime==83 & meanPersTimeProt>93);
idxG2 = idValid & lifeTime<=thresLT;

% now showinig the tracks
if nargin<3
    ltMax = 0.6*MD.nFrames_;
end
ltMaxFrame = floor(ltMax*60/MD.timeInterval_);
hG1 = drawTracksColormap(tracksNA2(idxG1), iFrame, 'lifeTime',[1 ltMaxFrame],spring);
%% short ones
% axes(ax1)
hG2 = drawTracksColormap(tracksNA2(idxG2), iFrame, 'lifeTime',[1 ltMaxFrame],spring);
%% scale bar
scaleBarLength = 5000; %nm
line([5 5+scaleBarLength/MD.pixelSize_],[10 10],'Color','w','LineWidth',3)
text(5+scaleBarLength/MD.pixelSize_/2,20,[num2str(scaleBarLength/1000) ' \mum'],...
    'HorizontalAlignment','center','Color','w','FontSize',7);
%% colorbar for lifetime
nLT = numel(1:ltMaxFrame);
colorMap2 = @spring;
frameCols2 = colorMap2(nLT);
if askPlotColorBar
    disp('Draw a rectangle to show a colorbar for lifetime')
    roi = drawrectangle;
    fWidth = floor(ax1.XLim(2));
    fHeight = floor(ax1.YLim(2));
    myPos = [roi.Position(1)/fWidth, (fHeight-sum(roi.Position([2 4])))/fHeight, ...
        0.05 roi.Position(4)/fHeight];
    roi.Visible = 'off';
    
    ax3  = axes('Position',myPos);
else
    ax3  = axes('Position',[0.90  0.2 0.05  0.6]);
end
cBarRGD2 = ones(length(frameCols2),5,3);
cBarRGD2(:,:,1) = repmat(frameCols2(:,1),1,5);
cBarRGD2(:,:,2) = repmat(frameCols2(:,2),1,5);
cBarRGD2(:,:,3) = repmat(frameCols2(:,3),1,5);
imshow(flip(cBarRGD2))
% tick
hold on
text(6,ltMaxFrame,num2str(0*MD.timeInterval_/60),'Color','w','FontSize',7);
text(6,1,num2str((ltMaxFrame-1)*MD.timeInterval_/60),'Color','w','FontSize',7);
text(10,ltMaxFrame/2,'Lifetime (min)','Color','w',...
    'HorizontalAlignment','center','Rotation',90,...
    'FontSize',7);
pp=0;
%% bring back the bar for the edges
axes(ax2)
axes(ax3)
%% save the figure
savefig(figHan,[MD.getPath filesep 'edge_adhesions.fig']);
print(figHan,[MD.getPath filesep 'edge_adhesions.eps'], '-depsc')
close(figHan)
%% Now let's look at one by one 
indexG1 = find(idxG1');
pp=pp+1;
ii = indexG1(pp);
% for ii=find(idxG1')
    sF = tracksNA2(ii).startingFrameExtra;
    eF = tracksNA2(ii).endingFrameExtra;
    iEdgeFrames = sF:eF;
%     [figHanAll(ii),ax1all(ii)] = ...
    [figHanAll,ax1all] = ...
        makeColoredEdgeOverlayFigure(MD,iFrame,iEdgeFrames,iChan,colorMap);
    drawTracksColormap(tracksNA2(ii), iFrame, 'lifeTime',[1 ltMax],hot);
%     plot(tracksNA2(ii).closestBdPoint(sF:eF,1),...
%         tracksNA2(ii).closestBdPoint(sF:eF,2),'w.');
    windows = wProc.loadChannelOutput(iFrame);
    % Now identify which window contains the adhesion location.
    curX  = tracksNA2(ii).xCoord(sF);
    curY  = tracksNA2(ii).yCoord(sF);
    iClosest = findClosestWindow(windows,[curX curY]);
    plotWindows(windows{iClosest}{1});
    %Get the distance
    distToW = cellfun(@(x) dist2Pt(x',[curX curY]), windows{iClosest}{1},'unif',false);
    [~,indMinDist]=cellfun(@min,distToW);
    % we want to choose the distance to the very edge of the cell mask,
    % which will be not the minimum
    [distBds,indDistBds] = max(cellfun(@min,distToW));
    closestBDPoint = windows{iClosest}{1}{indDistBds}(:,indMinDist(indDistBds));
    plot(closestBDPoint(1),closestBDPoint(2),'w.');
    text(tracksNA2(ii).xCoord(eF)+2,tracksNA2(ii).yCoord(eF),...
        {['LT: ' num2str(tracksNA2(ii).lifeTime)];
        ['Protrusion vel: ' num2str(meanEdgeProtVelAll(ii))];
        ['Persistent time: ' num2str(meanPersTimeProt(ii))];
        ['Distance to edge: ' num2str(distBds)];
        ['ID: ' num2str(ii)]}, 'Color','w');
    %% save the figure
    savefig(figHanAll,[MD.getPath filesep 'one_specific_adhesion.fig']);
    print(figHanAll,[MD.getPath filesep 'one_specific_adhesion.eps'], '-depsc')
    close(figHanAll)
% end

end

