function [meanLF,meanEV,meanEdgeAdvance] = getLifeTimeAndEdgeVelocity(MD)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
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

%% Get edge velocity
% I'll use mean(edgeAdvanceDist)/lifetime
eV = arrayfun(@(x) mean(x.edgeAdvanceDist)/x.lifeTime,tracksNA);
meanDist = arrayfun(@(x) mean(x.edgeAdvanceDist),tracksNA);
% Get adhesion life time
lifeTime = arrayfun(@(x) x.lifeTime,tracksNA);
% Plot between the two
% Scatter plot
% plot(lifeTime,eV,'ko')
% 2D histogram
%% refinement
% Some filterring is required
% 1. We are not interested in adhesions in static or retracting edge
idValid1 = meanDist>abs(0.5*std(meanDist));
%  2. We are not interested in adhesions that has movie-long lifetime
idValid2 = lifeTime <0.8*MD.nFrames_;
idValid = idValid1 & idValid2;

%% All
figure('Position',[100,100,400,900]), subplot(3,1,1)
plot(lifeTime(idValid),meanDist(idValid),'ko')
hold on
plot(mean(lifeTime(idValid)),mean(meanDist(idValid)),'bx','MarkerSize',20)
% Add barplot by binning data by 10 in x-axis. I can use discretize
% function!
% edges = [0:20:80 100:50:ceil(max(lifeTime(idValid))/50)*50];
edges = 0:20:ceil(max(lifeTime(idValid))/50)*50;
Y = discretize(lifeTime(idValid),edges);
mDValid = meanDist(idValid);
meanDistPerBin = arrayfun(@(x) mean(mDValid(Y==x)),1:numel(edges));
pB = bar(edges+10,meanDistPerBin);
pB.FaceAlpha=0;
pB.EdgeColor='r';
pB.LineWidth = 2;
stdDistPerBin = arrayfun(@(x) stdErrMean(mDValid(Y==x)),1:numel(edges));
eH = errorbar(edges+10,meanDistPerBin,...
    stdDistPerBin,'Color', 'r');
eH.YNegativeDelta = [];
eH.LineStyle='none';
eH.LineWidth = 2;
mPath = MD.getPath;
title(mPath(end-17:end))
xlabel('Lifetime (frames)')
ylabel('Edge advance distance (px)')

subplot(3,1,2)
meanDistBinArray = arrayfun(@(x) mDValid(Y==x),1:numel(edges),'unif',false);
boxPlotCellArray(meanDistBinArray,cellstr(num2str(edges'+10))',1,0,0);
xlabel('Lifetime (frames)')
ylabel('Edge advance distance (px)')


subplot(3,1,3)
plot(lifeTime(idValid),eV(idValid),'ko')
hold on
plot(mean(lifeTime(idValid)),mean(eV(idValid)),'rx','MarkerSize',20)
xlabel('Lifetime (frames)')
ylabel('Edge velocity (px/frame)')

meanLF = mean(lifeTime(idValid));
meanEdgeAdvance = mean(meanDist(idValid));

%% All - eV
meanEV = mean(eV(idValid));
end

