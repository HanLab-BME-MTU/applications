function globStats=testTrajectories(trajectoryData)
%TESTTRAJECTORIES produces discrimination matrices for the comparison of trajectoryData
%
% INPUT trajectoryData: Structure with multiple runs from trajectoryAnalysis
%
% OUTPUT globStats: discrimination matrices with p-values
%
% help created 3/05
% c: jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%============
% TEST INPUT
%============

% get info from input

% number of groups to compare
nGroups = length(trajectoryData(:));


 
% %         % currently, we do not need distanceSigma - therefore allow taking
% %         % from individualStatistics
% %         if isfield(trajectoryData(iGroup).individualStatistics(ti),'addedStats')
% %             compStruct(ti).distance = trajectoryData(iGroup).individualStatistics(ti).addedStats.distance(:,1);
% %         else
% %             % do not worry about the last distance we're not taking into
% %             % account. This is just a hack, anyway
% %             compStruct(ti).distance = trajectoryData(iGroup).individualStatistics(ti).dataListSeed(:,11);
% %         end
% % 
% %         % distance variation - distance minus mean of distance
% %         compStruct(ti).distanceSigma = ...
% %             (compStruct(ti).distance - trajectoryData(iGroup).individualStatistics(ti).summary.distanceMean(1));

%====================
% COMPARE GLOBALLY
%====================

% read values. Transform speeds for kstests
for i = 1:nGroups

    globStats(i).growthSpeeds = repeatEntries(...
        trajectoryData(i).overallDistribution.antipolewardSpeed(:,1),...
        trajectoryData(i).overallDistribution.antipolewardSpeed(:,2));
    
    globStats(i).growthSpeedsT = ...
        (globStats(i).growthSpeeds - mean(globStats(i).growthSpeeds));
    
    globStats(i).shrinkageSpeeds = abs(repeatEntries(...
        trajectoryData(i).overallDistribution.polewardSpeed(:,1),...
        trajectoryData(i).overallDistribution.polewardSpeed(:,2)));
    globStats(i).shrinkageSpeedsT = ...
        (globStats(i).shrinkageSpeeds - mean(globStats(i).shrinkageSpeeds));
    
    globStats(i).growthTimes = trajectoryData(i).overallDistribution.growthTime(:,1);
    globStats(i).shrinkageTimes = trajectoryData(i).overallDistribution.shrinkageTime(:,1);
    globStats(i).distance = trajectoryData(i).overallDistribution.distance(:,1);
%     globStats(i).distanceSigma = catStruct(1,'cs.distanceSigma');
end


    

globCompare = discriminationMatrix(globStats);

% compare g/s speeds
for i=1:nGroups
    [dummy,pValue] = ttest2(...
        globStats(i).growthSpeeds,globStats(i).shrinkageSpeeds,...
        0.05,'both','unequal');
    globCompare.speeds(i,1) = pValue;
    [dummy,pValue] = kstest2(...
        globStats(i).growthSpeeds,globStats(i).shrinkageSpeeds);
    globCompare.speeds(i,2) = pValue;
end

[cmap,cLimits]=logColormap;

% display
uH = uiViewPanel;
set(uH,'Name','Growth')
imshow(-log10(globCompare.growthSpeeds));
set(uH,'Colormap',cmap);
set(gca,'CLim',cLimits)
for i=1:nGroups
    for j=1:nGroups
        text(j,i,sprintf('%3.3f',globCompare.growthSpeeds(i,j)),...
            'HorizontalAlignment','center')
    end
end
uH = uiViewPanel;
set(uH,'Name','Shrinkage')
imshow(-log10(globCompare.shrinkageSpeeds));
set(uH,'Colormap',cmap);
set(gca,'CLim',cLimits)
for i=1:nGroups
    for j=1:nGroups
        text(j,i,sprintf('%3.3f',globCompare.shrinkageSpeeds(i,j)),...
            'HorizontalAlignment','center')
    end
end
uH = uiViewPanel;
set(uH,'Name','Speeds')
imshow(-log10(globCompare.speeds));
set(uH,'Colormap',cmap);
set(gca,'CLim',cLimits)
for i=1:nGroups
    for j=1:nGroups
        text(j,i,sprintf('%3.3f',globCompare.speeds(i,j)),...
            'HorizontalAlignment','center')
    end
end
uH = uiViewPanel;
set(uH,'Name','Catastrophe')
imshow(-log10(globCompare.growthTimes));
set(uH,'Colormap',cmap);
set(gca,'CLim',cLimits)
for i=1:nGroups
    for j=1:nGroups
        text(j,i,sprintf('%3.3f',globCompare.growthTimes(i,j)),...
            'HorizontalAlignment','center')
    end
end
uH = uiViewPanel;
set(uH,'Name','Catastrophe')
imshow(-log10(globCompare.shrinkageTimes));
set(uH,'Colormap',cmap);
set(gca,'CLim',cLimits)
for i=1:nGroups
    for j=1:nGroups
        text(j,i,sprintf('%3.3f',globCompare.shrinkageTimes(i,j)),...
            'HorizontalAlignment','center')
    end
end
uH = uiViewPanel;
set(uH,'Name','GlobalDistance')
imshow(-log10(globCompare.distance));
set(uH,'Colormap',cmap);
set(gca,'CLim',cLimits)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cmap,cLimits]=logColormap

colors = [1 1 1;... % boundary
    0   0   0;... % 1
    1   1   1  ;... % not sig
    1   1   0.6;... % p<0.1
    1   0.6 0  ;... % p<0.05
    1   0   0  ];   % p<0.001
cmap = repeatEntries(colors,[1;1;20;6;34;2]);
cLimits = [-0.1,3.05];

