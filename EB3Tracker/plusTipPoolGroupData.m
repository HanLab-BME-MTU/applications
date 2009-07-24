function [groupData] = plusTipPoolGroupData

[projGroupDir,projGroupName]=plusTipPickGroups;

histDir = uigetdir(pwd,'Please select output directory.');

[grpNames,m,movGroupIdx] = unique(projGroupName);


mCount=1; allGrwthSpdCell=cell(1,length(projGroupName));
for iGroup = 1:length(grpNames)
    growthSpeeds=[];
    meanStdGrowth=[];
    growthLifes=[];
    pauseSpeeds=[];
    pauseLifes=[];
    shrinkSpeeds=[];
    shrinkLifes=[];
    growthSpeeds=[];
    meanStdGrowth=[];
    data=[];
    
    % these movies are in iGroup - get info matrix into cell array
    tempIdx=find(movGroupIdx==iGroup);
    for iMov = 1:length(tempIdx)
        temp = load([projGroupDir{tempIdx(iMov)} filesep 'meta' filesep 'projData']);
        data{iMov,1} = temp.projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix;
    end

    growthSpeeds = cellfun(@(x) x(x(:,5)==1,4),data,'uniformoutput',0);
    % get movie-specific mean/std for growth speeds
    meanStdGrowth = cell2mat(cellfun(@(x) [mean(x) std(x)],growthSpeeds,'uniformoutput',0));
    % collect growth speeds for all movies (not just this group)
    allGrwthSpdCell(1,mCount:mCount+length(tempIdx)-1) = growthSpeeds';
    % turn cell array into matrix
    growthSpeeds = cell2mat(growthSpeeds);
    growthLifes = cell2mat(cellfun(@(x) x(x(:,5)==1,6),data,'uniformoutput',0));

    pauseSpeeds = cell2mat(cellfun(@(x) x(x(:,5)==2,4),data,'uniformoutput',0));
    pauseLifes = cell2mat(cellfun(@(x) x(x(:,5)==2,6),data,'uniformoutput',0));

    shrinkSpeeds = abs(cell2mat(cellfun(@(x) x(x(:,5)==3,4),data,'uniformoutput',0)));
    shrinkLifes = cell2mat(cellfun(@(x) x(x(:,5)==3,6),data,'uniformoutput',0));


    groupData.(grpNames{iGroup,1}).growthSpeeds = growthSpeeds;
    groupData.(grpNames{iGroup,1}).growthLifetimes = growthLifes.*temp.projData.secPerFrame;
    groupData.(grpNames{iGroup,1}).meanStdGrowth = meanStdGrowth;
    
    groupData.(grpNames{iGroup,1}).pauseSpeeds = pauseSpeeds;
    groupData.(grpNames{iGroup,1}).pauseLifetimes = pauseLifes.*temp.projData.secPerFrame;

    groupData.(grpNames{iGroup,1}).shrinkSpeeds = shrinkSpeeds;
    groupData.(grpNames{iGroup,1}).shrinkLifetimes = shrinkLifes.*temp.projData.secPerFrame;

    % create x-axis bins spanning all costs in sample
    n=linspace(min([growthSpeeds;pauseSpeeds;shrinkSpeeds]),max([growthSpeeds;pauseSpeeds;shrinkSpeeds]),25);

    % bin the samples
    [x1,nbins1] = histc(growthSpeeds,n); % forward
    [x2,nbins2] = histc(pauseSpeeds,n); % backward
    [x3,nbins3] = histc(shrinkSpeeds,n); % backward

    M=nan(max([length(x1) length(x2) length(x3)]),3);
    M(1:length(x1),1)=x1;
    M(1:length(x2),2)=x2;
    M(1:length(x3),3)=x3;

    % make a stacked histogram of growth, pause, and shrinkage
    figure
    bar(n,M,'stack')
    colormap([1 0 0; 0 0 1; 0 1 0])
    legend('growth','pause','shrinkage','Location','best')
    title(['Sub-Track Speed Distribution: ' grpNames{iGroup,1}])
    xlabel('speed (um/min)');
    ylabel('frequency of tracks');
    saveas(gcf,[histDir filesep 'stackedHist_' grpNames{iGroup,1} '.fig'])
    saveas(gcf,[histDir filesep 'stackedHist_' grpNames{iGroup,1} '.tif'])

    % make a growth speed histogram
    [x1 x2]=hist(growthSpeeds,25);
    bar(x2,x1,'r')
    title('growth speed distribution')
    xlabel('speed (um/min)');
    ylabel('frequency of tracks');

    saveas(gcf,[histDir filesep 'growthHist_' grpNames{iGroup,1} '.fig'])
    saveas(gcf,[histDir filesep 'growthHist_' grpNames{iGroup,1} '.tif'])
    
    % update the counter
    mCount = mCount+length(tempIdx);

end



% get nMovie-vector with number of growth trajectories 
maxSize=cellfun(@(x) length(x),allGrwthSpdCell);
% convert cell array of growth speeds into nTraj x nMovie matrix
allGrwthSpdMatrix=nan(max(maxSize'),length(allGrwthSpdCell));
for i=1:length(allGrwthSpdCell)
   allGrwthSpdMatrix(1:maxSize(i),i) = allGrwthSpdCell{1,i}; 
end
% for each value, put group name into matrix
condition=repmat(projGroupName',[max(maxSize'),1]);
% for each value, put movie number (within group) into matrix
[temp,nOccur]=countEntries(movGroupIdx);
movNum=cell2mat(arrayfun(@(x) [1:x]',nOccur,'uniformOutput',0));
movNumMat=repmat(movNum',[max(maxSize'),1]);

% make box plot comparing all the movies
figure
boxplot(allGrwthSpdMatrix(:),{condition(:) movNumMat(:)},'notch','on','orientation','horizontal');
title('Growth Speeds by Movie')
set(gca,'YDir','reverse')
xlabel('microns/minute')
saveas(gcf,[histDir filesep 'boxplotByMovie.fig'])
saveas(gcf,[histDir filesep 'boxplotByMovie.tif'])

% pool data at the level of the group and make box plot
figure
boxplot(allGrwthSpdMatrix(:),{condition(:)},'notch','on','orientation','horizontal');
title('Growth Speeds by Group')
set(gca,'YDir','reverse')
xlabel('microns/minute')
saveas(gcf,[histDir filesep 'boxplotByGroup.fig'])
saveas(gcf,[histDir filesep 'boxplotByGroup.tif'])
movIdx=[];
movIdx=cell(length(projGroupName),3);
movIdx(:,1)=projGroupName;
movIdx(:,2)=num2cell(movNum,length(projGroupName));
movIdx(:,3)=projGroupDir;
groupData.movIdx=movIdx;
save([histDir filesep 'groupData'],'groupData');
