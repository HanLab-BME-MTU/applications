function [groupData] = plusTipPoolGroupData

[projGroupDir,projGroupName]=plusTipPickGroups;

histDir = uigetdir(pwd,'Please select output directory.');

[grpNames,m,movGroupIdx] = unique(projGroupName);

for iGroup = 1:length(grpNames)
    tempIdx=find(movGroupIdx==iGroup);

    for iMov = 1:length(tempIdx)
        temp = load([projGroupDir{tempIdx(iMov)} filesep 'meta' filesep 'projData']);
        data{iMov,1} = temp.projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix;
    end

    growthSpeeds = cell2mat(cellfun(@(x) x(x(:,5)==1,4),data,'uniformoutput',0));
    growthLifes = cell2mat(cellfun(@(x) x(x(:,5)==1,6),data,'uniformoutput',0));

    pauseSpeeds = cell2mat(cellfun(@(x) x(x(:,5)==2,4),data,'uniformoutput',0));
    pauseLifes = cell2mat(cellfun(@(x) x(x(:,5)==2,6),data,'uniformoutput',0));

    shrinkSpeeds = abs(cell2mat(cellfun(@(x) x(x(:,5)==3,4),data,'uniformoutput',0)));
    shrinkLifes = cell2mat(cellfun(@(x) x(x(:,5)==3,6),data,'uniformoutput',0));


    groupData.(grpNames{iGroup,1}).growthSpeeds = growthSpeeds;
    groupData.(grpNames{iGroup,1}).growthLifetimes = growthLifes.*temp.projData.secPerFrame;

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

    % make the plot
    figure
    bar(n,M,'stack')
    colormap([1 0 0; 0 0 1; 0 1 0])
    legend('growth','pause','shrinkage','Location','best')
    title(['Sub-Track Speed Distribution: ' grpNames{iGroup,1}])
    xlabel('speed (um/min)');
    ylabel('frequency of tracks');
    saveas(gcf,[histDir filesep 'stackedHist_' grpNames{iGroup,1} '.fig'])
    saveas(gcf,[histDir filesep 'stackedHist_' grpNames{iGroup,1} '.tif'])

end
save([histDir filesep 'groupData'],'groupData');