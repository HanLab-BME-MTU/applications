function plusTipBatchQuadPlot(groupList,speedLims,speedDiv,lifeLims,lifeDiv,remBegEnd,timeRange,doPlot)

close all

[grpNames,b,m]=unique(groupList(:,1));

saveDir=uigetdir(pwd,'Choose directory for storing Quadrant Plot data');

for iGroup=1:length(grpNames)
    saveDirSub=[saveDir filesep grpNames{iGroup}];
    if ~isdir(saveDirSub)
        mkdir(saveDirSub)
    end

    grpIdx=find(m==iGroup);

    % initialize stuff
    speedPrct=zeros(length(grpIdx),1);
    speedThresh=zeros(length(grpIdx),1);
    lifePrct=zeros(length(grpIdx),1);
    lifeThresh=zeros(length(grpIdx),1);
    popRYGB=zeros(length(grpIdx),4);

    for iSub=1:length(grpIdx)
        projDir=groupList{grpIdx(iSub),2};
        load([projDir filesep 'meta' filesep 'projData'])

        [speedPrct(iSub),speedThresh(iSub),lifePrct(iSub),lifeThresh(iSub),...
            popRYGB(iSub,:)]=plusTipParamPlot('growthSpeed',[],speedDiv,'growthLifetime',[],lifeDiv,...
            projData,remBegEnd,timeRange,doPlot,speedLims,lifeLims);

        if doPlot==1
            saveas(1,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_all' '.fig'])
            saveas(1,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_all' '.png'])

            saveas(2,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_fastLong' '.fig'])
            saveas(2,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_fastLong' '.png'])

            saveas(3,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_slowLong' '.fig'])
            saveas(3,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_slowLong' '.png'])

            saveas(4,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_fastShort' '.fig'])
            saveas(4,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_fastShort' '.png'])

            saveas(5,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_slowShort' '.fig'])
            saveas(5,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_slowShort' '.png'])

            saveas(6,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_scatter' '.fig'])
            saveas(6,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_scatter' '.png'])
            
            saveas(7,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_prctBar' '.fig'])
            saveas(7,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_prctBar' '.png'])

            close all
        end

    end
    % make the percentage bar for each sub-roi in the same plot
    [nPrctRYGB]=plusTipQuadColorbar(popRYGB);
    titleStr=strrep(groupList{grpIdx(iSub),1},'_','-');
    title(titleStr)
    saveas(gcf,[saveDir filesep groupList{grpIdx(iSub),1} '_prctBarAll' '.fig'])
    saveas(gcf,[saveDir filesep groupList{grpIdx(iSub),1} '_prctBarAll' '.png'])
    close(gcf)

    % sum the populations to make the percentage bar for the whole group
    [nPrctRYGB]=plusTipQuadColorbar(sum(popRYGB,1));
    titleStr=strrep(groupList{grpIdx(iSub),1},'_','-');
    title([titleStr ', N=' num2str(length(grpIdx))])
    saveas(gcf,[saveDir filesep groupList{grpIdx(iSub),1} '_prctBarAll_Merged' '.fig'])
    saveas(gcf,[saveDir filesep groupList{grpIdx(iSub),1} '_prctBarAll_Merged' '.png'])
    close(gcf)
            
    save([saveDir filesep groupList{grpIdx(iSub),1} '_groupQuadInfo'],'speedPrct','speedThresh','lifePrct','lifeThresh',...
        'popRYGB','nPrctRYGB','speedLims','lifeLims','remBegEnd','timeRange');

end
