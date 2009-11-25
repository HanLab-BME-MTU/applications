function plusTipBatchQuadPlot(groupList,speedLims,speedDiv,lifeLims,lifeDiv,remBegEnd,timeRange,doPlot)
% plusTipBatchQuadPlot creates speed/lifetime quad plots for groups of projects in groupList
%
% SYNOPSIS: plusTipBatchQuadPlot(groupList,speedLims,speedDiv,lifeLims,lifeDiv,remBegEnd,timeRange,doPlot)
%
% INPUT:
%
%   groupList : nx2 cell array containing group identifiers in the first
%               column (see plusTipPickGroups) and project paths in the
%               second
%   speedLims : [lowerSpeedLimit, upperSpeedLimit] (microns/min) for the
%               x-axis
%   speedDiv  : number for division mark along the x-axis for all the quad
%               plots
%   lifeLims  : [lowerLifetimeLimit, upperLifetimeLimit] (sec) for the
%               y-axis
%   lifeDiv   : number for division mark along the y-axis for all the quad
%               plots
%   remBegEnd : 1 to remove data from the beginning/end of the movie, 0 to
%               keep (this option may bias the results since lifetimes at
%               the beginning/end are incomplete)
%   timeRange : [startingFrame endingFrame]
%   doPlot    : 1 to create the quad plots, 0 to make the percentage bars
%               only
%
%
% OUTPUT:
%
%   Sub-directories for each group represented in groupList, containing
%   individual quad plots and overlays for each project in the group. Also,
%   plots showing the percent of red/green/yellow/blue tracks for each
%   movie or for all tracks from the whole group are saved as .fig or .png
%   images.  groupQuadInfo variables are saved containing the input
%   parameters (speedLims, speedThresh=speedDiv,
%   lifeLims,lifeThresh=lifeDiv, remBegEnd) as well as some new data:
%       speedPrct: percentile for each movie of the speed distribution
%                  corresponding to speedThresh
%       lifePrct : percentile for each movie of the lifetime distribution
%                  corresponding to lifeThresh
%       popRYGB  : nx4 matrix containing the number of tracks in the
%                  red/yellow/green/blue populations for all n projects in
%                  the group
%       nPrctRYGB: nx4 matrix containing the corresponding percentages



close all

if ispc
    fileExt='.emf';
else
    fileExt='.jpg';
end

% get the group names
[grpNames,b,m]=unique(groupList(:,1));

% get the output directory
saveDir=uigetdir(pwd,'Choose directory for storing Quadrant Plot data');

% iterate through each group
for iGroup=1:length(grpNames)
    
    saveDirSub=[saveDir filesep grpNames{iGroup}];
    if ~isdir(saveDirSub)
        mkdir(saveDirSub)
    end

    grpIdx=find(m==iGroup);

    % initialize percentile, threshold, and population arrays
    speedPrct=zeros(length(grpIdx),1);
    speedThresh=zeros(length(grpIdx),1);
    lifePrct=zeros(length(grpIdx),1);
    lifeThresh=zeros(length(grpIdx),1);
    popRYGB=zeros(length(grpIdx),4);

    % load data and make plots for each project iSub in group iGroup
    for iSub=1:length(grpIdx)
        % load data for iSub
        projDir=formatPath(groupList{grpIdx(iSub),2});
        load([projDir filesep 'meta' filesep 'projData'])

        % runt the quad plot generator
        [speedPrct(iSub),speedThresh(iSub),lifePrct(iSub),lifeThresh(iSub),...
            popRYGB(iSub,:)]=plusTipParamPlot('growthSpeed',[],speedDiv,'growthLifetime',[],lifeDiv,...
            projData,remBegEnd,timeRange,doPlot,speedLims,lifeLims);

        if doPlot==1
            saveas(1,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_all' '.fig'])
            saveas(1,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_all' fileExt])

            saveas(2,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_fastLong' '.fig'])
            saveas(2,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_fastLong' fileExt])

            saveas(3,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_slowLong' '.fig'])
            saveas(3,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_slowLong' fileExt])

            saveas(4,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_fastShort' '.fig'])
            saveas(4,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_fastShort' fileExt])

            saveas(5,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_slowShort' '.fig'])
            saveas(5,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_slowShort' fileExt])

            saveas(6,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_scatter' '.fig'])
            saveas(6,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_scatter' fileExt])
            
            saveas(7,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_prctBar' '.fig'])
            saveas(7,[saveDirSub filesep groupList{grpIdx(iSub),1} '_' num2str(grpIdx(iSub)) '_prctBar' fileExt])

            close all
        end

    end
    % make the percentage bar for each sub-roi in the same plot
    [nPrctRYGB]=plusTipQuadColorbar(popRYGB);
    titleStr=strrep(groupList{grpIdx(iSub),1},'_','-');
    title(titleStr)
    saveas(gcf,[saveDir filesep groupList{grpIdx(iSub),1} '_prctBarAll' '.fig'])
    saveas(gcf,[saveDir filesep groupList{grpIdx(iSub),1} '_prctBarAll' fileExt])
    close(gcf)

    % sum the populations to make the percentage bar for the whole group
    [nPrctRYGB]=plusTipQuadColorbar(sum(popRYGB,1));
    titleStr=strrep(groupList{grpIdx(iSub),1},'_','-');
    title([titleStr ', N=' num2str(length(grpIdx))])
    saveas(gcf,[saveDir filesep groupList{grpIdx(iSub),1} '_prctBarAll_Merged' '.fig'])
    saveas(gcf,[saveDir filesep groupList{grpIdx(iSub),1} '_prctBarAll_Merged' fileExt])
    close(gcf)
            
    
    save([saveDir filesep groupList{grpIdx(iSub),1} '_groupQuadInfo'],'speedPrct','speedThresh','lifePrct','lifeThresh',...
        'popRYGB','nPrctRYGB','speedLims','lifeLims','remBegEnd','timeRange');

end
