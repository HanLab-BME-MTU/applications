function [groupListSpeedLifeDisp,grpData]=plusTipGroupListSummary(groupList)
% returns cell array with grouped projects and GROWTH speed, lifetime, and displacement
% data for each

groupListSpeedLifeDisp=cell(size(groupList,1),3);
groupListSpeedLifeDisp(:,1:2)=groupList;
for i=1:length(groupList)

    % load the sub-project data
    projDir=[groupList{i,2} filesep 'meta'];
    cd(projDir)
    load projData
    aT=plusTipMergeSubtracks(projData); % use merged data, don't crop beginning/end
    aT(aT(:,5)~=1,:)=[]; % use growths only
    aT(:,6)=aT(:,6).*projData.secPerFrame;
    aT(:,7)=aT(:,7).*(projData.pixSizeNm/1000);
    
    groupListSpeedLifeDisp{i,3}=aT(:,[4,6:7]);
end

[grpNames,b,m]=unique(groupList(:,1));
grpNames=cellfun(@(x) ['grp_' x],grpNames,'uniformoutput',0);

grpData.allGroups=cell2mat(groupListSpeedLifeDisp(:,3));
for iGroup=1:length(grpNames)
    grpIdx=find(m==iGroup);
    
    grpData.(grpNames{iGroup})=cell2mat(groupListSpeedLifeDisp(grpIdx,3));
end