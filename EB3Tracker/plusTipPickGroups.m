function [movDataSet]=plusTipPickGroups
% allow user to group movies that have been analyzed

movDataSet=[];

[result,notDone]=plusTipCheckIfDone;

finalList=result(setdiff(1:size(result,1),notDone),1);

% have user select groups of projects
projGroupDir=cell(size(finalList));
projGroupName=cell(size(finalList));
countMovie=1;
countLabel=1;
pickAgain='yes';
h = msgbox('Please select first group from the list');
uiwait(h)
while strcmpi(pickAgain,'yes')
    
    % user selection of projects for iGroup
    [selection,selectionList]=listSelectGUI(finalList,[],'move');
    
    % get name of group for the legend
    temp=inputdlg({'Enter group name:'},'Input for legend label',1);
    if isempty(temp) % user clicked 'x'
        legendLabel=num2str(iGroup);
    else
        if isempty(temp{1,1}) % user clicked ok but entered nothing
            legendLabel=num2str(countLabel);
        else
            legendLabel=temp{1,1}; % user entered a string
        end
    end

    % record project selection directories and legend names
    projGroupDir(countMovie:countMovie+length(selection)-1,1)=finalList(selection,1);
    for iMov=countMovie:countMovie+length(selection)-1
        projGroupName{iMov}=legendLabel;
    end
    countMovie=countMovie+length(selection);

    % get rid of already-picked ones for next round
    finalList=finalList(setdiff(1:length(finalList),selection),1);
    % ask whether to select another group
    pickAgain=questdlg('Select another group?');
    countLabel=countLabel+1;
end
projGroupDir(countMovie:end)=[];
projGroupName(countMovie:end)=[];

nProj=length(projGroupDir);

for i=1:nProj
    temp=load([projGroupDir{i,1} filesep 'meta' filesep 'projData.mat']);
    temp=temp.projData.typeStats;
    if i==1
        dataNames=fieldnames(temp);
        groupData=cell(nProj,length(dataNames));
    end
    temp=struct2cell(temp)';
    groupData(i,:)=temp;
end
groupData=abs(cell2mat(groupData));
movNum = strcat({'Movie'},num2str((1:nProj)','%d'));
movDataSet=dataset({projGroupName,'groupName'},{groupData(:,1:end-2),'growthSpeedMean','growthSpeedStd','Ppause','Pcat','pauseSpeedMean','pauseSpeedStd','shrinkSpeedMean','shrinkSpeedStd'},{projGroupDir,'directory'},'ObsNames',movNum);

