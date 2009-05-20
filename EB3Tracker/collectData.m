function collectData

[result,notDone]=plusTipCheckIfDone;

finalList=result(setdiff(1:size(result,1),notDone),1);

% have user select groups of projects
iGroup=1;
pickAgain='yes';
h = msgbox('Please select first group from the list');
uiwait(h)
while strcmpi(pickAgain,'yes')
    % user selection of projects for iGroup
    [selection,selectionList]=listSelectGUI(finalList,[],'move');
    % get name of group for the legend
    temp=inputdlg({'Enter group name:'},'Input for legend label',1);

    if isempty(temp) % user clicked 'x'
        legendLabel{iGroup,1}=num2str(iGroup);
    else
        if isempty(temp{1,1}) % user clicked ok but entered nothing
            legendLabel{iGroup,1}=num2str(iGroup);
        else
            legendLabel{iGroup,1}=temp{1,1}; % user entered a string
        end
    end
    % record project selection
    projList{iGroup,1}=finalList(selection,1);
    % get rid of already-picked ones for next round
    finalList=finalList(setdiff(1:length(finalList),selection),1);
    % ask whether to select another group
    pickAgain=questdlg('Select another group?');
    iGroup=iGroup+1;
end



% load first project to get field names
temp=load([finalList{1,1} filesep 'meta' filesep 'projData.mat']);
temp=temp.projData.typeStats;
names=fieldnames(temp);



[selection,selectionList]=listSelectGUI(vertcat(names,namesTypeStats),[],'move');



% loop through groups to collect data
for iGroup=1:nGroups
    subList=projList{iGroup,1};
    nProj=length(subList);
    
    for i=1:nProj
        temp=load([subList{i,1} filesep 'meta' filesep 'projData.mat']);
        temp=temp.projData;
        if i==1
            data=cell(nProj,length(names));
        end
        temp=struct2cell(temp)';
        data(i,:)=temp;
    end

    groupData{iGroup,1}=data;

end




cd(homeDir)