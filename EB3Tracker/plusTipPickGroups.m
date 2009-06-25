function [projGroupDir,projGroupName]=plusTipPickGroups
% plusTipPickGroups allows user to select groups of movies

% INPUT : user is asked to select the projList file(s) containing projects
%         to be used for group selection.  
% OUTPUT: projGroupDir : cell array containing paths to chosen projects
%         projGroupName: cell array containing the group name for each
%                        project




% ask user to select projList file and check which movies have been tracked
[allProjects,notDone]=plusTipCheckIfDone(2);

% show only the ones that have been tracked in the selection box
allProjects(notDone,:)=[];
allProjects=allProjects(:,1);

% have user select groups of projects
projGroupDir=cell(size(allProjects));
projGroupName=cell(size(allProjects));
countMovie=1;
countLabel=1;
pickAgain='yes';
h = msgbox('Please select first group from the list');
uiwait(h)
while strcmpi(pickAgain,'yes')
    
    % user selection of projects for iGroup
    [selection,selectionList]=listSelectGUI(allProjects,[],'move');
    
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
    projGroupDir(countMovie:countMovie+length(selection)-1,1)=allProjects(selection,1);
    for iMov=countMovie:countMovie+length(selection)-1
        projGroupName{iMov}=legendLabel;
    end
    countMovie=countMovie+length(selection);

    % get rid of already-picked ones for next round
    allProjects=allProjects(setdiff(1:length(allProjects),selection),1);
    % ask whether to select another group
    pickAgain=questdlg('Select another group?');
    countLabel=countLabel+1;
end
projGroupDir(countMovie:end)=[];
projGroupName(countMovie:end)=[];


