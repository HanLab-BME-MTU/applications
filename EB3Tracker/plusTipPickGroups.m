function [projGroupDir,projGroupName,plusTipInfo]=plusTipPickGroups(saveResult,autoPick)
% plusTipPickGroups allows user to select groups of movies
%
% SYNOPSIS: [projGroupDir,projGroupName]=plusTipPickGroups
%
% INPUT : user is asked to select the projList file(s) containing projects
%         to be used for group selection.
%
% OUTPUT: projGroupDir : cell array containing paths to chosen projects
%         projGroupName: cell array containing the group name for each
%                        project



if nargin<1 || isempty(saveResult) || saveResult~=1
    saveResult=0;
    saveDir=[];
else
    saveDir=uigetdir(pwd,'Select output directory for groupList.');
end

if nargin<2 || isempty(autoPick) || autoPick<1
    autoPick=0;
else
    autoPick=autoPick+3; % add to index for extra columns
end

% ask user to select projList file and check which movies have been tracked
[allProjects,notDone]=plusTipCheckIfDone(2);

% show only the ones that have been tracked in the selection box
allProjects(notDone,:)=[];
allProjects=allProjects(:,1);

% have user select groups of projects
projGroupDir=cell(size(allProjects));
projGroupName=cell(size(allProjects));
plusTipInfo=[];

if autoPick==0
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

else
    nMovs = size(allProjects,1);

    % add some extra rows for extra fields
    nCols=11;
    statNamesInfo=cell(1,nCols);
    plusTipInfo=cell(nMovs,nCols);

    for iMov=1:nMovs
        
        % GET PROJECT INFO
        plusTipInfo{iMov,1}='autoGrp'; % group label
        statNamesInfo{1}='label1';
        plusTipData{iMov,2}=''; % secondary label
        statNamesInfo{2}='label2';

        currentROI=formatPath(allProjects{iMov,1});
        % parse the path to get "words" used to identify target, oligo,
        % movie, and roi
        nChar=length(currentROI);
        if ispc
            filesepLoc=regexp(currentROI,'\\');
        else
            filesepLoc=regexp(currentROI,'\/');
        end
        wordStart=[1 filesepLoc+1]; wordEnd=[filesepLoc-1 nChar];
        words=cell(length(wordStart),1);
        for iWord=1:length(wordStart)
            words{iWord,1}=currentROI(wordStart(iWord):wordEnd(iWord));
        end
        % index of the cell which contains parts of the directory name
        roiIdx=find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'roi')),words,'uniformoutput',0)));

        % add entry for subroi if there is none
        if length(words)==roiIdx
            words{roiIdx+1}='';
        end

        % add last 8 bits of directory path info to plusTipInfo
        c=1;
        for i=3:10
            statNamesInfo{i}=['dir' num2str(i-2)];
            if i>length(words)
                plusTipInfo{iMov,i}='';
            else
                plusTipInfo{iMov,i}=words{end-c+1};
            end
            c=c+1;
        end

        statNamesInfo{11}='anDir';
        plusTipInfo{iMov,11}=allProjects{iMov,1};

    end
    
    projGroupDir=allProjects;
    projGroupName=plusTipInfo(:,autoPick);
end

if saveResult==1
    fileName='groupList';
    if autoPick==0
    nameList=unique(projGroupName);
    for i=1:length(nameList)
        fileName=[fileName '_' nameList{i}];
    end
    else
        fileName=[fileName '_' num2str(autoPick-3)];
    end

    save([saveDir filesep fileName],'projGroupName','projGroupDir','plusTipInfo')
end



