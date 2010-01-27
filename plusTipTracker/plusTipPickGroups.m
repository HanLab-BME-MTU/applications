function [groupList]=plusTipPickGroups(autoGrp,relDirs,projList,saveResult)
% plusTipPickGroups allows user to select groups of movies
%
% SYNOPSIS: [groupList]=plusTipPickGroups(autoGrp,relDirs,projList,saveResult)
%
% INPUT : user is asked to select the projList file(s) containing projects
%         to be used for group selection.

%         autoGrp (opt)   : if [] (def), user picks the groups and gives
%                           them unique names.  (leave names empty to use
%                           1,2,3... as the group names)
%                           if 1, user is asked to select one or more
%                           categories from the parsed file path of the first
%                           project in projList to define how groups should
%                           be created.  for this to work all projects in
%                           projList need to have corresponding categories
%                           at the same folder levels.
%         relDirs (opt)   : vector containing the directory number(s)
%                           relative to the roi_x folder (roi is 1, one
%                           level up is 2, two levels up is 3, etc.) that
%                           should be used for making the groups 
%                           (eg. [5 4]). this option only works if autoGrp
%                           is 1.
%        saveResult       : 1 to prompt to save, 0 to not do so
%
% OUTPUT: groupList : n x 2 cell array, where n is the number of projects
%                     chosen for all groups. the first column contains the
%                     group name.  the second contains the project path.

groupList=[];

if nargin<1 || isempty(autoGrp) || autoGrp~=1
    autoGrp=[];
end

if nargin<2 || ~isvector(relDirs)
    relDirs=[];
end


% ask user to select projList file and check which movies have been tracked
if nargin<3 || isempty(projList)
    [projList]=combineProjListFiles;
    if isempty(projList)
        return
    end
end
[allProjects]=projList2Cell(projList);

if nargin<4 || isempty(saveResult)
    saveResult=0;
end


% have user select groups of projects
nProj=length(allProjects);
groupList=cell(nProj,2);

% user picks - no auto groups
if isempty(autoGrp)
    countMovie=1;
    countLabel=1;
    pickAgain='yes';
    h = msgbox('Please select first group from the list');
    uiwait(h)
    while strcmpi(pickAgain,'yes')

        % user selection of projects for iGroup
        selection=listSelectGUI(allProjects(:,1),[],'move');

        % get name of group for the legend
        temp=inputdlg({'Enter group name:'},'Input for legend label',1);
        if isempty(temp) % user clicked 'x'
            legendLabel=num2str(countLabel);
        else
            if isempty(temp{1,1}) % user clicked ok but entered nothing
                legendLabel=num2str(countLabel);
            else
                legendLabel=temp{1,1}; % user entered a string
            end
        end

        % record project selection directories and legend names
        for iProj=1:length(selection)
            groupList(countMovie,:)={legendLabel allProjects{selection(iProj),1}};
            countMovie=countMovie+1;
        end

        % get rid of already-picked ones for next round
        allProjects(selection,:)=[];
        
        % ask whether to select another group
        pickAgain=questdlg('Select another group?');
        countLabel=countLabel+1;
    end
    groupList(countMovie:end,:)=[];

else
    % auto groups
    for iProj=1:nProj
       
        currentROI=formatPath(allProjects{iProj,1});
        % parse the path
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

        % add entry for subROIs and subx folders if there are no sub-rois
        if length(words)==roiIdx
            words(roiIdx+1:roiIdx+2)={'noSubROIs';'noSubs'};
        end
        
        % invert order
        words=words(end:-1:1);
        
        % ask user for which categories to use to make the groups
        
        if iProj==1
            if isempty(relDirs)
                h = msgbox('Please select levels to be used for grouping');
                uiwait(h)
                autoGrp=listSelectGUI(words,[],'copy');
                autoGrp=autoGrp(end:-1:1);
            else
                % autoGrp is relative to the roi directory, so add 2 since
                % words is relative to sub directory
                autoGrp=relDirs+2;
            end
        end


        macroWord=[];
        for iNum=1:length(autoGrp)
            macroWord=[macroWord '_' words{autoGrp(iNum)}];
        end
        macroWord(1)=[];
        
        groupList(iProj,:)={macroWord allProjects{iProj,1}};

    end
    
end
if saveResult==1
    temp=inputdlg({'Enter file name:'},'',1);
    dirName=uigetdir(pwd,'Select output directory for groupList.');
    save([dirName filesep temp{1}],'groupList')
end