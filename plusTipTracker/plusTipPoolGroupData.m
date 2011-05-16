function [groupData]=plusTipPoolGroupData(groupList,saveDir,doBtw,doWtn,doPlot,remBegEnd)
% plusTipPoolGroupData pools plus tip data from multiple projects in groups
%
% SYNOPSIS:  [groupData]=plusTipPoolGroupData(groupList,saveDir,doBtw,doWtn,doPlot,remBegEnd)
%
% INPUT:
% groupList : output of plusTipPickGroups, nProj x 2 cell array where the
%             first column contains the group identifier for each group,
%             and the second column contains the project path
% saveDir   : path to output directory
% doBtw     : 1 to pool data from projects in the same group, 0 to skip
% doWtn     : 1 to save text file of distributions from all projects in
%             each group in a "withinGroupComparisons" folder
% doPlot    : 1 to make histograms and boxplots for within and/or between
%             group data
% remBegEnd : 1 to remove tracks existing at the beginning
%             or end of the movie
%
% OUTPUT:
% groupData : structure containing group information and fields for 9
%             distributions: growth speed (gs), fgap speed (fs), bgap speed
%             (bs), growth lifetime (gl), fgap lifetime (fl), bgap lifetime
%             (bl), growth displacement (gd), fgap displacement (fd), and
%             bgap displacement (bd).


homeDir=pwd;

if nargin<1 || isempty(groupList)
    [groupList]=combineGroupListFiles;
end

if nargin<2 || isempty(saveDir)
    saveDir=uigetdir(pwd,'Select output directory for pooled group data.');
end

if nargin<3 || isempty(doBtw)
    doBtw=1;
end

if nargin<4 || isempty(doWtn)
    doWtn=1;
end

if nargin<5 || isempty(doPlot)
    doPlot=1;
end

% assume we should use all data
if nargin<6 || isempty(remBegEnd)
    remBegEnd=0;
end

projGroupName=groupList(:,1);
projGroupDir=cellfun(@(x) formatPath(x),groupList(:,2),'uniformoutput',0);


% fix the names if there are spaces or hyphens and append prefix 'grp'
projGroupName=cellfun(@(x) strrep(x,'-','_'),projGroupName,'uniformoutput',0);
projGroupName=cellfun(@(x) strrep(x,' ','_'),projGroupName,'uniformoutput',0);
projGroupName=cellfun(@(x) ['grp_' x],projGroupName,'uniformoutput',0);


% count unique groups and keep them in order of the original list
[btwGrpNames,m,projGroupIdx] = unique(projGroupName);
[b,idx]=sort(m);
btwGrpNames=btwGrpNames(idx);


projCount=1; % all-project counter
projNum=cell(length(projGroupName),1); % within-group counter
allDataCell=cell(length(btwGrpNames),1); % cell-array for data matrix
groupData=struct('info',{},...
    'gs',{},'fs',{},'bs',{},...
    'gl',{},'fl',{},'bl',{},...
    'gd',{},'fd',{},'bd',{});

for iGroup = 1:length(btwGrpNames)

    % indices of projects in iGroup
    tempIdx=strmatch(btwGrpNames(iGroup),projGroupName,'exact');

    dataByProject=cell(length(tempIdx),1);
    trkCount=1;
    for iProj = 1:length(tempIdx)

        temp = load([projGroupDir{tempIdx(iProj)} filesep 'meta' filesep 'projData']);
        if remBegEnd==1
            % this output has data at beginning/end removed and units
            
            % already converted
            [dummy1,dummy2,dataMat]=plusTipMergeSubtracks(temp.projData);
        else
            % this output just gives merged tracks without converting units
            % or removing beginning/end data
            [dataMat,dummy1,dummy2]=plusTipMergeSubtracks(temp.projData);
            dataMat(:,6)=dataMat(:,6).* temp.projData.secPerFrame; % convert lifetimes to seconds
            dataMat(:,7)=dataMat(:,7).*(temp.projData.pixSizeNm/1000); % convert displacements to microns
        end

        % reassign the track numbers so when combined from multiple projects they don't repeat
        trkIdx=unique(dataMat(:,1));
        dataMat(:,1)=swapMaskValues(dataMat(:,1),trkIdx,[trkCount:trkCount+length(trkIdx)-1]);
        trkCount=trkCount+length(trkIdx);

        % assign matrix to cell array
        dataByProject{iProj,1}=dataMat;

        projNum{projCount,1}=iProj;
        projNum{projCount,2}=formatPath(temp.projData.anDir);

        projCount=projCount+1;
    end

    % concat all the data 
    allData=cell2mat(dataByProject);
    [temp.projData,M]=plusTipDynamParam(allData,temp.projData,1,0); % keep this on 1
    % and do not attempt to remove fields because this will give an error 

    if doBtw==1
        % put data in cell array for bwt group box plot
        allDataCell{iGroup,1}=allData;

        % make structure containing the concatenated distributions
        groupData(iGroup,1).info.name=btwGrpNames{iGroup,1};
        groupData(iGroup,1).info.groupListIdx=tempIdx;
        groupData(iGroup,1).info.stats= temp.projData.stats;
        groupData(iGroup,1).gs=M(~isnan(M(:,1)),1);
        groupData(iGroup,1).fs=M(~isnan(M(:,2)),2);
        groupData(iGroup,1).bs=M(~isnan(M(:,3)),3);
        groupData(iGroup,1).gl=M(~isnan(M(:,4)),4);
        groupData(iGroup,1).fl=M(~isnan(M(:,5)),5);
        groupData(iGroup,1).bl=M(~isnan(M(:,6)),6);
        groupData(iGroup,1).gd=M(~isnan(M(:,7)),7);
        groupData(iGroup,1).fd=M(~isnan(M(:,8)),8);
        groupData(iGroup,1).bd=M(~isnan(M(:,9)),9);
    end

    if doWtn==1
        tempDir=[saveDir filesep 'withinGroupComparisons' filesep btwGrpNames{iGroup,1}];
        if isdir(tempDir)
            rmdir(tempDir,'s')
        end
        mkdir(tempDir);

        % write out speed/lifetime/displacement distributions into a text file
        dlmwrite([tempDir filesep 'gs_fs_bs_gl_fl_bl_gd_fd_bd_' btwGrpNames{iGroup,1} '.txt'], M, 'precision', 3,'delimiter', '\t','newline', 'pc');

        if doPlot==1
            % save histograms of pooled distributions from iGroup
            plusTipMakeHistograms(M,tempDir);

            % here are the names for each movie in iGroup
            wtnGrpNames=repmat(btwGrpNames(iGroup,1),[length(dataByProject) 1]);
            for iName=1:length(wtnGrpNames)
                wtnGrpNames{iName}=[wtnGrpNames{iName} '_' num2str(iName)];
            end

            % make within-group boxplots (show each movie in iGroup)
            plusTipMakeBoxplots(dataByProject,wtnGrpNames,tempDir);
        end

    end

    clear dataByProject allData

end
if doBtw==1
    
    tempDir=[saveDir filesep 'btwGroupComparisons'];
    if isdir(tempDir)
        rmdir(tempDir,'s')
    end
    mkdir(tempDir);

    if doPlot==1
        % make between-group boxplots (show pooled data)
        plusTipMakeBoxplots(allDataCell,btwGrpNames,tempDir);
    end

%     % save movie reference list
%     groupData.projIdx=cell(length(projGroupName),3);
%     groupData.projIdx(:,1)=projGroupName;
%     groupData.projIdx(:,2:3)=projNum(:,1:2);

    save([tempDir filesep 'groupData'],'groupData');
end


cd(homeDir)