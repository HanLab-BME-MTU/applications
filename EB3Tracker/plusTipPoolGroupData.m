function plusTipPoolGroupData(groupList)


% get output directory and load groupLIst
saveDir=uigetdir(pwd,'Select output directory for pooled group data.');


projGroupName=groupList(:,1);
projGroupDir=groupList(:,2);


% fix the names if there are spaces or hyphens and append prefix 'grp'
projGroupName=cellfun(@(x) strrep(x,'-','_'),projGroupName,'uniformoutput',0);
projGroupName=cellfun(@(x) strrep(x,' ','_'),projGroupName,'uniformoutput',0);
projGroupName=cellfun(@(x) ['grp_' x],projGroupName,'uniformoutput',0);


% count unique groups and keep them in order of the original list
[grpNames,m,projGroupIdx] = unique(projGroupName);
[b,idx]=sort(m);
grpNames=grpNames(idx);


projCount=1; % all-project counter
projNum=cell(length(projGroupName),1); % within-group counter
allDataCell=cell(length(grpNames),1); % cell-array for data matrix
for iGroup = 1:length(grpNames)

    % indices of projects in iGroup
    tempIdx=strmatch(grpNames(iGroup),projGroupName,'exact');

    data=cell(length(tempIdx),1);
    trkCount=1;
    for iProj = 1:length(tempIdx)

        temp = load([projGroupDir{tempIdx(iProj)} filesep 'meta' filesep 'projData']);
        [dummy1,dummy2,dataMatCrpSecMic]=plusTipMergeSubtracks(temp.projData);

        % reassign the track numbers so when combined from multiple projects they don't repeat
        trkIdx=unique(dataMatCrpSecMic(:,1));
        dataMatCrpSecMic(:,1)=swapMaskValues(dataMatCrpSecMic(:,1),trkIdx,[trkCount:trkCount+length(trkIdx)-1]);
        trkCount=trkCount+length(trkIdx);
        
        % assign matrix to cell array
        data{iProj,1}=dataMatCrpSecMic;

        projNum{projCount,1}=iProj;
        projNum{projCount,2}=formatPath(temp.projData.anDir);

        projCount=projCount+1;
    end

    % concat all the data
    allData=cell2mat(data);
    allDataCell{iGroup,1}=allData;

    % get stats for the whole group
    [groupData.(grpNames{iGroup,1}).stats,M]=plusTipDynamParam(allData);

    % store speed/lifetime/displacement lists for growth/fgap/bgap
    groupData.(grpNames{iGroup,1}).gs=M(~isnan(M(:,1)),1);
    groupData.(grpNames{iGroup,1}).fs=M(~isnan(M(:,2)),2);
    groupData.(grpNames{iGroup,1}).bs=M(~isnan(M(:,3)),3);
    groupData.(grpNames{iGroup,1}).gl=M(~isnan(M(:,4)),4);
    groupData.(grpNames{iGroup,1}).fl=M(~isnan(M(:,5)),5);
    groupData.(grpNames{iGroup,1}).bl=M(~isnan(M(:,6)),6);
    groupData.(grpNames{iGroup,1}).gd=M(~isnan(M(:,7)),7);
    groupData.(grpNames{iGroup,1}).fd=M(~isnan(M(:,8)),8);
    groupData.(grpNames{iGroup,1}).bd=M(~isnan(M(:,9)),9);
   
    
    tempDir=[saveDir filesep 'withinGroupComparisons' filesep grpNames{iGroup,1}];
    if ~isdir(tempDir)
        mkdir(tempDir);
    end
    plusTipMakeHistograms(M,tempDir);

    % write out speed/lifetime/displacement distributions into a text file
    dlmwrite([tempDir filesep 'gs_fs_bs_gl_fl_bl_gd_fd_bd_' grpNames{iGroup,1} '.txt'], M, 'precision', 3,'delimiter', '\t','newline', 'pc');

    % make within-group boxplots
    setName=repmat(grpNames(iGroup,1),[length(data) 1]);
    for iName=1:length(setName)
        setName{iName}=[setName{iName} '_' num2str(iName)];
    end
    plusTipMakeBoxplots(data,setName,tempDir);


    clear data allData

end

% make between-group boxplots
tempDir=[saveDir filesep 'btwGroupComparisons'];
if ~isdir(tempDir)
    mkdir(tempDir);
end
plusTipMakeBoxplots(allDataCell,grpNames,tempDir);

% save movie reference list
groupData.projIdx=cell(length(projGroupName),3);
groupData.projIdx(:,1)=projGroupName;
groupData.projIdx(:,2:3)=projNum(:,1:2);

save([tempDir filesep 'groupData'],'groupData');
