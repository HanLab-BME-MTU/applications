function [plusTipInfo,plusTipData]=plusTipGetGroupDataset(saveResult,useSavedGrp)
% plusTipGetGroupDataset makes dataset array with data from groupList

% INPUT
% saveResult:   1 to save the result, 0 to not save
% useSavedGrp:  1 to pick a previously-created groupList file, 0 to pick
%               groups on the fly
%
% OUTPUT
% plusTipInfo:  nProjects x 11 dataset array containing labels and
%               directory info that can be used to group data for statistics
% plusTipData:  nProjects x nStats dataset array containing statistics from
%               projData, created during post-processing
%
% Kathryn Applegate, 09/2009



% get output directory
if nargin<1 || isempty(saveResult) || saveResult~=1
    saveResult=0;
    saveDir=pwd;
else
    saveDir = uigetdir(pwd,'Select output directory for plusTipData');
end


if nargin<2 || isempty(useSavedGrp)
    % ask user to pick groups
    [projGroupDir,projGroupName]=plusTipPickGroups(1);
else
    homeDir=pwd;
    cd([saveDir filesep '..'])
    [fileName,pathName] = uigetfile('*.mat','Select saved groups file');
    if fileName==0
        return
    end
    load([pathName filesep fileName]);
    cd(homeDir)
end

% load first project to get fieldnames
temp=load([projGroupDir{1} filesep 'meta' filesep 'projData']);
projData=temp.projData;

% pick which data to extract
statNames1={'numTracks';'pair2pairDiffMicPerMinStd';'meanDisp2medianNNDistRatio';'percentFgapsReclass'};
statNames2=fieldnames(projData.stats);

% number of movies in the set
nMovs = size(projGroupName,1);

% add some extra rows for extra fields
nCols=length([statNames1; statNames2])+50;
statNamesInfo=cell(1,nCols);
statNamesData=cell(1,nCols);
plusTipInfo=cell(nMovs,nCols);
plusTipData=cell(nMovs,nCols);

for iMov=1:nMovs
    temp=load([projGroupDir{iMov} filesep 'meta' filesep 'projData']);
    projData=temp.projData;
    
    % GET PROJECT INFO
    plusTipInfo{iMov,1}=projGroupName{iMov}; % group label
    statNamesInfo{1}='label1';
    plusTipData{iMov,2}=''; % secondary label
    statNamesInfo{2}='label2';

    currentROI=formatPath(projData.anDir);
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
    plusTipInfo{iMov,11}=projData.anDir;
    
    
    % GET PROJECT DATA
    count=1;
    % add data pulled from projData
    for iName=1:length(statNames1)
        plusTipData{iMov,count}=projData.(statNames1{iName});
        statNamesData{count}=statNames1{iName};
        count=count+1;
    end

    % add data pulled from projData.stats
    for iName=1:length(statNames2)
        values=projData.stats.(statNames2{iName});
        tempName=statNames2{iName};
        % some measurements have more than one value - here we put each in
        % a separate column and label with 2,3,...
        for v=1:length(values)
            plusTipData{iMov,count}=values(v);
            if v==1
                statNamesData{count}=tempName;
            else
                statNamesData{count}=[tempName '_' num2str(v)];
            end
            count=count+1;
        end
    end
    
end

% get rid of extra rows
emptyIdxInfo=find(cell2mat(cellfun(@(x) sum(isempty(x)),statNamesInfo,'uniformoutput',0)),1,'first');
statNamesInfo(:,emptyIdxInfo:end)=[];
plusTipInfo(:,emptyIdxInfo:end)=[];

emptyIdxData=find(cell2mat(cellfun(@(x) sum(isempty(x)),statNamesData,'uniformoutput',0)),1,'first');
statNamesData(:,emptyIdxData:end)=[];
plusTipData(:,emptyIdxData:end)=[];

% contruct string to contain command for dataset construction of INFO
temp=plusTipInfo; 
clear plusTipInfo

NameObs = strcat({'Project '},num2str((1:size(temp,1))','%d'));
str='dataset(';
for i=1:length(statNamesInfo)
    s=['{temp(:,' num2str(i) '),statNamesInfo{' num2str(i) '}},'];
    str=[str s];
end
str=[str(1:end-1) ',''ObsNames'',NameObs);'];

plusTipInfo=eval(str);

% contruct string to contain command for dataset construction of DATA
temp=plusTipData;
clear plusTipData

NameObs = strcat({'Project '},num2str((1:size(temp,1))','%d'));
str='dataset(';
for i=1:length(statNamesData)
    s=['{cell2mat(temp(:,' num2str(i) ')),statNamesData{' num2str(i) '}},'];
    str=[str s];
end
str=[str(1:end-1) ',''ObsNames'',NameObs);'];

plusTipData=eval(str);


%a=grpstats([dataset({plusTipInfo.label1,'label1'})
%plusTipData],'label1',{'mean','sem'})