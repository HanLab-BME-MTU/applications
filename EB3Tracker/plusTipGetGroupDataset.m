function [plusTipDataset]=plusTipGetGroupDataset
% plusTipGetGroupDataset makes dataset array with dynamics parameters from movies in groupList

% INPUT
% user is asked for an output directory where dataset should be stored, as
% well as the projList file(s) from which data should be pulled
%
% OUTPUT
% plusTipDataset: dataset array containing directory info, tracking
%                 parameters used, and extracted MT dynamics parameters for
%                 all projects in projList.  if a project hasn't been
%                 processed, its row will contain NaNs and the "goodMovie"
%                 column will have a 0.
%
% Kathryn Applegate, 09/2009



homeDir=pwd;
saveDir=uigetdir(pwd,'Select output directory for plusTipDataset.');
cd(saveDir)

% ask user to select projList file and check which movies have been tracked
[allProjects,notDone]=plusTipCheckIfDone;

% keep only the ones that have been tracked
%allProjects(notDone,:)=[];

nProj=size(allProjects,1); % total number of columns
numDir=8; % number of columns devoted to proj path parts

for iProj=1:nProj
    try
        currentROI=formatPath(allProjects{iProj,1});
        

        if iProj==1
            % load first project to get fieldnames
            temp=load([currentROI filesep 'meta' filesep 'projData']);
            projData=temp.projData;
            
            % pick which data to extract
            dataNames={'numTracks';'pair2pairDiffMicPerMinStd';'meanDisp2medianNNDistRatio';'percentFgapsReclass'};
            statNames=fieldnames(projData.stats);
            trackParamNames=fieldnames(projData.trackingParameters);

            % the first 5+numDir variables will contain the project path and parts of
            % the path which can be used later to make grouping variables
            nVar=5+numDir+length(trackParamNames)+length(dataNames)+length(cell2mat(struct2cell(projData.stats)'));

            plusTipDataset=cell(nProj,nVar);
            varNames=cell(1,nVar);
        end

        % record project path, movie status placeholder (assumes all good for
        % now), and detection/tracking/postprocessing timestamps
        varNames{1}='projPath';
        plusTipDataset{iProj,1}=currentROI;

        varNames{2}='goodMovie';
        plusTipDataset{iProj,2}=1;

        % check whether a value exists for the timestamps
        varNames(3:5)={'detectTime','trackTime','postTime'};
        if ~isempty(allProjects{iProj,2})
            plusTipDataset{iProj,3}=allProjects{iProj,2}{1,1};
        else
            plusTipDataset{iProj,3}=[];
        end
        if ~isempty(allProjects{iProj,3})
            plusTipDataset{iProj,4}=allProjects{iProj,3}{1,1};
        else
            plusTipDataset{iProj,4}=[];
        end
        if ~isempty(allProjects{iProj,4})
            plusTipDataset{iProj,5}=allProjects{iProj,4}{1,1};
        else
            plusTipDataset{iProj,5}=[];
        end
        c=6; % next one (c is counter)
        
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
            words{roiIdx+1}=' ';
        end

        % invert order
        words=words(end:-1:1);
        for i=1:numDir
            varNames{c}=['dir' num2str(i)];
            if i>length(words)
                plusTipDataset{iProj,c}='';
            else
                plusTipDataset{iProj,c}=words{i};
            end
            c=c+1;
        end

        % load first project to get fieldnames
        temp=load([currentROI filesep 'meta' filesep 'projData']);
        projData=temp.projData;

        % add tracking parameters pulled from projData
        for iName=1:length(trackParamNames)
            varNames{c}=trackParamNames{iName};
            plusTipDataset{iProj,c}=projData.trackingParameters.(trackParamNames{iName});
            c=c+1;
        end

        % add data pulled from projData
        for iName=1:length(dataNames)
            varNames{c}=dataNames{iName};
            plusTipDataset{iProj,c}=projData.(dataNames{iName});
            c=c+1;
        end

        % add data pulled from projData.stats
        for iName=1:length(statNames)
            values=projData.stats.(statNames{iName});
            tempName=statNames{iName};
            % some measurements have more than one value (SEs) - here we put
            % each in a separate column and label with 2,3,...
            for v=1:length(values)
                if v==1
                    varNames{c}=tempName;
                else
                    varNames{c}=[tempName '_' num2str(v)];
                end
                plusTipDataset{iProj,c}=values(v);
                c=c+1;
            end
        end
    catch
        if iProj==1
            error('Problem with first project in projList - tracking or post-processing may be out of date')
        else
            % change goodMovie to 0 and all dynamics parameters to nan
            plusTipDataset{iProj,2}=0;
            for iVar=numDir+5+1:nVar
                plusTipDataset{iProj,iVar}=nan;
            end
        end
    end

end

% contruct string to contain command for dataset construction of dataset
% array
temp=plusTipDataset;
clear plusTipDataset

NameObs = strcat({'Project '},num2str((1:nProj)','%d'));
str='dataset(';
for i=1:nVar
    if i<=numDir+5 && i~=2
        s=['{temp(:,' num2str(i) '),varNames{' num2str(i) '}},'];
    else
        s=['{cell2mat(temp(:,' num2str(i) ')),varNames{' num2str(i) '}},'];
    end
    str=[str s];
end
str=[str(1:end-1) ',''ObsNames'',NameObs);'];

plusTipDataset=eval(str);
plusTipDataset.Properties.DimNames{1,1}='Projects';

% save result
fileName='plusTipDataset';
save([saveDir filesep fileName],'plusTipDataset')
if ispc
    extstr='.xls';
else
    extstr='.txt';
end
exportDatasetArray(plusTipDataset,'file',[saveDir filesep fileName extstr])


if any(plusTipDataset.goodMovie==0)
    msgbox('Check output: one ore more projects likely need to be retracked','Potential Problem')
end

cd(homeDir)