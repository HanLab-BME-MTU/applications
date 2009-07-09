function [movDataSet]=plusTipExtractGroupData(projGroupDir,projGroupName)
% plusTipExtractGroupData extracts various kinds of data from movies grouped together
% this is used by plusTipScatterplot

% INPUT : projGroupDir : cell array containing paths to chosen projects
%         projGroupName: cell array containing the group name for each
%                        project
% OUTPUT: movDataSet: dataset array containing groups labeled with
% % user-entered group name and the following parameters extracted for each:
%                     growthSpeedMean 
%                     growthSpeedStd
%                     Ppause 
%                     Pcat
%                     pauseSpeedMean
%                     pauseSpeedStd 
%                     shrinkSpeedMean 
%                     shrinkSpeedStd

movDataSet=[];

% collect the data
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

