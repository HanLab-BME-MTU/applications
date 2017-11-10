function [ output_args ] = makeCombinedRosePlot( groupList, saveDir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checkInput
% get projData in correct format
if nargin<1 || isempty(groupList)
    % if not given as input, ask user for ROI directory
    % assume images directory is at same level
  [fileName,pathName]=uigetfile('*.mat','Please select groupList');
   if isequal(fileName,0)
        return
    end
   groupList=load([pathName filesep fileName]);
   % projData=projData.projData;
end


% get output directory from the user
if ~nargin<2 || isempty(saveDir)
    saveDir=uigetdir(pwd,'Please select OUTPUT directory');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Group Stuff
projGroupName=groupList(:,1);
projGroupDir=cellfun(@(x) formatPath(x),groupList(:,2),'uniformoutput',0);


% fix the names if there are spaces or hyphens and append prefix 'grp'
projGroupName=cellfun(@(x) strrep(x,'-','_'),projGroupName,'uniformoutput',0);
projGroupName=cellfun(@(x) strrep(x,' ','_'),projGroupName,'uniformoutput',0);
projGroupName=cellfun(@(x) ['grp_' x],projGroupName,'uniformoutput',0);

[btwGrpNames, m, projGroupIdx] = unique(groupList(:,1)); 
[b,idx]=sort(m);
btwGrpNames=btwGrpNames(idx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% body 
for iGroup = 1: length(btwGrpNames)
    
    tempList = groupList(projGroupIdx == iGroup,:); 
 
 allAngles = [];
  
  for iProj = 1:length(tempList)
      projData = load([formatPath(char(tempList(iProj,2))) filesep 'meta' filesep 'projData.mat']); 
      projData = projData.projData; 
      allAngles = [allAngles ; projData.anglesFrame2Frame_INSIDE_REGION] ; 
  end
  saveFig1 = figure;
  [tout rout] = rose(allAngles); 
  polar(tout,rout./max(rout)); 
  groupName = char(btwGrpNames(iGroup));
  groupNameTitle = regexprep(groupName,'_',' '); 
  title({groupNameTitle; 'Rose Plot of All Frame-to-Frame Growth Displacement (Converted to Theta)'});
  saveas(saveFig1,[saveDir filesep [groupName 'angles_histogram.eps']], 'psc2');
  

     

end

