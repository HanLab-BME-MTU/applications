%% Yukako's guide to amazing sub-region analysis!!!!


%% 
% run getProj to create a list of all the projects. you will have to choose
% a list of projects to use for batch sub-roi detection
getProj


%% HERE WE SET UP THE SUB-REGION SELECTION SETTINGS
% 1 if you want to exclude some tracks, like those coming out of the
% centrosome
excludeRegion=1;

% microns from the cell edge for automatic region selection
micFromEdge=10; 

% 1 if you want only tracks which spend half their lives or more in the
% sub-roi to count, 0 to include any with a min number of frames
midPoint=1; 

% 3 is default, but it doesn't get used if midPoint=1
minFrames=[];

% now run the function.  be prepared to draw a polygon around each cell and
% exclude a small region if you chose that option.  
plusTipBatchSubRoi(excludeRegion,micFromEdge,midPoint,minFrames);


%% MAKE LISTS OF PROJECTS IN CELL CENTER VS. PERIPHERY
% now run getProj to get lists of the sub-projects from the cell center
% (those named 'sub_1') and those from the cell periphery ('sub_2)
% [projList_cellCenter]=getProj('sub_1');
% [projList_cellPeriph]=getProj('sub_2');

projList_all=getProj;

%% MAKE GROUPLIST FILES NEEDED TO RUN BATCH QUAD PLOTS
% now run the group picking function.  pick the movies corresponding to
% each group for the whole cell, cell center, and cell periphery, and name
% them appropriately. 

%[groupList_Center]=plusTipPickGroups([],[],projList_cellCenter);
% now do the same for the cell periphery groups.
%[groupList_Periph]=plusTipPickGroups([],[],projList_cellPeriph);
% if you need to re-do them, that's ok - just delete the old ones and save
% the new ones.

% Now, pick the groups (avoid using a dash or space in the names).
% These should correspond to the groups:
%   roi_1 : whole cell
%   sub_1 : cell center
%   sub_2 : cell periphery
[groupList]=plusTipPickGroups([],[],projList_all);


% if you need to, you can combine multiple projList files
[groupList]=combineGroupListFiles(1);


% figure out good parameters from the whole data set
% grpData gives you a structure with [speed lifetime displacement] data
% from all the groups
[groupListSpeedLifeDisp,grpData]=plusTipGroupListSummary(groupList);

allSpeeds=grpData.allGroups(:,1);
allLifetimes=grpData.allGroups(:,2);

% these might be good values to use for the divisions - mean speed and
% lifetime from all groups in the whole dataset
meanSpeed=mean(allSpeeds);
meanLife=mean(allLifetimes);

figure; hist(allSpeeds,25)
figure; hist(allLifetimes,25)

%% MAKE THE BATCH QUAD PLOTS

% first set the input parameters
speedLims=[0 40]; % speed limits for x-axis
speedDiv=9.9; % speed division mark in microns/min
lifeLims=[0 40]; % lifetime limits for y-axis
lifeDiv=9; % lifetime division mark in seconds
remBegEnd=0; % 1 to remove data at beginning/end of the movie, 0 to keep it
timeRange=[1 30]; % frame range to use

% now run the quadrant scatter plot fuction
plusTipBatchQuadPlot(groupList,speedLims,speedDiv,lifeLims,lifeDiv,remBegEnd,timeRange,1);


% for getting out mean/std and frequencies from groupData (result of
% plusTipPoolGroupData
for i=1:length(groupData)
    nGrowths=groupData(i,1).info.stats.nGrowths;
    temp=groupData(i,1).info.stats.growth_speed_mean_SE;
    meanStd(i,1:2)=[temp(1) temp(2)*sqrt(nGrowths)]
    
    fgapFreq(i,1)=groupData(i,1).info.stats.fgap_freq_time_mean_SE(1);
    
end

popRYGB,nPrctRYGB]=plusTipGroupDataQuadPlot(groupData,'growthSpeed',speedDiv,...
    'growthLifetime',lifeDiv,speedLims,lifeLims);



