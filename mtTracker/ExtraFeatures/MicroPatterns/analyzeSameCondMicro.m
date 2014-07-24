function [ output_args ] = analyzeSameCondMicro(statDir,doWtn,doPlot,removeBeginEnd,HPattern,numWindows,windowSize)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin<1 || isempty(statDir)
    statDir=uigetdir(pwd,'Select Directory Where Statistics Can Be Found');
end



groupListDir = [statDir filesep 'groupLists']; 


% Make Directories for Analysis

analSaveDir1 = [ statDir filesep 'CompareBtwRegionTypes'];
if ~isdir(analSaveDir1); mkdir(analSaveDir1); end 

analSaveDir2 = [statDir filesep 'CompareWindowsWithinRegionType'];
if ~isdir(analSaveDir2); mkdir(analSaveDir2); end 




if HPattern== 1
    regionTypes{1,1} = 'Ad';
    regionTypes{2,1} = 'AdCorn';
    regionTypes{3,1} = 'NonAd';
    
else 
    regionTypes{1,1} = 'Ad';
    regionTypes{2,1} = 'NonAd';
end

% Pool the data from the groupList

for iWindow = 1:numWindows

    saveDir = [analSaveDir1 filesep num2str(iWindow*windowSize) '_uM'];
    if ~isdir(saveDir);  mkdir(saveDir); end
    
    
    
    %Load GroupList
    groupList = load([groupListDir filesep 'groupListCompareBtwRegionTypes_' num2str(iWindow*windowSize) 'uM.mat']);
    
    groupList = groupList.groupList;
    %Perform Pooling of Data
    groupData = plusTipExtractGroupData(groupList,removeBeginEnd); 
    
  plusTipPoolGroupData(groupData,saveDir, doWtn, doPlot); 
    
    
    % Collect Data
   
   
    
    [statsCellGS,statsCellFG, statsCellBG] = plusTipGetStats(saveDir,'stats',groupData,[],1,1,1,0);
    
  
    
    
    
    
    
    meanValueGS = statsCellGS(2:(length(regionTypes)+1),4);
    meanValueGL = statsCellGS(2:(length(regionTypes)+1),6);
    meanValueGD = statsCellGS(2:(length(regionTypes)+1),8);
    
    if iWindow == 1
       forPlot.meanValueGS = [regionTypes meanValueGS];
       forPlot.meanValueGL = [regionTypes meanValueGL];
       forPlot.meanValueGD = [regionTypes meanValueGD]; 
       
    else
        
       forPlot.meanValueGS = [forPlot.meanValueGS meanValueGS];
       forPlot.meanValueGL = [forPlot.meanValueGL meanValueGL]; 
       forPlot.meanValueGD = [forPlot.meanValueGD meanValueGD]; 
    
    end
    
    
    plusTipTestDistrib(groupData,saveDir,0.005,20,21); 
    
    
    
    

end % end for iWindows

saveDir = [analSaveDir1 filesep 'GreaterThan_' num2str(numWindows*windowSize) 'uM' ];
if ~isdir(saveDir);
    mkdir(saveDir) ;end ;

 %Load GroupList
    groupList = load([groupListDir filesep 'groupListCompareBtwRegionTypes_GreaterThan' num2str(numWindows*windowSize) 'uM.mat']);
    
    groupList = groupList.groupList;
    %Perform Pooling of Data
    groupData = plusTipExtractGroupData(groupList,removeBeginEnd ); 
    
    plusTipPoolGroupData(groupData, saveDir, doWtn, doPlot); 
    
   
    
    % Collect Data
   
    
    
    [ statsCellGS, statsCellFG, statsCellBG ] = plusTipGetStats(saveDir,'stats',groupData,[],1,1,1,0);
    
  plusTipTestDistrib(groupData,saveDir,0.005,20,21); 
    
    
    meanValueGS = statsCellGS(2:(length(regionTypes)+1),4);
    meanValueGL = statsCellGS(2:(length(regionTypes)+1),6);
    meanValueGD = statsCellGS(2:(length(regionTypes)+1),8);
    
       forPlot.meanValueGS = [forPlot.meanValueGS meanValueGS];
       forPlot.meanValueGL = [forPlot.meanValueGL meanValueGL]; 
       forPlot.meanValueGD = [forPlot.meanValueGD meanValueGD]; 
    
    
    
    save([statDir filesep 'forPlotting'],'forPlot');
    
    
    
    
    
    
%Do Same For Spatial Comparisons Within Region Type
    
    % Pool the data from the groupList


for iRegionType= 1:length(regionTypes)
    
    % Collect Data
   
    
    saveDir = [analSaveDir2 filesep char(regionTypes{iRegionType})];
    if ~isdir(saveDir) ; 
    mkdir(saveDir); end
    
    


    %Load GroupList
    groupList = load([groupListDir filesep 'groupListCompareWindows_' char(regionTypes{iRegionType,1}) '.mat']);
    groupList = groupList.groupList;
    
    
    %Perform Pooling of Data
    groupData = plusTipExtractGroupData(groupList,removeBeginEnd); 
    
    plusTipPoolGroupData(groupData, saveDir, doWtn, doPlot); 
    
    
    [ statsCellGS, statsCellFG, statsCellBG] = plusTipGetStats(saveDir,'stats',groupData,[],1,1,1,0);
    
    
    %Perform Statistical Tests
    
   plusTipTestDistrib(groupData,saveDir,0.005,20,21); 
   
   
    
    
end % end for iRegionType

        
     


    
end % 


     







