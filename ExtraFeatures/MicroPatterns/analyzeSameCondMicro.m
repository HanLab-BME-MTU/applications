function [ output_args ] = analyzeSameCondMicro(statDir,doBtw,doWtn,doPlot,removeBeginEnd,HPattern,numWindows,windowSize)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin<1 || isempty(statDir)
    statDir=uigetdir(pwd,'Select Directory Where Statistics Can Be Found');
end



groupListDir = [statDir filesep 'groupLists']; 


% Make Directories for Analysis

analSaveDir1 = [ statDir filesep 'CompareBtwRegionTypes'];
mkdir(analSaveDir1);

analSaveDir2 = [statDir filesep 'CompareWindowsWithinRegionType'];
mkdir(analSaveDir2);



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
    mkdir(saveDir);
    
    
    
    %Load GroupList
    groupList = load([groupListDir filesep 'groupListCompareBtwRegionTypes_' num2str(iWindow*windowSize) 'uM.mat']);
    
    groupList = groupList.groupList;
    %Perform Pooling of Data
    [groupData] = plusTipPoolGroupData(groupList,saveDir, doBtw, doWtn, doPlot, removeBeginEnd); 
    
    
    % Collect Data
   
   
    
    [statsCellGS,statsCellFG, statsCellBG] = plusTipGetStats(saveDir,'stats',groupData,[],1,1,1,0);
    
    
    %Perform Statistical Tests
    [nRows nCols] = size(statsCellBG);
    
    zeroBGaps = 0;
    
    for i = 1:(nRows-1)
          if cell2mat(statsCellBG(1+i,2)) == 0 
              zeroBGaps = 1 ;
          else 
          end 
    end 
    
  zeroFGaps = 0 ;
    for i = 1: nRows-1
        if cell2mat(statsCellFG(1+i,2)) == 0 
            zeroFGaps = 1; 
        else 
        end 
    end 
    
    if zeroBGaps == 1  && zeroFGaps == 0
            params{1,1} = 'gs';
            params{1,2} = 'fs';
            params{1,3} = 'gl';
            params{1,4} = 'fl';
            params{1,5} = 'gd';
            params{1,6} = 'fd';
        
            [discrimMat] = plusTipTestDistrib(saveDir,groupData,params);
    
        elseif zeroFGaps ==1 && zeroBGaps == 0
        
            params{1,1} = 'gs';
            params{1,2} = 'bs';
            params{1,3} = 'gl';
            params{1,4} = 'bl';
            params{1,5} = 'gd';
            params{1,6} = 'bd';
        
            [discrimMat] = plusTipTestDistrib(saveDir,groupData,params);
        
        elseif zeroFGaps == 1 && zeroBGaps == 1 
        
            params{1,1} = 'gs';
            params{1,2} = 'gl';
            params{1,3} = 'gd';
        
        else 
            [discrimMat] = plusTipTestDistrib(saveDir,groupData,[]);
    end 
    
    
    
    
    
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
    
    
    
    
    
    
    

end % end for iWindows

saveDir = [analSaveDir1 filesep 'GreaterThan_' num2str(numWindows*windowSize) 'uM' ];
mkdir(saveDir);

 %Load GroupList
    groupList = load([groupListDir filesep 'groupListCompareBtwRegionTypes_GreaterThan' num2str(numWindows*windowSize) 'uM.mat']);
    
    groupList = groupList.groupList;
    %Perform Pooling of Data
    [groupData] = plusTipPoolGroupData(groupList, saveDir, doBtw, doWtn, doPlot, removeBeginEnd); 
    
   
    
    % Collect Data
   
    
    
    [ statsCellGS, statsCellFG, statsCellBG, stats ] = plusTipGetStats(saveDir,'stats',groupData,[],1,1,1,0);
    
    
    %Perform Statistical Tests
    
    [nRows nCols] = size(statsCellBG);
    
    zeroBGaps = 0;
    
    for i = 1: nRows-1
          if cell2mat(statsCellBG(1+i,2)) == 0 
              zeroBGaps = 1 ;
          else 
          end 
    end 
    
   zeroFGaps = 0 ;
    for i = 1: nRows-1
        if cell2mat(statsCellFG(1+i,2)) == 0 
            zeroFGaps = 1; 
        else 
        end 
    end 
    
    if zeroBGaps == 1  && zeroFGaps == 0
            params{1,1} = 'gs';
            params{1,2} = 'fs';
            params{1,3} = 'gl';
            params{1,4} = 'fl';
            params{1,5} = 'gd';
            params{1,6} = 'fd';
        
            [discrimMat] = plusTipTestDistrib(saveDir,groupData,params);
    
        elseif zeroFGaps ==1 && zeroBGaps == 0
        
            params{1,1} = 'gs';
            params{1,2} = 'bs';
            params{1,3} = 'gl';
            params{1,4} = 'bl';
            params{1,5} = 'gd';
            params{1,6} = 'bd';
        
            [discrimMat] = plusTipTestDistrib(saveDir,groupData,params);
        
        elseif zeroFGaps == 1 && zeroBGaps == 1 
        
            params{1,1} = 'gs';
            params{1,2} = 'gl';
            params{1,3} = 'gd';
        
        else 
            [discrimMat] = plusTipTestDistrib(saveDir,groupData,[]);
    end 
    
    
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
    mkdir(saveDir);
    
    


    %Load GroupList
    groupList = load([groupListDir filesep 'groupListCompareWindows_' char(regionTypes{iRegionType,1}) '.mat']);
    groupList = groupList.groupList;
    
    
    %Perform Pooling of Data
    [groupData] = plusTipPoolGroupData(groupList, saveDir, doBtw, doWtn, doPlot, removeBeginEnd); 
    
    
    [ statsCellGS, statsCellFG, statsCellBG] = plusTipGetStats(saveDir,'stats',groupData,[],1,1,1,0);
    
    
    %Perform Statistical Tests
    
    [nRows nCols] = size(statsCellBG);
    
    zeroBGaps = 0;
    
    for i = 1: nRows-1
          if cell2mat(statsCellBG(1+i,2)) == 0 
              zeroBGaps = 1 ;
          else 
          end 
    end
    
    zeroFGaps = 0 ;
    for i = 1: nRows-1
        if cell2mat(statsCellFG(1+i,2)) == 0 
            zeroFGaps = 1; 
        else 
        end 
    end 
    
    if zeroBGaps == 1  && zeroFGaps == 0
            params{1,1} = 'gs';
            params{1,2} = 'fs';
            params{1,3} = 'gl';
            params{1,4} = 'fl';
            params{1,5} = 'gd';
            params{1,6} = 'fd';
        
            [discrimMat] = plusTipTestDistrib(saveDir,groupData,params);
    
        elseif zeroFGaps ==1 && zeroBGaps == 0
        
            params{1,1} = 'gs';
            params{1,2} = 'bs';
            params{1,3} = 'gl';
            params{1,4} = 'bl';
            params{1,5} = 'gd';
            params{1,6} = 'bd';
        
            [discrimMat] = plusTipTestDistrib(saveDir,groupData,params);
        
        elseif zeroFGaps == 1 && zeroBGaps == 1 
        
            params{1,1} = 'gs';
            params{1,2} = 'gl';
            params{1,3} = 'gd';
        
        else 
            [discrimMat] = plusTipTestDistrib(saveDir,groupData,[]);
    end 
    
    
    
end % end for iRegionType

        
     


    
end % 


     







