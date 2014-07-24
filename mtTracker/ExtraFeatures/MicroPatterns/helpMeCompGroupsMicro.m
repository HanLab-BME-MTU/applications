function [ output_args ] = helpMeCompGroupsMicro
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[groupList, saveDir] = combineGroupListFiles(1); 

[groupData] = plusTipPoolGroupData(groupList,saveDir,1,1,1,1); 

plusTipTestDistrib(saveDir,groupData,[]); 

plusTipGetStats(saveDir,'stats',groupData,[],1,1,1,0); 

end

