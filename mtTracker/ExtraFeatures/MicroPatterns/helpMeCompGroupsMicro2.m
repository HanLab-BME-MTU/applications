function [ output_args ] = helpMeCompGroupsMicro2( statDir1,statDir2,saveDir)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if nargin<1 || isempty(statDir1)
    statDir1=uigetdir(pwd,'Select Statistics Directory 1.');
end

if nargin<2 || isempty(statDir2) 
    statDir2= uigetdir(pwd,'Select Statistics Directory 2.'); 
end 

if nargin<3 || isempty(saveDir) 
    saveDir = uigetdir(pwd, 'Select SaveDir'); 
end 



groupLists = cell(9,1) ; 

groupLists{1,1} = 'groupListAd_3uM.mat';
groupLists{2,1} = 'groupListAd_6uM.mat';
groupLists{3,1} = 'groupListAd_GreaterThan_6uM.mat';
groupLists{4,1} = 'groupListAdCorn_3_uM.mat'; 
groupLists{5,1} = 'groupListAdCorn_6_uM.mat'; 
groupLists{6,1} = 'groupListAdCornGreaterThan_6uM.mat'; 
groupLists{7,1} = 'groupListNonAd3_uM.mat'; 
groupLists{8,1} = 'groupListNonAd6_uM.mat'; 
groupLists{9,1} = 'groupListNonAd_GreaterThan6uM.mat'; 

subRegions = cell(3,1);
subRegions{1} = 'Adhesion';
subRegions{2} = 'AdhesionCorner';
subRegions{3} = 'NonAdhesion'; 

values = cell(3,1); 
values{1} = '3uM'; 
values{2} = '6uM'; 
values{3} = 'Greater6uM';

if ~isdir([saveDir filesep 'SubRegions']) 
    mkdir([saveDir filesep 'SubRegions'])
end 

for i = 1:length(subRegions)
  for   j = 1:length(values) 
      dir = [saveDir filesep 'SubRegions' filesep subRegions{i} filesep values{j}]; 
    if ~isdir(dir)
    mkdir(dir)
    end 
  end 
end 
    

saveList = cell(9,1); 
saveList{1,1} = ['SubRegions' filesep 'Adhesion' filesep '3uM']; 
saveList{2,1} = ['SubRegions' filesep 'Adhesion' filesep '6uM']; 
saveList{3,1} = ['SubRegions' filesep 'Adhesion' filesep 'Greater6uM']; 
saveList{4,1} = ['SubRegions' filesep 'AdhesionCorner' filesep '3uM'];
saveList{5,1} = ['SubRegions' filesep 'AdhesionCorner' filesep '6uM']; 
saveList{6,1} = ['SubRegions' filesep  'AdhesionCorner' filesep 'Greater6uM']; 
saveList{7,1} = ['SubRegions' filesep 'NonAdhesion' filesep '3uM'];  
saveList{8,1} = ['SubRegions' filesep 'NonAdhesion' filesep '6uM']; 
saveList{9,1} = ['SubRegions' filesep 'NonAdhesion' filesep 'Greater6uM']; 



for i= 1:length(groupLists)
    group1 = load([statDir1 filesep 'groupLists' filesep groupLists{i}]); 
    group2 = load([statDir2 filesep 'groupLists' filesep groupLists{i}]); 
groupList = [group1.groupList;group2.groupList]; 



saveDirFinal = [saveDir filesep saveList{i}];

save([saveDirFinal filesep groupList{i}],'groupList'); 

groupData = plusTipExtractGroupData(groupList,1); 

plusTipPoolGroupData(groupData,saveDirFinal,0,1); 

plusTipTestDistrib(groupData,saveDirFinal,0.005,1,20); 

plusTipGetStats(saveDirFinal,'stats',groupData,[],1,1,1,0);
end 



end

