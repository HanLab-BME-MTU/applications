function [ output_args ] = collectParametersForModeling(dir,micropatternOutput)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if isempty(dir)
    dir = uigetdir(pwd,'Select a Directory Where the Data is HPattern Data is Stored'); 
end 


% set up region types
regionTypes{1,1} = 'NonAdhesion';
regionTypes{2,1} = 'Adhesion';
regionTypes{3,1} = 'AdhesionCorner';

maskCurrent = micropatternOutput{1}.maskParams; % for now just do the first
numWindows = maskCurrent.numWindows; 

windows{1,1} = '3uM'; 
windows{2,1} = '6uM'; 
windows{3,1} = 'GreaterThan_6uM'; 

 % non directional folder 
  zoneCount = 1;
% go into each folder
for iRegion = 1:3
    for iWindow = 1:numWindows+1
      s= load([dir filesep 'nonDirectional' filesep 'SubRegions' filesep regionTypes{iRegion,1} filesep windows{iWindow} filesep 'groupData']);
        
    
        for iGroup = 1:numel(micropatternOutput{1}.groupData)
            
            values(1,zoneCount,iGroup) =  s.groupDataCurrent.pooledStats{iGroup}.growth_speed_mean_INSIDE_REGION;
            values(2,zoneCount,iGroup) = s.groupDataCurrent.pooledStats{iGroup}.bgap_speed_mean;
            values(3,zoneCount,iGroup) = 1/ s.groupDataCurrent.pooledStats{iGroup}.growth_lifetime_mean_INSIDE_REGION;
            values(4,zoneCount,iGroup) = 1/s.groupDataCurrent.pooledStats{iGroup}.bgap_lifetime_mean;
        
 
            
            %     params{iGroup}.growthVel = growthVel;
            %     params{iGroup}.bgapVel = bgapVel;
            %     params{iGroup}.catFreq = catFreq;
            %     params{iGroup}.rescueFreq = rescueFrequ;
            %     params{iGroup}.numNucl = numNucl;
            %     params{iGroup}.termEvents = termEvents;
            
            
            
        end
        zoneCount = zoneCount+1 ;
    end
end


zoneCount = 1; 

% SubTrackEnd Param
for iRegion = 1:3
    for iWindow = 1:numWindows+1
      s= load([dir filesep 'subTrackEnd' filesep 'SubRegions' filesep regionTypes{iRegion,1} filesep windows{iWindow} filesep 'groupData']);
        
    
        for iGroup = 1:numel(micropatternOutput{1}.groupData)
            
            if (iWindow ==1 && iRegion == 1)
                regions = numel(s.groupDataCurrent.stats{iGroup});
                numCells(iGroup) = regions./2;
            end
            
             %norm by region number 
            if (iWindow == 1 || iWindow ==2) 
                mult = 2; % two regions ad or non ad regions
            else 
                mult = 4; % 4 ending in ad corn regions
            end
            
            values(5,zoneCount,iGroup) =  s.groupDataCurrent.pooledStats{iGroup}.nGrowths/numCells(iGroup);
            values(6,zoneCount,iGroup) = s.groupDataCurrent.pooledStats{iGroup}.nGrowthTermEvents/numCells(iGroup); 
            
            
            
            subTrackDist = arrayfun(@(x) s.groupDataCurrent.stats{iGroup}{x}.percentWholeCellTracksTermInSubRegion,1:numel(s.groupDataCurrent.stats{iGroup}),'uniformOutput',0);
            
            values(7,zoneCount,iGroup) = mult*nanmean(cell2mat(subTrackDist)); 
            
        end 
          
            
 
            
            %     params{iGroup}.growthVel = growthVel;
            %     params{iGroup}.bgapVel = bgapVel;
            %     params{iGroup}.catFreq = catFreq;
            %     params{iGroup}.rescueFreq = rescueFrequ;
            %     params{iGroup}.numNucl = numNucl;
            %     params{iGroup}.termEvents = termEvents;
            
            
            zoneCount = zoneCount+1 ;
        end
        
end

zoneCount = 1; 
% Nucleation Param
for iRegion = 1:3
    for iWindow = 1:numWindows+1
      s= load([dir filesep 'nucleation' filesep 'SubRegions' filesep regionTypes{iRegion,1} filesep windows{iWindow} filesep 'groupData']);
        
    
        for iGroup = 1:numel(micropatternOutput{1}.groupData)
            
            values(8,zoneCount,iGroup) =  s.groupDataCurrent.pooledStats{iGroup}.numNucleationEvents/numCells(iGroup);
            
            nucDensityDist = arrayfun(@(x) s.groupDataCurrent.stats{iGroup}{x}.nucleationDensity,1:numel(s.groupDataCurrent.stats{iGroup}),'uniformOutput',0);
            
            values(9,zoneCount,iGroup) = nanmean(cell2mat(nucDensityDist)); 
            
            
            %norm by region number 
            if (iWindow == 1 || iWindow ==2) 
                mult = 2; % two regions ad or non ad regions
            else 
                mult = 4; % 4 ending in ad corn regions
            end
            ratioTracksNucDist = arrayfun(@(x) s.groupDataCurrent.stats{iGroup}{x}.percentWholeCellTracksNucInSubRegion,1:numel(s.groupDataCurrent.stats{iGroup}),'uniformOutput',0); 
            
            values(10,zoneCount,iGroup) = mult*nanmean(cell2mat(ratioTracksNucDist)); 
            
            
            
        end 
          
            
 
            
            %     params{iGroup}.growthVel = growthVel;
            %     params{iGroup}.bgapVel = bgapVel;
            %     params{iGroup}.catFreq = catFreq;
            %     params{iGroup}.rescueFreq = rescueFrequ;
            %     params{iGroup}.numNucl = numNucl;
            %     params{iGroup}.termEvents = termEvents;
            
            zoneCount = zoneCount+1 ;
            
        end
        
end

% get region names
count = 1; 
for iRegion = 1:3
    for iWindow = 1:numWindows+1
zone{count,1} = [regionTypes{iRegion} windows{iWindow}];  
count = count +1;
    end 
end

%  set variable names 
varNames{1,1} = 'Growth Velocity (um per min)'; 
varNames{2,1} = 'Bgap (Shrinkage) Velocity (um per min)'; 
varNames{3,1} = 'Catastrophe Freq (Inverse Sec)';
varNames{4,1} = 'Rescue Frequency (Inverse Sec) '; 
varNames{5,1} = 'Number Of Growth Subtracks Ending In Region (Number per Cell) '; 
varNames{6,1} = 'Number Of Microtubule Trajectories Terminally Ending In Region (ie Terminal Catastrophes) (Number per Cell)'; 
varNames{7,1} = 'Ratio of Microtubule Trajectories Ending in Region to Total Microtubule Trajectories Observed In Whole Cell ';
varNames{8,1} = 'Number Of Nucleation Events (Microtubule Trajectory Initiation Sites) In Region Per Cell (Number Per Cell)';
varNames{9,1} = 'Nucleation Density (Number of Nucleation Events per um^2 per Sec)';
varNames{10,1} = 'Ratio of Microtubule Trajectories in Region to Total Microtubule Trajectories Observed In Whole Cell'; 

% get groupNames
 s= load([dir filesep 'nonDirectional' filesep 'SubRegions' filesep regionTypes{1,1} filesep windows{1,1} filesep 'groupData']);
        
groupNames = s.groupDataCurrent.names;
idxSort = [6,9,3,5,8,2,4,7,1]'; 
zone = zone(idxSort); 
values = values(:,idxSort,:); 

for iGroup = 1:numel(micropatternOutput{1}.groupData)
  
groupName = strrep(groupNames{iGroup},'NonAdhesion_3uM','');    
modelParams  = dataset({values(:,:,iGroup),zone{:}},'ObsNames',varNames); 
export(modelParams,'file',['Model_Parameters', groupName]);  
end 


