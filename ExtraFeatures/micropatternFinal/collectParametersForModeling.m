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
        
    
        for iGroup = 1:2
            
            values(1,zoneCount,iGroup) =  s.groupDataCurrent.pooledStats{iGroup}.growth_speed_mean_INSIDE_REGION;
            values(2,zoneCount,iGroup) = s.groupDataCurrent.pooledStats{iGroup}.bgap_speed_mean;
            values(3,zoneCount,iGroup) = s.groupDataCurrent.pooledStats{iGroup}.growth_lifetime_mean_INSIDE_REGION;
            values(4,zoneCount,iGroup) = s.groupDataCurrent.pooledStats{iGroup}.bgap_lifetime_mean;
            
 
            
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
        
    
        for iGroup = 1:2
            
            values(5,zoneCount,iGroup) =  s.groupDataCurrent.pooledStats{iGroup}.nGrowths;
            values(6,zoneCount,iGroup) = s.groupDataCurrent.pooledStats{iGroup}.nGrowthTermEvents; 
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
        
    
        for iGroup = 1:2
            
            values(7,zoneCount,iGroup) =  s.groupDataCurrent.pooledStats{iGroup}.numNucleationEvents;
            
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
varNames{2,1} = 'Bgap Velocity (um per min)'; 
varNames{3,1} = 'Catastrophe Freq (Inverse Sec)';
varNames{4,1} = 'Rescue Frequency (Inverse Sec) '; 
varNames{5,1} = 'Number Of Growth Subtracks Ending In Region (Total Per Movie)'; 
varNames{6,1} = 'Number Of Compound Tracks Ending In Region (ie Terminal Catastrophes) (Total Per Movie)'; 
varNames{7,1} = 'Number Of Nucleation Events In Region (Total Per Movie) '; 

% get groupNames
 s= load([dir filesep 'nonDirectional' filesep 'SubRegions' filesep regionTypes{1,1} filesep windows{1,1} filesep 'groupData']);
        
groupNames = s.groupDataCurrent.names;
idxSort = [6,9,3,5,8,2,4,7,1]'; 
zone = zone(idxSort); 
values = values(:,idxSort,:); 

for iGroup = 1:2
  
groupName = strrep(groupNames{iGroup},'NonAdhesion_3uM','');    
modelParams  = dataset({values(:,:,iGroup),zone{:}},'ObsNames',varNames); 
export(modelParams,'file',['Model_Parameters', groupName]);  
end 


