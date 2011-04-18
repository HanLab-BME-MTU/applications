function [  statsCellGS, statsCellFG, statsCellBG,stats] = plusTipGetStats(saveDir,statFileName,groupData,params2Extract,gsStats,fgStats,bgStats,percentStats)
%Extracts Data from groupDataSet for plotting in Excel etc.  
%
%INPUT:
% statFileName: Name of new file to be saved locally
%
% groupData: Output structure from plusTipPoolGroupData
% 
% params2Extract: Cell Array of parameters the user would like to extract
%                   if like to pick and choose different parameters
%                
%
% gsStats: 1 if you want to extract all growth stats to an excel file
% 
% fgStats: 1 if you want to extract all forward gap related stats
%
% bgStats: 1 if you want to extract all back gap related stats 
% 
% percentStats: 1 if you want to extract percent related stats

if nargin<1 || isempty(saveDir)
    saveDir=uigetdir(pwd,'Select Output Directory.');
end





% Extract input data (need in this format to get group names)
temp=struct2cell(groupData);
 
%Get number of groups in Group Data Set
[nGroups x] = size(groupData);

% Initiate Output
grpNames=cell(nGroups,1);

stats = zeros(nGroups,length(params2Extract));

%Extract Group Names
for i=1:nGroups
    grpNames{i,1}=temp{1,i}.name;
end

% If user included a custom cell with parameters they want to extract 
% Collect this data and write the output
if ~isempty(params2Extract);
    statNames = cell(1,length(params2Extract)+1);
    for iParam = 1:length(params2Extract)
        param = params2Extract{iParam};
        statNames{iParam+1} = param;
        for iGroup = 1:nGroups        
            stats(iGroup,iParam) = groupData(iGroup).info.stats.(param)(1);
       
        end     
    end
    stats1 = num2cell(stats);
    statsCellCus = [grpNames stats1];
    statsCellCus = [statNames; statsCellCus]; 
    filename = [statFileName,'_SelectedStats'];
    save([saveDir filesep filename],'stats_cellCus');
end

%
if gsStats == 1
    statNames{1,1} = [];   
    statNames{1,2} = 'nGrowths';
    statNames{1,3} = 'growth_speed_median';
    statNames{1,4} = 'growth_speed_mean_std';
    statNames{1,5} = 'growth_lifetime_median';
    statNames{1,6} = 'growth_lifetime_mean_std';
    statNames{1,7} = 'growth_length_median';
    statNames{1,8} = 'growth_length_mean_std';
    statNames{1,9} = 'percentTimeGrowth';
    statNames{1,10} = 'percentGrowthTerminal';
    
    for iParam = 1:9
        param = statNames{1,iParam+1};
        for iGroup = 1:nGroups        
            stats(iGroup,iParam) = groupData(iGroup).info.stats.(param)(1);
       
        end     
    end
    stats2 = num2cell(stats);
    statsCellGS = [grpNames stats2];
    statsCellGS = [statNames; statsCellGS]; 
    filename = [statFileName,'_GrowthStats'];
    save([saveDir filesep filename],'statsCellGS');
end     

if fgStats == 1 
  
    statNamesFG{1,1} = [];   
    statNamesFG{1,2} = 'nFgaps';
    statNamesFG{1,3} = 'fgap_speed_median';
    statNamesFG{1,4} = 'fgap_speed_mean_std';
    statNamesFG{1,5} = 'fgap_lifetime_median';
    statNamesFG{1,6} = 'fgap_lifetime_mean_std';
    statNamesFG{1,7} = 'fgap_length_median';
    statNamesFG{1,8} = 'fgap_length_mean_std';
    %statNamesFG{1,9} = 'fgap_freq_time';
    %statNamesFG{1,9} = 'fgap_freq_length';
    statNamesFG{1,9} = 'percentTimeFgap';
    statNamesFG{1,10} = 'avgIndivPercentTimeFgap';
    statNamesFG{1,11} = 'percentGrowthLinkedForward';
    
    for iParam = 1:10
        param = statNamesFG{1,iParam+1};
        for iGroup = 1:nGroups        
            stats(iGroup,iParam) = groupData(iGroup).info.stats.(param)(1);
       
        end     
    end
    stats3 = num2cell(stats);
    statsCellFG = [grpNames stats3];
    statsCellFG = [statNamesFG; statsCellFG]; 
    filename = [statFileName,'_FGapStats'];
    save([saveDir filesep filename],'statsCellFG');
end    

if bgStats ==1 
  
    statNames{1,1} = [];   
    statNames{1,2} = 'nBgaps';
    statNames{1,3} = 'bgap_speed_median';
    statNames{1,4} = 'bgap_speed_mean_std';
    statNames{1,5} = 'bgap_lifetime_median';
    statNames{1,6} = 'bgap_lifetime_mean_std';
    statNames{1,7} = 'bgap_length_median';
    statNames{1,8} = 'bgap_length_mean_std';
    %statNames{1,9} = 'bgap_freq_time';
    %statNames{1,9} = 'bgap_freq_length';
    statNames{1,9} = 'percentTimeBgap'; 
    statNames{1,10} = 'avgIndivPercentTimeBgap';
    statNames{1,11} = 'percentGrowthLinkedBackward';
    
    for iParam = 1:10
        param = statNames{1,iParam+1};
        for iGroup = 1:nGroups        
            stats(iGroup,iParam) = groupData(iGroup).info.stats.(param)(1);
       
        end     
    end
    stats = num2cell(stats);
    statsCellBG = [grpNames stats];
    statsCellBG = [statNames; statsCellBG]; 
    filename = [statFileName,'_BGapStats'];
    save([saveDir filesep filename],'statsCellBG');
end
    
if percentStats == 1
    statNamesPer{1,1} = [];
    statNamesPer{1,2} = 'percentGrowthTerminal';
    statNamesPer{1,3} = 'percentGrowthLinkedForward';
    statNamesPer{1,4} = 'percentGrowthLinkedBackward';
    statNamesPer{1,5} = 'percentTimeGrowth';
    statNamesPer{1,5} = 'percentTimeFgap';
    statNamesPer{1,6} = 'percentTimeBgap';

        for iParam = 1:5
        param = statNamesPer{1,iParam+1};
            for iGroup = 1:nGroups        
                stats(iGroup,iParam) = groupData(iGroup).info.stats.(param)(1);
       
            end     
        end
    stats = num2cell(stats);
    statsCellPer = [grpNames stats];
    statsCellPer = [statNamesPer; statsCellPer]; 
    filename = [statFileName,'_Percentages'];
    save([saveDir filesep filename],'statsCellPer');
end
        




end

