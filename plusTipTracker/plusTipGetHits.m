function [hits1 hits2] = plusTipGetHits(groupData,saveDir,stringency,testID1,testID2)
% plusTipGetHits analyze statistics fields from multiple projects in groups
%
% SYNOPSIS:  plusTipPoolGroupData(groupList,saveDir,stringency,testID1,testID2
%
% INPUT:
% groupData : output of plusTipExtractGroupData
% saveDir   : path to output directory
% stringency: stringency to be applied on the result of the statistical
%             tests
% testID1   : id of the first statistical test to be performed on the data
%             (see discriminationMatrix)
% testID2   : id of the second statistical test to be performed on the data
%             (see discriminationMatrix)
%
% OUTPUT:
% hits1      : result of the first statistical test
% hits2      : result of the second statistical test
%
% Maria Bagonis, April 2011
% Sebastien Besson, July 2011

% Input check;
ip=inputParser;
ip.addRequired('groupData',@(x) isstruct(x) || iscell(x) || isempty(x));
ip.addRequired('saveDir',@(x) ischar(x) || isempty(x));
ip.addRequired('stringency',@isscalar);
ip.addRequired('testID1',@isscalar);
ip.addRequired('testID2',@isscalar);
ip.parse(groupData,saveDir,stringency,testID1,testID2);

% Call extract groupData function if empty or string input
if isempty(groupData) || iscell(groupData)
   groupData=plusTipExtractGroupData(groupData); 
end

% Launch interface if empty directory
if isempty(saveDir)
    saveDir=uigetdir(pwd,'Select output directory for group data.');
end

if isdir(saveDir), rmdir(saveDir,'s'); end
mkdir(saveDir);

%% Set Up DataStruct That Can Be Fed Into the Statistical Analysis

% Open a Project Data file to get fields from projData.stats
statsNames = fieldnames(groupData.stats{1}{1});
nGroups=numel(groupData.names);
dataStruct(nGroups)=struct();

%Set up the test structure required as input into the discriminationMatrix
% function
for i = 1:length(statsNames)
    testStruct.(statsNames{i}) = [testID1 testID2];
end

for iGroup = 1:nGroups   
    % Initiate A Structure That Will Hold All The Data for each group  
    for i = 1:length(statsNames)
        dataStruct(iGroup).(statsNames{i}) = ...
                cellfun(@(x) x.(statsNames{i}),groupData.stats{iGroup});
    end
    
    %% Calculate the Average and STD of the Per Cell Parameter Over All Projects (ie Cells) in Group  
    for i = 1:length(statsNames)
        mean_std.(statsNames{i})(iGroup,1) = mean(dataStruct(iGroup).(statsNames{i}));
        mean_std.(statsNames{i})(iGroup,2) = std(dataStruct(iGroup).(statsNames{i})); 
    end 
    
end % Finished Collecting Data for All Projects and all Groups

%% Save the Lists of Mean and Std Over all Cells in Group for Plotting.

title = cell(1,3);
title{1,1} = [];
title{1,2} = 'Mean Value Of Per Project Parameter: Calculated From The List of Value Corresponding to Each Project In Group (ie N = numProjects)';
title{1,3} = 'STD Of Per Project Parameter: Calculated From The List of Values Corresponding to Each Project In Group (ie N = numProjects)';

for i = 1:length(statsNames)
    mean_stdOverAllCellsInGroup.(statsNames{i}) = [groupData.names' num2cell(mean_std.(statsNames{i}))];
    mean_stdOverAllCellsInGroup.(statsNames{i}) = [title ; mean_stdOverAllCellsInGroup.(statsNames{i})];   
end

save([saveDir filesep 'dataStruct'],'dataStruct');
save([saveDir filesep 'MEAN_STD_OVER_ALL_CELLS_IN_GROUP'],'mean_stdOverAllCellsInGroup');

%% Feed This Data Into discrimination matrix
[discrimMats]  = discriminationMatrix(dataStruct,testStruct);

for i = 1:length(statsNames)
    pValues.(statsNames{i}) = num2cell(discrimMats.(statsNames{i}));
    
    %ADD TITLES
    pValues.(statsNames{i}) = [groupData.names' pValues.(statsNames{i})];
    pValues.(statsNames{i}) = [{[] groupData.names{:}} ; pValues.(statsNames{i})];
    
end

save([saveDir filesep 'discrimMat_PerCell'],'pValues');

%% Get Index of Hits And Generate Hits List
for i = 1:length(statsNames)
    hitsIdx.(statsNames{i}) = find(discrimMats.(statsNames{i})(:,1) < stringency ); % hits along
end

% initiate hitsList
titleHits = cell(1,4);
titleHits{1,1} = 'Group Name Of Hit';
titleHits{1,2} = 'Group Mean Of Per Project Parameter';
titleHits{1,3} = 'Group Mean Of Hit Condition / Group Mean Of Control Condition';
titleHits{1,4} = ['P-Value for Stats Test Number ' num2str(testID1)];

hits1=struct();
for iParam = 1:length(statsNames)
    
    nameOfField = statsNames{iParam};
    hitsList = cell(length(hitsIdx.(nameOfField)),4);
    %Write
    for iHit = 1:length(hitsIdx.(nameOfField))
        hitsList{iHit,1} = groupData.names{hitsIdx.(nameOfField)(iHit)};
        hitsList{iHit,2} = mean_std.(nameOfField)(hitsIdx.(nameOfField)(iHit),1);
        hitsList{iHit,3} = mean_std.(nameOfField)(hitsIdx.(nameOfField)(iHit),1)/mean_std.(nameOfField)(1,1);
        hitsList{iHit,4} = discrimMats.(nameOfField)(hitsIdx.(nameOfField)(iHit),1);
    end % iHit
    
    %Save in Larger Structure (but only if there is a hit)
    if ~isempty(hitsList), hits1.(nameOfField) = [titleHits ; hitsList]; end
end

save([saveDir filesep 'hitsTest1'],'hits1');

%% Do Same for Second Test
for iParam = 1:length(statsNames)
    hitsIdx.(statsNames{iParam}) = find(discrimMats.(statsNames{iParam})(1,:) < stringency );
end

% initiate hitsList

titleHits = cell(1,4);
titleHits{1,1} = 'Group Name Of Hit';
titleHits{1,2} = 'Group Mean Of Per Project Parameter';
titleHits{1,3} = 'Group Mean Of Hit Condition / Group Mean Of Control Condition';
titleHits{1,4} = ['P-Value for Stats Test Number ' num2str(testID2)];

hits2=struct();
for iParam = 1:length(statsNames)   
    nameOfField = statsNames{iParam};
    hitsList = cell(length(hitsIdx.(nameOfField)),4);
    %Write
    for iHit = 1:length(hitsIdx.(nameOfField))
        hitsList{iHit,1} = groupData.names{hitsIdx.(nameOfField)(iHit)};
        hitsList{iHit,2} = mean_std.(nameOfField)(hitsIdx.(nameOfField)(iHit),1);
        hitsList{iHit,3} = mean_std.(nameOfField)(hitsIdx.(nameOfField)(iHit),1)/mean_std.(nameOfField)(1,1);
        hitsList{iHit,4} = discrimMats.(nameOfField)(1,hitsIdx.(nameOfField)(iHit));
    end % iHit
    
    %Save in Larger Structure (but only if there is a hit)
    
    if isempty(hitsList) ~=  1
        hits2.(nameOfField) = [titleHits ; hitsList];
    end
    
end

save([saveDir filesep 'hitsTest2'],'hits2');


%% Plot results

% Rerieve the hit field of the datastruct only
nonHitFields = setxor(fieldnames(dataStruct(1)),fieldnames(hits1));
hitsDataStruct =arrayfun(@(x) rmfield(x,nonHitFields),dataStruct);
hitsStatsNames=fieldnames(hitsDataStruct);

% Plot barplots
for i=1:numel(hitsStatsNames)
    rawData={hitsDataStruct(:).(hitsStatsNames{i})};
    plotData= cellfun(@(x) mean(x),rawData);
    steData= cellfun(@(x) std(x)/sqrt(size(x,1)),rawData);
    f = figure;
    barplot2(plotData,steData,'YLabel',strrep(hitsStatsNames{i},'_',' '),...
        'XLabels',strrep(groupData.names,'_',' ' ),'Interpreter','none');
    print(f,'-dtiff', '-r300',[saveDir filesep 'histogram_' hitsStatsNames{i} '.tif']);
    close(f);
    
    plusTipPerCellPlot(hitsDataStruct,mean_stdOverAllCellsInGroup,hitsStatsNames{i},...
        groupData.names)
    print(gcf,'-dtiff', '-r300',[saveDir filesep 'plot_' hitsStatsNames{i} '.tif']);
    close(gcf);
    
end

end

