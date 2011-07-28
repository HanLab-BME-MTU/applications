function [hits1 hits2] = plusTipGetHits(groupList,saveDir,stringency,testID1,testID2,varargin)
% plusTipGetHits analyze statistics fields from multiple projects in groups
%
% SYNOPSIS:  plusTipPoolGroupData(groupList,saveDir,stringency,testID1,testID2
%
% INPUT:
% groupList : cell array of structures containing the output of the
%             plusTipTracker postprocessing
%             output of plusTipExtractGroupData
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

if nargin<1 || isempty(saveDir)
    saveDir=uigetdir(pwd,'Select output directory for hit files.');
end

projGroupName=groupList(:,1);
projGroupDir= groupList(:,2);

% fix the names if there are spaces or hyphens and append prefix 'grp'
projGroupName=cellfun(@(x) strrep(x,'-','_'),projGroupName,'uniformoutput',0);
projGroupName=cellfun(@(x) strrep(x,' ','_'),projGroupName,'uniformoutput',0);
projGroupName=cellfun(@(x) ['grp_' x],projGroupName,'uniformoutput',0);


% count unique groups and keep them in order of the original list
[btwGrpNames,m,projGroupIdx] = unique(projGroupName);
[b,idx]=sort(m);
btwGrpNames=btwGrpNames(idx);

groupNames = cell(length(btwGrpNames),1);
for iGroup = 1:length(btwGrpNames)
    groupNames(iGroup,1)= btwGrpNames(iGroup,1);
end

groupNamesHor = cell(1,length(btwGrpNames));
for iGroup = 1:length(btwGrpNames)+1
    if iGroup == 1
        groupNamesHor{1,1} = [];
    else
        groupNamesHor(1,iGroup) = btwGrpNames(iGroup-1,1);
    end
end
%% Set Up DataStruct That Can Be Fed Into the Statistical Analysis


% Open a Project Data file to get fields from projData.stats
temp = load([char(projGroupDir(1)) filesep 'meta' filesep 'projData']);

paramNames = fieldnames(temp.projData.stats);

%Set up the test structure required as input into the discriminationMatrix
% function

for iParam = 1:length(paramNames)
    testStruct.(char(paramNames(iParam))) = [testID1 testID2];
end

for iGroup = 1:length(btwGrpNames)
    
    % get indices of projects in iGroup
    tempIdx= strmatch(btwGrpNames(iGroup),projGroupName,'exact');
    
    % Initiate A Structure That Will Hold All The Data for each group
    
    for iParam = 1:length(paramNames)
        dataStruct(iGroup).(char(paramNames(iParam))) = zeros(length(tempIdx),1);
    end
    
    % Loop through each project in group and collect values for each parameter
    % from the projData.stats file
    
    for iProj = 1:length(tempIdx)
        temp = load([projGroupDir{tempIdx(iProj)} filesep 'meta' filesep 'projData']);
        
        for iParam = 1:length(paramNames)
            dataStruct(iGroup).(char(paramNames(iParam)))(iProj,1) = getfield(temp.projData.stats, char(paramNames(iParam)));
        end
    end % Finished collecting data from all projects in a given group
    
    
    %% Calculate the Average and STD of the Per Cell Parameter Over All Projects (ie Cells) in Group
    
    for iParam = 1:length(paramNames)
        
        mean_std.(char(paramNames(iParam)))(iGroup,1) = mean(dataStruct(iGroup).(char(paramNames(iParam))));
        mean_std.(char(paramNames(iParam)))(iGroup,2) = std(dataStruct(iGroup).(char(paramNames(iParam))));
        
    end % iParam
    
end % Finished Collecting Data for All Projects and all Groups
%% Save the Lists of Mean and Std Over all Cells in Group for Plotting.

title = cell(1,3);
title{1,1} = [];
title{1,2} = 'Mean Value Of Per Project Parameter: Calculated From The List of Value Corresponding to Each Project In Group (ie N = numProjects)';
title{1,3} = 'STD Of Per Project Parameter: Calculated From The List of Values Corresponding to Each Project In Group (ie N = numProjects)';

for iParam = 1:length(paramNames)
    mean_stdOverAllCellsInGroup.(char(paramNames(iParam))) = [groupNames num2cell(mean_std.(char(paramNames(iParam))))];
    mean_stdOverAllCellsInGroup.(char(paramNames(iParam))) = [title ; mean_stdOverAllCellsInGroup.(char(paramNames(iParam)))];
    
end

save([saveDir filesep 'dataStruct'],'dataStruct');
save([saveDir filesep 'MEAN_STD_OVER_ALL_CELLS_IN_GROUP'],'mean_stdOverAllCellsInGroup');

%% Feed This Data Into discrimination matrix
[discrimMats]  = discriminationMatrix(dataStruct,testStruct);

for iParam = 1:length(paramNames)
    pValues.(char(paramNames(iParam))) = num2cell(discrimMats.(char(paramNames(iParam))));
    
    %ADD TITLES
    pValues.(char(paramNames(iParam))) = [groupNames pValues.(char(paramNames(iParam)))];
    pValues.(char(paramNames(iParam))) = [groupNamesHor ; pValues.(char(paramNames(iParam)))];
    
end

save([saveDir filesep 'discrimMat_PerCell'],'pValues');

%% Get Index of Hits And Generate Hits List
for iParam = 1:length(paramNames)
    
    hitsIdx.(char(paramNames(iParam))) = find(discrimMats.(char(paramNames(iParam)))(:,1) < stringency ); % hits along
end

% initiate hitsList

titleHits = cell(1,4);
titleHits{1,1} = 'Group Name Of Hit';
titleHits{1,2} = 'Group Mean Of Per Project Parameter';
titleHits{1,3} = 'Group Mean Of Hit Condition / Group Mean Of Control Condition';
titleHits{1,4} = ['P-Value for Stats Test Number ' num2str(testID1)];


for iParam = 1:length(paramNames)
    
    nameOfField = char(paramNames(iParam));
    hitsList = cell(length(hitsIdx.(nameOfField)),4);
    %Write
    for iHit = 1:length(hitsIdx.(nameOfField))
        hitsList{iHit,1} = char(groupNames(hitsIdx.(nameOfField)(iHit),1));
        hitsList{iHit,2} = mean_std.(nameOfField)(hitsIdx.(nameOfField)(iHit),1);
        hitsList{iHit,3} = mean_std.(nameOfField)(hitsIdx.(nameOfField)(iHit),1)/mean_std.(nameOfField)(1,1);
        hitsList{iHit,4} = discrimMats.(nameOfField)(hitsIdx.(nameOfField)(iHit),1);
    end % iHit
    
    %Save in Larger Structure (but only if there is a hit)
    
    
    
    if isempty(hitsList) ~= 1
        hits1.(nameOfField) = [titleHits ; hitsList];
    end
end

save([saveDir filesep 'hitsTest1'],'hits1');

%% Do Same for Second Test
for iParam = 1:length(paramNames)
    
    hitsIdx.(char(paramNames(iParam))) = find(discrimMats.(char(paramNames(iParam)))(1,:) < stringency );
end

% initiate hitsList

titleHits = cell(1,4);
titleHits{1,1} = 'Group Name Of Hit';
titleHits{1,2} = 'Group Mean Of Per Project Parameter';
titleHits{1,3} = 'Group Mean Of Hit Condition / Group Mean Of Control Condition';
titleHits{1,4} = ['P-Value for Stats Test Number ' num2str(testID2)];


for iParam = 1:length(paramNames)
    
    nameOfField = char(paramNames(iParam));
    hitsList = cell(length(hitsIdx.(nameOfField)),4);
    %Write
    for iHit = 1:length(hitsIdx.(nameOfField))
        hitsList{iHit,1} = char(groupNames(hitsIdx.(nameOfField)(iHit),1));
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

nonHitFields = setxor(fieldnames(dataStruct(1)),fieldnames(hits1));
plotDataStruct =arrayfun(@(x) rmfield(x,nonHitFields),dataStruct);
paramNames=fieldnames(plotDataStruct);
for i=1:numel(paramNames)
    rawData={plotDataStruct(:).(paramNames{i})};
    plotData= cellfun(@(x) mean(x),rawData);
    steData= cellfun(@(x) std(x)/sqrt(size(x,1)),rawData);
    f = figure;
    barplot2(plotData,steData,'YLabel',paramNames{i},...
        'XLabels',groupNames,'Interpreter','none');
    print(f,'-dtiff', '-r300',[saveDir filesep 'histogram_' paramNames{i} '.tif']);
    close(f);
end

end

