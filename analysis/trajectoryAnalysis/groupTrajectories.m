function [nClusters,clusterList,comparison,dendrogramHandles] = groupTrajectories(trajectoryData,weight)
%GROUPTRAJECTORIES tries to find groups of trajectories within a condition
%
% SYNOPSIS
%
% INPUT taData: Output from trajectoryAnalysis from a single condition
%
% OUTPUT nClusters : number of different clusters
%        clusterList : for every trajectory: which cluster it belongs to
%        comparison : discrimination matrices
%        dendrogramHandles : handles to the lines in the dendrogram
%
% c: jonas 11/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
% constants
significanceLevel = 0.05;

if ~isstruct(trajectoryData)
    if nargin == 1 || isempty(weight)
        weight = [1,1];
    end
    nTrajectories = length(trajectoryData);
    upperIndices = find(triu(ones(nTrajectories),1));
lowerIndices = find(tril(ones(nTrajectories),-1));
pValues = zeros(2,length(upperIndices));
pValues(1,:) = trajectoryData(lowerIndices)'/weight(1);
pValues(2,:) = trajectoryData(upperIndices)'/weight(2);

else

% count individual trajectories
nTrajectories = length(trajectoryData.individualStatistics);
if nTrajectories == 2
    error('You need at least 3 trajectories to group!')
end

% collect data
compStruct = buildCompStruct(trajectoryData, nTrajectories);
toc
% compare data. Use default test (ttest for mean, KS for corrected
% distributions)
individualComparisons = discriminationMatrix(compStruct);
toc
% read individualComparisons. Build indexMatrices first by reading the indices of
% upper and lower triangular matrices
upperIndices = find(triu(ones(nTrajectories),1));
lowerIndices = find(tril(ones(nTrajectories),-1));

pValues = zeros(6,nTrajectories * (nTrajectories-1)/2);

pValues(1,:) = individualComparisons.growthSpeeds(lowerIndices)';
pValues(2,:) = individualComparisons.growthSpeeds(upperIndices)';
pValues(3,:) = individualComparisons.shrinkageSpeeds(lowerIndices)';
pValues(4,:) = individualComparisons.shrinkageSpeeds(upperIndices)';
pValues(5,:) = individualComparisons.distanceSigma(lowerIndices)';
pValues(6,:) = individualComparisons.distanceSigma(upperIndices)';
pValues(7,:) = individualComparisons.growthTimes(lowerIndices)';
pValues(8,:) = individualComparisons.growthTimes(upperIndices)';
pValues(9,:) = individualComparisons.shrinkageTimes(lowerIndices)';
pValues(10,:) = individualComparisons.shrinkageTimes(upperIndices)';


% use only comparison of means - distributions need more data to make
% sense!
pValues = pValues([1,3,7,9],:);

end

% large difference = small p-value. Take 1-p for distances
pValues = max(max(pValues(:)),1)-pValues;


% calculate euclidean distance
pDistance = sqrt(sum(pValues.^2,1));

% make links. Use shortest distance
links = linkage(pDistance,'average');

% visualize
z = squareform(pDistance);
minDist = min(pDistance);
maxDist = max(pDistance);
uiViewPanel,imshow(-z,[-maxDist,-minDist]);
figure
[dendrogramHandles,t,p]=dendrogram(links,0);
uiViewPanel,imshow(-z(p,p),[-maxDist,-minDist]);

% cluster into two groups. Check for difference
nClusters = 2;
clusterList = cluster(links,nClusters);


% compare
groupComparisons = compareGroups(compStruct,clusterList);

% if different: continue. If not, return
if ~isDifferent(groupComparisons, significanceLevel)
    oldComparison = [];
    oldClusterList = [];
else
    % different. Loop till we find no more difference. Whenever a cluster
    % could not be splitted anymore, it becomes locked
    allLocked = 0;
    lockedList = [];
    
    while ~allLocked && nClusters < nTrajectories
        % remember last round's result
        oldClusterList = clusterList;
        oldComparison = groupComparisons;
        % try one more cluster
        nClusters = nClusters + 1;
        clusterList = cluster(links,nClusters);
        % overwrite locked indices
        for ll = lockedList
            clusterList(clusterList == ll) = ll;
        end
        % try to find new cluster
        newIdx = find(clusterList == nClusters);
        % if it was in a locked cluster, newIdx is empty
        if isempty(newIdx)
            % try one more cluster (there is still an unlocked cluster!
        else
            % compare the new groups
            groupComparisons = compareGroups(compStruct,clusterList);
            % if different: goto next. Else, add the cluster we just split
            % to the locked list
            if isDifferent(groupComparisons, significanceLevel)
                % continue (still check for allLocked - there could be only
                % singles left)
            else
                % find cluster to lock and remove new cluster from list
                oldNumber = oldClusterList(newIdx(1));
                lockedList = [lockedList,oldNumber]; 
                clusterList(newIdx) = oldNumber;
            end
            
            % check for singles
            [idx, count] = countEntries(clusterList);
            countOne = count == 1;
            lockedList = unique([lockedList,idx(countOne)']);
            
            % check for allLocked
            allLocked = all(ismember(clusterList,lockedList));
        end
            
    end
    
end
toc
% return output
nClusters = nClusters - 1;
comparison = oldComparison;
clusterList = oldClusterList;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function groupComparisons = compareGroups(compStruct,clusterList)
%===================================================================
% COMPARE GROUPS
%===================================================================

groupNumbers = unique(clusterList);
nGroups = length(groupNumbers);
compNames = fieldnames(compStruct);
nNames = length(compNames);

compStructNew(1:nGroups) = cell2struct(cell(size(compNames)),compNames,1);
% loop to fill new groupComparisons matrix
for iGroup = 1:nGroups
    % find indices belonging to cluster
    groupIdx = clusterList == groupNumbers(iGroup);
    % loop through fields to catenate
    for iName = 1:nNames
        compStructNew(iGroup).(compNames{iName}) = ...
            cat(1,compStruct(groupIdx).(compNames{iName}));
    end
end

% compare
groupComparisons = discriminationMatrix(compStructNew);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = isDifferent(groupComparisons, significanceLevel)
%===================================================================
% CHECK WHETHER CONDITIONS ARE DIFFERENT
%===================================================================

% for out to become 1, no two conditions can be the same in all conditions

% get basic info about input
tests = fieldnames(groupComparisons);
nTests = length(tests);
nConditions = size(groupComparisons.(tests{1}),1);
upperIndices = find(triu(ones(nConditions),1));
lowerIndices = find(tril(ones(nConditions),-1));

pValues = zeros(2 * nTests, nConditions * (nConditions-1)/2);

% read pValues into matrix
for i = 1:nTests
    ip = (i-1)*2;
    pValues(ip+1,:) = groupComparisons.(tests{i})(lowerIndices)';
    pValues(ip+2,:) = groupComparisons.(tests{i})(upperIndices)';
end

% check for significance
pValues = pValues < significanceLevel;

% test only means?

% check whether there is no groupComparisons that is never significant
out = all(pValues,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compStruct = buildCompStruct(trajectoryData, nTrajectories)
%====================================================================
% BUILD COMPARISON STRUCTURE
%====================================================================


% collect growth speeds, shrinkage speeds and distanceSigma
% (abs(dist-mean))
compStruct(1:nTrajectories) = struct('growthSpeeds',[],'shrinkageSpeeds',[],...
    'distanceSigma',[]);

for ti = 1:nTrajectories
    % read data. Calculate indices first, though
    [growthIdx,growthGroupIdx] = ...
        trajectoryAnalysisMainCalcStatsFindGroups(...
        trajectoryData.individualStatistics(ti).dataListGroup,1);
    [shrinkageIdx,shrinkageGroupIdx] = ...
        trajectoryAnalysisMainCalcStatsFindGroups(...
        trajectoryData.individualStatistics(ti).dataListGroup,2);


    % speeds. Speeds persisting over longer periods of time are
    % repeated per interval
    if ~isempty(growthIdx)
        compStruct(ti).growthSpeeds = repeatEntries(trajectoryData.individualStatistics(ti).dataListGroup(...
            growthIdx,4), diff(trajectoryData.individualStatistics(ti).dataListGroup(...
            growthIdx,[1:2]),1,2))*60;
    end
    if ~isempty(shrinkageIdx)
        compStruct(ti).shrinkageSpeeds = repeatEntries(trajectoryData.individualStatistics(ti).dataListGroup(...
            shrinkageIdx,4), diff(trajectoryData.individualStatistics(ti).dataListGroup(...
            shrinkageIdx,[1:2]),1,2))*60;
    end
    
    if ~isempty(growthGroupIdx)
        compStruct(ti).growthTimes = diff([trajectoryData.individualStatistics(ti).dataListGroup(...
            growthGroupIdx(:,1),1),trajectoryData.individualStatistics(ti).dataListGroup(...
            growthGroupIdx(:,2),2)],1,2);
    end

    if ~isempty(shrinkageGroupIdx)
        compStruct(ti).shrinkageTimes = diff([trajectoryData.individualStatistics(ti).dataListGroup(...
            shrinkageGroupIdx(:,1),1),trajectoryData.individualStatistics(ti).dataListGroup(...
            shrinkageGroupIdx(:,2),2)],1,2);
    end
    
    

    % don't compare frequencies
%     % growth/ shrinkage group times
%     if ~isempty(growthGroupIdx)
%         compStruct(ti).growthTimes = diff([trajectoryData.individualStatistics(ti).dataListGroup(...
%             growthGroupIdx(:,1),1),trajectoryData.individualStatistics(ti).dataListGroup(...
%             growthGroupIdx(:,2),2)],1,2);
%     end
%      if ~isempty(shrinkageGroupIdx)
%         compStruct(ti).shrinkageTimes = diff([trajectoryData.individualStatistics(ti).dataListGroup(...
%             shrinkageGroupIdx(:,1),1),trajectoryData.individualStatistics(ti).dataListGroup(...
%             shrinkageGroupIdx(:,2),2)],1,2);
%     end
 
    % distance mean - all distance


    % distance variation - distance minus mean of distance
    compStruct(ti).distanceSigma = ...
        (trajectoryData.individualStatistics(ti).addedStats.distance(:,1) - ...
        trajectoryData.individualStatistics(ti).summary.distanceMean(1)).^2;
end % loop trajectories