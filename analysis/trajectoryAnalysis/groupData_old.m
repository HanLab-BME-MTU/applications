function [finalComparison, clusterIdxOut, links] = groupData(data,experiment)
%GROUPDATA hierarchically groups data into significantly different clusters
%
% groupData provides the framework for the hierarchical clustering and the
% significance testing. Since tests and significance levels are strongly
% problem-dependent, reading out the data from the input, the combination
% of data and the statistical comparison will be contained within
% switch-case statements that you will have to customize for your own
% experiment.
%
% SYNOPSIS [finalComparison, clusterIdx] = groupData(data,experiment)
%
% INPUT    data: data you want to compare
%          experiment (scalar): selects your experiment. Please add any
%               new experiment to the list here:
%                1: find groups of trajectories within one condition
%                   for the chromdyn project
%                2: find groups of trajectories within one condition for
%                   chromdyn project. Use instantaneous velocities
%
% OUTPUT   finalComparison: discrimination matrix for the comparison of the
%               significant clusters
%          clusterIdxOut  : index of cluster for every entry in data
%          links          : links matrix (the way you'd create it with
%                           linkage
%
% c: jonas 12/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%====================
% TEST INPUT
%====================

% check nargin (isscalar returns 0 if empty)
if nargin < 2 || isempty(data) || ~isscalar(experiment)
    error('insufficent or wrong input arguments!')
end

% read experiment-dependent parameters
switch experiment
    case 1
        significanceLevel = 0.05;
    case 2
        significanceLevel = 0.01;
    otherwise
        error('experiment %i is not defined yet!',experiment)
end

%==================


%==================
% READ DATA
%==================

% read data into compStruct (input for discriminationMatrix). compStruct is
% a structure array with a field for every aspect of the data that should
% be compared, for example: compStruct(1:n).growthSpeed,
% compStruct(1:n).shrinkageSpeed etc., where n is the number of groups of
% data to compare.
% If you don't provide the data as a compStruct, edit the subfunction that
% does the transformation.

switch experiment
    case {1,2} % if you use the subfunction: add your experimentID to the cell
        compStruct = readDataIntoCompStruct(data,experiment);
    otherwise
        compStruct = data;
end

% free memory
data = [];

%=================


%=================
% CLUSTER DATA
%=================

% if we have n data to compare, there will be n-1 groupings. Loop. Start by
% calculating the discrimination matrix for the input. Calculate the lowest
% distance between clusters (if multiple simultaneous comparisons) and
% select the clusters to be joined. Also calculate whether the current
% clusters were all significantly different. Store in links. Rinse, repeat.
% Variables used:
%   compStruct(1:n)   : input for discriminationMatrix
%   comparison(1:n-1) : result from discriminationMatrix. Basically, it is
%                       a matrix where at (1,2) there is the p-value of the
%                       comparison between group 1 and 2 (you don't need to
%                       use discriminationMatrix, of course)
%   clusterIdx        : for every original group: To which cluster in the
%                       current comparison it belongs
%   clusterConversion : cluster number in links for every row/col in the
%                       comparison-matrix
%   links             : clustering matrix analoguous to the output from
%                       linkage.m, with the difference that it is always
%                       indicated whether the clusters were all
%                       significantly different in col. 4 (links normally
%                       has 3 cols: cluster 1, cluster 2, distance)




% --- preassign matrices

nComparisons = length(compStruct)-1;

% comparison
compNames = fieldnames(compStruct);
nFields = length(compNames);
comparison(1:nComparisons) = cell2struct(cell(nFields,1),compNames,1);

% compareIdx, clusterIdx
clusterIdx = cell(1,nComparisons + 1); % in the loop, we always assign the next entry
clusterIdx{1} = (1:nComparisons + 1)';
clusterConversion = clusterIdx; % they are the same only for the first element!

% links
links = zeros(nComparisons,4);
tic
% --- loop
for iComparison = 1:nComparisons

    % assign some general values
    currentNComparisons = nComparisons - iComparison + 1;
    % upper/lowerIndices are indices to the p-values in the comparison
    % matrices
    upperIndices = find(triu(ones(currentNComparisons+1),1));
    lowerIndices = find(tril(ones(currentNComparisons+1),-1));
    % upper/lowerComparisons are the indices of the groups that are
    % compared at upperIndices. For example, upperComparisons(:,1) = [1;2],
    % while lowerComparisons(:,1) = [2;1]. Careful, though: the indices are
    % ordered differently
    [rows, cols] = ndgrid(1:currentNComparisons+1);
    upperComparisons = [rows(upperIndices)';cols(upperIndices)'];
    lowerComparisons = [rows(lowerIndices)';cols(lowerIndices)'];
    % to compensate for the different ordering, we sort here
    [dummy,upperSortIdx] = ...
        ismember(lowerComparisons',upperComparisons([2,1],:)','rows');
    upperIndices = upperIndices(upperSortIdx);
    upperComparisons = upperComparisons(:,upperSortIdx);
    % now, upperIndices and lowerIndices are still different, but upper and
    % lower comparisons are the same, which is good. Assign them to
    % comparisonIndices
    comparisonIndices = upperComparisons;
    

    % compare data
    switch experiment
        case {1,2}
            % default comparison: ttest for means, kstest for distributions
            comparison(iComparison) = ...
                discriminationMatrix(compStruct);
    end

    % create distance list (as from pDist.m) - store as variable "distance"
    % Since we're dealing with p-values here, have distance be 1-p
    % check for whether all is different (define significanceLevel at the
    % beginning of the program, please).
    switch experiment
        case 1
            % for distance: we use only the means. Calculate euclidean
            % distance.
            % we assume there to be a difference if at least one of the
            % tests reports a difference. Don't use distributions if less
            % than 5 trajectories grouped, don't compare single
            % trajectories at all
            pValuesUpper = ones(nFields,length(upperIndices));
            pValuesLower = ones(nFields,length(lowerIndices));

            % read pValues
            for iField = 1:nFields % counting down requires no preassignment
                pValuesLower(iField,:) = ...
                    comparison(iComparison).(compNames{iField})(lowerIndices);
            end
            for iField = 1:2
                pValuesUpper(iField,:) = ...
                    comparison(iComparison).(compNames{iField})(upperIndices);
            end
            % calculate distance on comparisons of the means only
            distance = sqrt(sum((1-pValuesLower(1:2,:)).^2,1));

            % find whether all are different.
            
             % uniqueEntries = 1:currentNComparisons
            [uniqueEntries,numberOfOccurences] = ...
                countEntries(clusterIdx{iComparison});
            
            % mask uppers (distribution comparison) wherever there are less
            % than five trajectories - set 1, because these don't count for
            % check whether significant.
            notFiveEntries = uniqueEntries(numberOfOccurences < 5);
            notFiveIndicesUpper = ...
                find(sum(ismember(upperComparisons,notFiveEntries),1));
            pValuesUpper(:,notFiveIndicesUpper) = 1;

            % mask comparisons with single entries by setting p=0 (making
            % them significant no matter what)

            % singleIndices: uniqueEntries where numberOfOccurences = 1
            singleEntries = uniqueEntries(numberOfOccurences == 1);
            singleIndicesUpper =...
                find(sum(ismember(upperComparisons,singleEntries),1));
            singleIndicesLower =...
                find(sum(ismember(lowerComparisons,singleEntries),1));
            pValuesUpper(:,singleIndicesUpper) = 0;
            pValuesLower(:,singleIndicesLower) = 0;

            

            % check whether all are different
            allDifferent = all(...
                any(...
                pValuesUpper < significanceLevel |...
                pValuesLower < significanceLevel,1));
        case 2
            % use p-values for means AND distributions for everything
             pValuesUpper = ones(nFields,length(upperIndices));
            pValuesLower = ones(nFields,length(lowerIndices));

            % read pValues
            for iField = 1:nFields % counting down requires no preassignment
                pValuesLower(iField,:) = ...
                    comparison(iComparison).(compNames{iField})(lowerIndices);
            end
            for iField = 1:nFields
                pValuesUpper(iField,:) = ...
                    comparison(iComparison).(compNames{iField})(upperIndices);
            end
            % calculate distance on comparisons of the means only
            distance = sqrt(sum((1-[pValuesLower;pValuesUpper]).^2,1));

            % find whether all are different.
            
             % uniqueEntries = 1:currentNComparisons
            [uniqueEntries,numberOfOccurences] = ...
                countEntries(clusterIdx{iComparison});
            

            % mask comparisons with single entries by setting p=0 (making
            % them significant no matter what)

            % singleIndices: uniqueEntries where numberOfOccurences = 1
            singleEntries = uniqueEntries(numberOfOccurences == 1);
            singleIndicesUpper =...
                find(sum(ismember(upperComparisons,singleEntries),1));
            singleIndicesLower =...
                find(sum(ismember(lowerComparisons,singleEntries),1));
            pValuesUpper(:,singleIndicesUpper) = 0;
            pValuesLower(:,singleIndicesLower) = 0;

            

            % check whether all are different
            allDifferent = all(...
                any(...
                pValuesUpper < significanceLevel |...
                pValuesLower < significanceLevel,1));
        otherwise           
           
    end % switch for distances and allDifferent

    %--- write links
    [smallestDistance,smallestDistanceIdx] = min(distance);
    linkedGroups = comparisonIndices(:,smallestDistanceIdx);
    cc = clusterConversion{iComparison};
    links(iComparison,:) = [cc(linkedGroups(1)),...
        cc(linkedGroups(2)),...
        smallestDistance,...
        allDifferent];

    % --- update compStruct, clusterConversion, clusterIdx
    % The two groups to be combined (linkedGroups) will be joined and
    % appended at the end of  compStruct, and the entries corresponding to
    % the two separate linkedGroups will be removed. This joining can be
    % problem-specific.
    % A new cluster-number will be added at the end of clusterConversion,
    % and linkedGroups will be removed.
    % All entries in clusterIdx belonging to linkedGroups will be
    % replaced with currentNComparisons-1 (because they have been joined at
    % the end of the list). From the other entries, 0, 1 or 2 will be
    % subtracted, depending on wheter no, one, or two of the linkedGroups
    % had smaller indices.
    switch experiment
        case {1,2}
            % for trajectories, joining is simple
            % (we're reassigning compStruct here. Since we only need to
            % reassign once per big loop, the loss of time should not be
            % severe)
            for iField = 1:nFields
                compStruct(currentNComparisons + 2).(compNames{iField}) = ...
                    [compStruct(linkedGroups(1)).(compNames{iField});...
                    compStruct(linkedGroups(1)).(compNames{iField})];
            end
            % remove old entries
            compStruct(linkedGroups) = [];
    end

    % update clusterConversion
    cc = clusterConversion{iComparison};
    cc(end+1) = max(cc) + 1;
    cc(linkedGroups) = [];
    clusterConversion{iComparison + 1} = cc;

    % update clusterIdx
    ci = clusterIdx{iComparison};
    linkedIdx = ismember(ci,linkedGroups);
    notLinkedIdx = ~linkedIdx;
    ci(linkedIdx) = currentNComparisons - 1; % number of comparisons next time
    delta = (ci(notLinkedIdx) > linkedGroups(1)) + ...
        (ci(notLinkedIdx) > linkedGroups(2));
    ci(notLinkedIdx) = ci(notLinkedIdx) - delta;
    clusterIdx{iComparison + 1} = ci;

    
    toc
end % loop iComparisons: cluster data

%================================
% FIND SIGNIFICANT CLUSTERS
%================================

% Because we noted in the links whether all data is significantly
% different, this is easy. Just move up links(:,4) until we find the first
% 0. This is the last link we do.

insignificant = find(~links(:,4));
if ~isempty(insignificant)
    lastGoodLink = insignificant(end);
else
    lastGoodLink = 0;
end

%================================
% OUTPUT
%================================

% assign output
clusterIdxOut = clusterIdx{lastGoodLink + 1};
finalComparison = comparison(lastGoodLink + 1);
%links = links(:,1:3);

% plot
figure;
h = dendrogram(links(:,1:3),0);
set(h(lastGoodLink + 1 : end),'Color','r')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compStruct = readDataIntoCompStruct(data,experiment)
%====================================================================
% BUILD COMPARISON STRUCTURE
%====================================================================
switch experiment
    case 1
        nTrajectories = length(data.individualStatistics);
        % collect growth speeds, shrinkage speeds and distanceSigma
        % (abs(dist-mean))
        compStruct(1:nTrajectories) = struct('growthSpeeds',[],'shrinkageSpeeds',[],...
            'growthTimes',[],'shrinkageTimes',[]);

        for ti = 1:nTrajectories
            % read data. Calculate indices first, though
            [growthIdx,growthGroupIdx] = ...
                trajectoryAnalysisMainCalcStatsFindGroups(...
                data.individualStatistics(ti).dataListGroup,1);
            [shrinkageIdx,shrinkageGroupIdx] = ...
                trajectoryAnalysisMainCalcStatsFindGroups(...
                data.individualStatistics(ti).dataListGroup,2);


            % speeds. Speeds persisting over longer periods of time are
            % repeated per interval
            if ~isempty(growthIdx)
                compStruct(ti).growthSpeeds = repeatEntries(data.individualStatistics(ti).dataListGroup(...
                    growthIdx,4), diff(data.individualStatistics(ti).dataListGroup(...
                    growthIdx,[1:2]),1,2))*60;
            end
            if ~isempty(shrinkageIdx)
                compStruct(ti).shrinkageSpeeds = repeatEntries(data.individualStatistics(ti).dataListGroup(...
                    shrinkageIdx,4), diff(data.individualStatistics(ti).dataListGroup(...
                    shrinkageIdx,[1:2]),1,2))*60;
            end

            if ~isempty(growthGroupIdx)
                compStruct(ti).growthTimes = diff([data.individualStatistics(ti).dataListGroup(...
                    growthGroupIdx(:,1),1),data.individualStatistics(ti).dataListGroup(...
                    growthGroupIdx(:,2),2)],1,2);
            end

            if ~isempty(shrinkageGroupIdx)
                compStruct(ti).shrinkageTimes = diff([data.individualStatistics(ti).dataListGroup(...
                    shrinkageGroupIdx(:,1),1),data.individualStatistics(ti).dataListGroup(...
                    shrinkageGroupIdx(:,2),2)],1,2);
            end
        end
    case 2
        nTrajectories = length(data.individualStatistics);
        % collect growth speeds, shrinkage speeds (inst. velocities)
        compStruct(1:nTrajectories) = struct('growthSpeeds',[],'shrinkageSpeeds',[]);
        
        for ti = 1:nTrajectories
            % read delta distance
            deltaDistance = data.individualStatistics(ti).dataListSeed(:,9);
            % remove deleted frames
            deltaDistance(...
                data.individualStatistics(ti).dataListSeed(:,3) < 0) = [];
            % put positive delta into growthSpeeds, negative into shrinkage
            compStruct(ti).growthSpeeds = ...
                deltaDistance(deltaDistance > 0);
            compStruct(ti).shrinkageSpeeds = ...
                deltaDistance(deltaDistance < 0);
            
        end
        
    otherwise
        error('reading data not defined')
end
