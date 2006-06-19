function [links, groupIdx, optData, linkData] = groupData(data, distanceFunction, groupingOptions, distanceFunctionParameters)
%GROUPDATA groups data with user-defined updates to the distance function
%
% SYNOPSIS: [links, groupIdx, optData, linkData] = ...
%           groupData(data,distanceFunction, groupingOptions, distanceFunctionParameters)
%
% INPUT data: data that is recognized by the distance function you specify
%
%       distanceFunction: name of the function that calculates the
%               distance between the objects you want to group. This
%               function also determines at what level the tree is to be
%               cut.
%               The distance function has to have the following inputs and
%               outputs:
%               [distance, associatedInfo, data, parameters] = distanceFunction(...
%                   data, parameters)
%               In: data - input data to groupData
%                   parameters - structure that contains at least the
%                    fields (created by groupData!)
%                      .group - either [], if the distance function has to be
%                           calculated for the entire data set, or [i,j,k],
%                           if sets i and j are grouped into a new set, k.
%                      .distanceOrder - ordered indices of active groups. If
%                           the last two of six data sets have been
%                           grouped (group = [5,6,7]), distanceOrder will
%                           be [1,2,3,4,7]. If group is empty, then
%                           distanceOrder will be empty, too!
%                      .groupIdx - as in output, but with zeros where no
%                           data is yet availale
%                      .groupInfo - associatedInfo for the current group
%                   any other parameters that are required for the
%                   operation of the distanceFunction can be passed as the
%                   structure distanceFunctionParameters. Importantly, the
%                   cutoff criterion has to be passed with this structure.
%
%               Out: distance - vector containing the distance between the
%                       sets. Element 1 is the distance between 1 and
%                         2, element 2 between 1 and 3 etc., where 1,2,3,
%                       refer to the order of sets in data (newly joined
%                       sets will be added at the end)
%                    associatedInfo - array with the same number of rows as
%                       distance contains elements. Place additional
%                       information about the potential links into the
%                       columns; they will be read into linkData.
%                       !! The first column of associatedInfo is 1 if the
%                       two grouped datasets are to be considered different
%                       (either because the distance is above a certain
%                       threshold, or because you want to cut off the tree
%                       at a specified level). Therefore
%                    data - the (updated) input data (this allows you to
%                       save calculations)
%                    parameters - same as input parameters
%
%       groupingOptions: optional structure with options for groupData
%               .verbose - whether or not to plot a dendrogram [{1}/0]
%               .labels  - cell array with labels for every data set
%
%		distanceFunctionParameters: optional argument that will be passed to the distance function
%
% OUTPUT links: link array as returned by linkage.m
%		 groupIdx: n-by-m array (m=n) that lists for each of the n
%			entries the group it belongs to at level m. Indices of level 1
%			are 1:n.
%        optData: additional output data depending on the groupingOptions
%           .collectedData(1:nGroups)
%                  .data(1:nData) - input data grouped according to the
%                       cutoff criterion
%                  .dataIdx(1:nData) - for every data set: index into
%                       original data
%           .collectedGrouping - how collectedData is grouped beyond the
%                                cutoff (links array)
%           .plotHandles [empty if verbose == 0]
%                   .figureH
%                   .axesH
%                   .lineH - handles of all the grouping lines of the
%                            dendrogram
%		 linkData: additional data describing the links that are returned
%			 by the distance function.
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 22-May-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%=====================
%% TEST INPUT
%=====================

% defaults
defaultOptions.verbose = true; % plot dendrogram
defaultOptions.labels = []; % no labels

if nargin < 2 || isempty(data) || isempty(distanceFunction)
    error('groupData needs at least two non-empty input arguments!')
end

if isempty(which(distanceFunction))
    error('The distance function %s needs to be a function on the matlab path',distanceFunction)
end

if nargin < 4
    distanceFunctionParameters = [];
end

% test options
if nargin < 3 || isempty(groupingOptions)
    groupingOptions = struct;
end
% fill fields
options = fieldnames(defaultOptions);
for fn = 1:length(options)
    if ~isfield(groupingOptions,options{fn}) || isempty(groupingOptions.(options{fn}))
        groupingOptions.(options{fn}) = defaultOptions.(options{fn});
    end
end

% preassign optData
optData = struct('collectedData',[],'collectedGrouping',[],...
    'plotHandles',[]);

%====================


%========================
%% PREPARE GROUPING LOOP
%========================

% get distances and preassign output

% run the distance function a first time. Don't make any assumptions about
% how the data look like - therefore don't supply distanceOrder
% set parameters
parameters = distanceFunctionParameters;
parameters.group = [];
parameters.distanceOrder = [];
parameters.groupIdx = [];
parameters.associatedInfo = [];
[distance, associatedInfo, data, parameters] = ...
    feval(distanceFunction, data, parameters);

% read number of sets
nDistances = length(distance(:));
nSets = (1+sqrt(1+8*nDistances))/2;
if floor(nSets) ~= nSets
    error('distanceFunction has to return one distance for every pair of sets in data!')
end

% calculate number of links
nLinks = nSets - 1;
links = zeros(nLinks,3);
parameters.groupIdx = zeros(nSets,nSets-1);
parameters.groupIdx(:,1) = (1:nSets)';


% read number of associatedInfo
% if isempty(associatedInfo)
%     addInfo = false;
%     linkData = [];
% else
%     addInfo = true;
if size(associatedInfo,1) ~= nDistances
    error('associatedInfo has to have as many rows as there are distances')
end
linkData = zeros(nLinks, size(associatedInfo,2));
% end


%==========================



%==========================
%% GROUPING LOOP
%==========================

% loop until all sets have been grouped, store info for all the linkages.
% Don't update distances in the last round

parameters.distanceOrder = (1:nSets)'; % create here because we never know
nGroups = nSets;
% remember original data
originalData = data;


for i = 1:nLinks

    % create list of indices
    [q,p] = find(tril(true(nGroups),-1));
    distanceIdx = parameters.distanceOrder([p,q]);
    % make sure distanceIdx always has the right orientation
    distanceIdx = reshape(distanceIdx,[],2);

    % group. Find minimum distance
    [minDist, minDistIdx] = min(distance);

    % store minimum distance, group
    newGroup = distanceIdx(minDistIdx,:);
    links(i,:) = [newGroup,minDist];
    %if addInfo
    linkData(i,:) = associatedInfo(minDistIdx,:);
    %end

    % if we aren't done yet: Update everything.
    if i < nLinks
        % remove the two grouped entries from distanceOrder, add new group
        nGroups = nGroups - 1;
        parameters.distanceOrder(ismember(parameters.distanceOrder,distanceIdx(minDistIdx,:))) = [];
        newGroupNumber = nSets + i; % we will always make a link!!
        parameters.distanceOrder(end+1) = newGroupNumber;
        parameters.groupIdx(:,i+1) = parameters.groupIdx(:,i);
        parameters.groupIdx(...
            ismember(parameters.groupIdx(:,i+1),newGroup),i+1) = ...
            newGroupNumber;

        % call distance function
        parameters.group = [newGroup,newGroupNumber];
        parameters.associatedInfo = linkData(i,:);
        [distance, associatedInfo, data, parameters] = ...
            feval(distanceFunction, data, parameters);


    end

end % grouping loop

% assign output
groupIdx = parameters.groupIdx;


%=====================
%% CUTOFF
%=====================

% check whether to cutoff at all - removed
%if ~isempty(groupingOptions.cutoff)
% check whether to cut according to distance criterion, or whether to
% make a set number of groups -- moved into distance function
%     if round(groupingOptions.cutoff) == -abs(groupingOptions.cutoff)
%         % take a given number of groups
%         nGroups = -groupingOptions.cutoff;
%
%         % get a pseudo-threshold
%         threshold = mean(links(end-nGroups+(1:2),3));
%
%     else
%         % read threshold
%         threshold = groupingOptions.cutoff;
%
%         % find nGroups
%         aboveThresholdIdx = links(:,3)>threshold;
%         nGroups = sum(aboveThresholdIdx) + 1;
%     end

% find number of links that are to be cut. If there is a non-cut link in
% between cut ones, consider that link cut, too.
cutIdx = find(linkData(:,1),1,'first');
if isempty(cutIdx)
    % if no significant difference, everything belongs to the same group
    nGroups = 1;
else
    % if 10 data sets, there are 9 links, and if the very last one is
    % significant, there are two groups
    nGroups = nLinks - cutIdx + 2;
end


if nGroups > 1
    % find n labels of groups
    groupLabels = unique(groupIdx(:,end-nGroups+2));

    % cut links
    optData.collectedGrouping = links(end-nGroups+2:end,:);

    % collectedData : puts data into nGroups groups by reading from groupIdx
    for iGroup = 1:nGroups
        % find which data we're putting into the group
        optData.collectedData(iGroup).dataIdx = ...
            find(groupIdx(:,end-nGroups+2) == groupLabels(iGroup));
        % collect the data
        optData.collectedData(iGroup).data = ...
            originalData(optData.collectedData(iGroup).dataIdx);

        % replace groupLabel with iGroup in groupLinks
        optData.collectedGrouping(...
            optData.collectedGrouping == groupLabels(iGroup)) = ...
            iGroup;
    end
    % color-threshold: look up smallest cut distance
    threshold = links(cutIdx,3);
else
    % only one group
    optData.collectedData.dataIdx = 1:nSets;
    optData.collectedData.data = originalData;
    optData.collectedGrouping = [];
    threshold = -1;
end


% removed because we can always attempt to groups
% else
%     % define threshold for dendrogram
%     threshold = -1;
%
%     % only one group
%         optData.collectedData.dataIdx = 1:nSets;
%         optData.collectedData.data = originalData;
%         optData.collectedGrouping = [];
%
% end


% show dendrogram if selected
if groupingOptions.verbose
    optData.plotHandles.figureH = figure;
    if ~isempty(links)
        if isempty(groupingOptions.labels)
            optData.plotHandles.lineH = dendrogram(links,0,...
                'orientation','right',...
                'COLORTHRESHOLD',threshold);
        else
            optData.plotHandles.lineH = dendrogram(links,0,...
                'labels',groupingOptions.labels,'orientation','right',...
                'COLORTHRESHOLD',threshold);
        end
    end
    optData.plotHandles.axesH = gca;
end






