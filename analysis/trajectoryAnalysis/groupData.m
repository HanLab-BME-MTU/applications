function [links, groupIdx, linkData] = groupData(data, distanceFunction, distanceFunctionParameters)
%GROUPDATA groups data with user-defined updates to the distance function
%
% SYNOPSIS: [links, groupIdx, linkData] = groupData(data, distanceFunction, distanceFunctionParameters)
%
% INPUT data: data that is recognized by the distance function you specify
%       distanceFunction: name of the function that calculates the
%               distance between the objects you want to group.
%               The distance function has to have the following inputs and
%               outputs:
%               [distance, associatedInfo, data, parameters] = distanceFunction(...
%                   data, parameters)
%               In: data - input data to groupData
%                   parameters - structure that contains at least the
%                   fields (created by groupData!)
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
%                   any other parameters that are required for the
%                   operation of the distanceFunction can be passed as the
%                   structure distanceFunctionParameters
%                   
%               Out: distance - vector containing the distance between the
%                       sets. Element 1 is the distance between 1 and
%                       2, element 2 between 1 and 3 etc., where 1,2,3,
%                       refer to the order of sets in data (newly joined
%                       sets will be added at the end)
%                    associatedInfo - array with the same number of rows as
%                       distance contains elements. Place additional
%                       information about the potential links into the
%                       columns; they will be read into linkData
%                    data - the (updated) input data (this allows you to
%                       save calculations)
%                    parameters - same as input parameters
%		distanceFunctionParameters: optional argument that will be passed to the distance function
%
% OUTPUT links: link array as returned by linkage.m
%		 groupIdx: n-by-m array (m=n) that lists for each of the n
%			entries the group it belongs to at level m. Indices of level 1
%			are 1:n.
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

if nargin < 2 || isempty(data) || isempty(distanceFunction)
    error('groupData needs at least two non-empty input arguments!')
end

if isempty(which(distanceFunction))
    error('The distance function %s needs to be a function on the matlab path',distanceFunction)
end

if nargin < 3
    distanceFunctionParameters = [];
end

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
if isempty(associatedInfo)
    addInfo = false;
    linkData = [];
else
    addInfo = true;
    if size(associatedInfo,1) ~= nDistances
        error('associatedInfo has to have as many rows as there are distances')
    end
    linkData = zeros(nLinks, size(associatedInfo,2));
end

%==========================



%==========================
%% GROUPING LOOP
%==========================

% loop until all sets have been grouped, store info for all the linkages.
% Don't update distances in the last round

parameters.distanceOrder = (1:nSets)'; % create here because we never know
nGroups = nSets;


for i = 1:nLinks

    % create list of indices
    [q,p] = find(tril(true(nGroups),-1));
    distanceIdx = parameters.distanceOrder([p,q]);

    % group. Find minimum distance
    [minDist, minDistIdx] = min(distance);

    % store minimum distance, group
    newGroup = distanceIdx(minDistIdx,:);
    links(i,:) = [newGroup,minDist];
    if addInfo
        linkData(i,:) = associatedInfo(minDistIdx,:);
    end

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
        [distance, associatedInfo, data, parameters] = ...
    feval(distanceFunction, data, parameters);
        

    end

end % linking loop

% assign output
groupIdx = parameters.groupIdx;

