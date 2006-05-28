function [distance, associatedInfo, data, parameters] = groupArma_distance_ARMA(data,parameters)
%GROUPARMA_DISTANCE_WNV is the distance function for grouping the white noise variance when comparing Arma models
%
% SYNOPSIS: [distance, associatedInfo, data] = groupArma_distance_ARMA(data,parameters)
%
% INPUT data - input data to groupData with at least the fields
%              xParamK
%              arParamK
%              maParamK
%              varCovarMatF
%              numObserve - nuber of observations
%              orderLen - [ar,ma] order length
%       parameters - structure with fields
%               group (required by groupData) - either [], if the distance
%                           function has to be 
%                           calculated for the entire data set, or [i,j,k],
%                           if sets i and j are grouped into a new set, k.
%               distanceOrder (req. by groupData) - ordered indices of
%                           active groups. If 
%                           the last two of six data sets have been
%                           grouped (group = [5,6,7]), distanceOrder will
%                           be [1,2,3,4,7]. If group is empty, then
%                           distanceOrder will be empty, too!
%               groupIdx (req. by groupData) 
%               initialProb (created here) - probabilities of the first
%                           comparison (1-pValue)
%               initialF (created here) - test statistics of first
%                           comparison
%               initialM (created here) - maxOrder, minDegFreedom for each
%                           pair in the first comparison
%               initialProbIdx (created here) - index into intialProb
%               mode (input) - distance that will be used to group data
%                       1  1-pValue (high p-value=short distance)
%                       2  -log(pValue)
%                       3  fcdf(1-pValue, 2,2000)
%                       4  testStatistic
%               
%
% OUTPUT distance - vector containing the distance between the
%                       sets. Element 1 is the distance between 1 and
%                       2, element 2 between 1 and 3 etc., where 1,2,3,
%                       refer to the order of sets in data (newly joined
%                       sets will be added at the end)
%        associatedInfo - array with additional data for each pair in cols
%                       1: prob
%                       2: logProb
%                       3: fStat
%                       4: test statistic
%                       5: minProb of individual sets
%                       6: maxProb of individual sets
%                       7: minLogProb of individual sets
%                       8: maxLogProb of individual sets
%                       9: minFStat of individual sets
%                      10: maxFStat of individual sets
%                      11: minDirectF of individual sets
%                      12: maxDirectF of individual sets
%        parameters - updated parameter structure. 
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 22-May-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first, combine if necessary, or write distanceOrder and set flag to store
% distances.
% Then, loop through all pairs and calculate the distance measures.

%=============
%% test Input
%=============
isInit = ~isfield(parameters,'group') || isempty(parameters.group);

if isInit
    % initialize
    nData = length(data);
    parameters.distanceOrder = 1:nData;
    
    % make sure all the data has a parameter vector of correct length
    
    
else
    % we need to combine two sets. As all comparisons are based on the test
    % statistics of the initial round of comparisons, there is nothing to
    % do here    
end

% turn off log0-warning
w=warning;
warning('off','MATLAB:log:logOfZero');
%======================


%==========================
%% DistanceMatrix: prepare
%==========================

% double-loop through the data and compare. If the code becomes too slow
% here, remember all the old distances, and only update the necessary
% distances! 

nSets = length(parameters.distanceOrder);
nComparisons = (nSets-1)*nSets/2;
associatedInfo = zeros(nComparisons,12);
ct = 0;

if isInit
    % save the distanceIdx
    parameters.initialProbIdx = zeros(nComparisons,2);
    parameters.initialM = zeros(nComparisons,2);
       
else
    % find the latest groupIdx
    lastIdx = find(parameters.groupIdx(1,:),1,'last');
    groupIdx = parameters.groupIdx(:,[1,lastIdx]);
end

%============================
%% DistanceMatrix: calculate
%============================

% in the first round, we need to calculate all the individual comparisons.
% later, we just look up the initial test statistics
for i = 1:nSets-1
    iIdx = parameters.distanceOrder(i);
    for j = i+1:nSets
        jIdx = parameters.distanceOrder(j);
        ct = ct + 1;
        
        if isInit
            % for every pair: get test statistic
            [data(iIdx).paramVec,data(jIdx).paramVec,...
                data(iIdx).varCovMat,data(jIdx).varCovMat] = ...
                armaxOrderMatch(data(iIdx),data(jIdx));
            associatedInfo(ct,4) = globalCoefTestStat(...
                data(iIdx).paramVec, data(jIdx).paramVec,...
                data(iIdx).varCovMat, data(jIdx).varCovMat);
            
            % fTest
            parameters.initialM(ct,1) = ...
                max(sum(data(iIdx).orderLen),sum(data(jIdx).orderLen));
            parameters.initialM(ct,2) = ...
                min(data(iIdx).numObserve - sum(data(iIdx).orderLen),...
                data(jIdx).numObserve - sum(data(jIdx).orderLen));
            
            associatedInfo(ct,1) = fcdf(associatedInfo(ct,4),...
                parameters.initialM(ct,1),parameters.initialM(ct,2));
            
            % min/max probabilities are the same as the actual
            % probabilities
            associatedInfo(ct,5:6) = associatedInfo(ct,1);
            associatedInfo(ct,11:12) = associatedInfo(ct,4);
            
            % do -log, fStats below.
            
            % remember test 
            parameters.initialProbIdx(ct,:) = [i,j];
            
        else
            
            % for every pair: check groupIdx, find the test statistics for
            % the comparison of every single test, then take the weighted
            % mean and calculate the probability from the combined
            % statistic
            
            % find sets (original data sets that belong to groups i and j
            iSets = groupIdx(groupIdx(:,2) == iIdx,1);
            jSets = groupIdx(groupIdx(:,2) == jIdx,1);
            
            % find comparisons between elements in iSets and jSets. Take if
            % in a row of initialProbIdx there is an entery of iSets and
            % an entry of jSets
            comparisonIdx = ...
                any(ismember(parameters.initialProbIdx,iSets),2) &...
                any(ismember(parameters.initialProbIdx,jSets),2);
            
            % weighted mean
            associatedInfo(ct,4) = sum(...
                parameters.initialF(comparisonIdx).*...
                parameters.initialM(comparisonIdx,1))/...
                sum(parameters.initialM(comparisonIdx,1));
            
            % test
            associatedInfo(ct,1) = fcdf(associatedInfo(ct,4),...
                parameters.initialM(ct,1),parameters.initialM(ct,2));
            
            % min/max
            associatedInfo(ct,5) = ...
                min(parameters.initialProb(comparisonIdx));
            associatedInfo(ct,6) = ...
                max(parameters.initialProb(comparisonIdx));
            
            associatedInfo(ct,11) = ...
                min(parameters.initialF(comparisonIdx));
            associatedInfo(ct,12) = ...
                max(parameters.initialF(comparisonIdx));
            
        end % if isInit
        
    end % loop j
end % loop i


%=======================
%% finish up
%=======================

% fill up associatedInfo, fill in distance, and, if necessary, initialProb

% -logProb. pValue is already 1-p
associatedInfo(:,[2,7,8]) = -log10(1-associatedInfo(:,[1,5,6]));

% fStats
associatedInfo(:,[3,9,10]) = finv(associatedInfo(:,[1,5,6]),2,2000);

% read distance
distance = associatedInfo(:,parameters.mode);

% store initialProb
if isInit
    parameters.initialProb = associatedInfo(:,1);
    parameters.initialF = associatedInfo(:,4);
end

% reset warnings
warning(w);