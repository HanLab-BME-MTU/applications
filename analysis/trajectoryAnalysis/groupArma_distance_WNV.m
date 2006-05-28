function [distance, associatedInfo, data, parameters] = groupArma_distance_WNV(data,parameters)
%GROUPARMA_DISTANCE_WNV is the distance function for grouping the white noise variance when comparing Arma models
%
% SYNOPSIS: [distance, associatedInfo, data] = groupARma_distance_WNV(data,group,distanceOrder, options)
%
% INPUT data - input data to groupData with at least the fields
%              numObserve - nuber of observations
%              wnVariance - white noise variance
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
%               initialF (created here) - variance ratios of the first
%                           comparison
%               initialProbIdx (created here) - index into intialProb
%               mode (input) - distance that will be used to group data
%                       1  1-pValue (high p-value=short distance)
%                       2  -log(pValue)
%                       3  fcdf(1-pValue, 2000,2000)
%                       4  variance ratios
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
%                       4: variance ratio
%                       5: minProb of individual sets
%                       6: maxProb of individual sets
%                       7: minLogProb of individual sets
%                       8: maxLogProb of individual sets
%                       9: minFStat of individual sets
%                      10: maxFStat of individual sets
%                      11: minDirectF of individual sets
%                      12: maxDirectF of individual sets
%        data - the (updated) input data with appended variances and n's
%        parameters - parameter structure. 
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
    
    % augment data now
    data(2*nData-1).wnVariance = [];
    
else
    % we need to combine two sets
    ijk = parameters.group;
    
    % the new variance is the weighted mean of the two individual variances
    % n: number of observations
    % m: number of fixed parameters (vector in the original data, but from
    % then on, all we use is the sum!)
    % nu: number of degrees of freedom
    n(1) = data(ijk(1)).numObserve;
    n(2) = data(ijk(2)).numObserve;
    m(1) = sum(data(ijk(1)).orderLen);
    m(2) = sum(data(ijk(2)).orderLen);
    nu(1) = n(1) - m(1);
    nu(2) = n(2) - m(2);
    data(ijk(3)).wnVariance = nu(1)/sum(nu) * data(ijk(1)).wnVariance + ...
        nu(2)/sum(nu) * data(ijk(2)).wnVariance;
    data(ijk(3)).numObserve = sum(n);
    data(ijk(3)).orderLen = sum(m);
    
end

% turn off log0-warning
w=warning;
warning('off','MATLAB:log:logOfZero');
%======================


%======================
%% DistanceMatrix
%======================

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
    
else
    % find the latest groupIdx
    lastIdx = find(parameters.groupIdx(1,:),1,'last');
    groupIdx = parameters.groupIdx(:,[1,lastIdx]);
end

for i = 1:nSets-1
    iIdx = parameters.distanceOrder(i);
    for j = i+1:nSets
        jIdx = parameters.distanceOrder(j);
        ct = ct + 1;
        
        % calculate probability. We want the one where we divide the
        % larger variance over the smaller variance (slightly more
        % numerically stable than the other, we eventually need it
        % for the distance criterion). Since that's annoying    
        % because of the index-shuffling involved, we just divide one by
        % the other and then make sure that p>=0.5
        F = data(iIdx).wnVariance/data(jIdx).wnVariance;
        associatedInfo(ct,1) = fcdf(...
            F,...
            data(iIdx).numObserve - sum(data(iIdx).orderLen),...
            data(jIdx).numObserve - sum(data(jIdx).orderLen));
        if associatedInfo(ct,1) < 0.5
            associatedInfo(ct,1) = 1-associatedInfo(ct,1);
        end
        
        % don't calculate log, fStat yet - we do that in parallel at the
        % end
        
        % store variance ratio
        if F<1
            associatedInfo(ct,4) = 1/F;
        else
            associatedInfo(ct,4) = F;
        end
        
        
        % if not init, find distances between groups.
        if isInit
            % remember distanceIdx
            parameters.initialProbIdx(ct,:) = [i,j];
            
            % min/max probabilities are the same as the actual
            % probabilities
            associatedInfo(ct,5:6) = associatedInfo(ct,1);
            associatedInfo(ct,11:12) = associatedInfo(ct,4);
            
        else
            % read from initialProb
            
            % find sets (original data sets that belong to groups i and j
            iSets = groupIdx(groupIdx(:,2) == iIdx,1);
            jSets = groupIdx(groupIdx(:,2) == jIdx,1);
            
            % find comparisons between elements in iSets and jSets. Take if
            % in a row of initialProbIdx there is an entery of iSets and
            % an entry of jSets
            comparisonIdx = ...
                any(ismember(parameters.initialProbIdx,iSets),2) &...
                any(ismember(parameters.initialProbIdx,jSets),2);
            
            associatedInfo(ct,5) = ...
                min(parameters.initialProb(comparisonIdx));
            associatedInfo(ct,6) = ...
                max(parameters.initialProb(comparisonIdx));
            
            % do the same for variance ratios
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
associatedInfo(:,[3,9,10]) = finv(associatedInfo(:,[1,5,6]),2000,2000);

% read distance
distance = associatedInfo(:,parameters.mode);

% store initialProb
if isInit
    parameters.initialProb = associatedInfo(:,1);
    parameters.initialF = associatedInfo(:,4);
end

% reset warnings
warning(w);