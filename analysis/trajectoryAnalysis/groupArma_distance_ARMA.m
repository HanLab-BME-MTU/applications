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
%              lengthSeries (if grouping by recalculating ARMA descriptors
%                            is selected)
%
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
%
%               cutoff - (passed down from groupArmaDescriptors)
%                        either a negative integer indicating how many
%                        groups to make, or a positive double with the
%                        cutoff probability
%               groupingMode - (passed down from groupArmaDescriptors)
%                        two-element array [mode, criterion]
%                        mode indicates which way the data should be
%                        grouped. 0: no recombinateion of data. 1:
%                        recombination up to a certain level of difference.
%                        2: recombination (=recalculating ARMA) all the way
%                        through.
%                        criterion: specifies the exact level until which
%                        the data should be grouped.
%
%               initialProb (created here) - -log10(p) of the first
%                           comparison
%               initialF (created here) - test statistics of first
%                           comparison
%               initialM (created here) - maxOrder, minDegFreedom for each
%                           pair in the first comparison
%                   initialF, initialM will be recalculated whenever new
%                   ARMA descriptors are generated
%               initialProbIdx (created here) - index into intialProb
%               initialFIdx (created here) - index into initial F/M
%               
%
% OUTPUT distance - vector containing the distance between the
%                       sets. Element 1 is the distance between 1 and
%                       2, element 2 between 1 and 3 etc., where 1,2,3,
%                       refer to the order of sets in data (newly joined
%                       sets will be added at the end)
%        associatedInfo - array with additional data for each pair in cols
%                       1: link incompatible according to cutoff? (required
%                           by groupData)
%                       2: -log10(prob)
%                       3: Fstat
%                       4: n1
%                       5: n2
%                       6: minLogProb of individual sets
%                       7: maxLogProb of individual sets
%                       8: 1 if extrapolated log
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
    
    % init recalc
    recalc = false;
    
    
else
    % we need to combine two sets. As all comparisons are based on the test
    % statistics of the initial round of comparisons, there is nothing to
    % do here as long as we don't recalculate ARMA descriptors
    
    ijk = parameters.group;

    % check for groupingMode
    switch parameters.mode(1)
        case 0
            % don't recalc
            recalc = false;
        case 1
            % recalc only if last p-value was below threshold
            if parameters.associatedInfo(2) < -log10(parameters.mode(2))
                recalc = true;
            else
                recalc = false;
            end
        case 2
            % recalc all the time
            recalc = true;
    end

    % if recalc, we're creating a new data set
    
    if recalc

        % make new length series
        data(ijk(3)).lengthSeries = ...
            [data(ijk(1)).lengthSeries,data(ijk(2)).lengthSeries];

        % set initial guesses for new data. Take from set with more
        % observations (if ijk2 has more, idx into ijk will become 2
        initialIdx = ijk((data(ijk(2)).numObserve > data(ijk(1)).numObserve)+1);
        % ar, ma Param have transformed values in the first col
        initialGuess.arParamP0 = data(initialIdx).arParamK(2,:);
        initialGuess.maParamP0 = data(initialIdx).maParamK(2,:);
        initialGuess.xParamP0 = data(initialIdx).xParamK;

        % recalculate
        fitResults = armaXFitKalman(data(ijk(3)).lengthSeries,[],initialGuess);

        % assign to data
        for fn=fieldnames(fitResults)'
            data(ijk(3)).(char(fn)) = fitResults.(char(fn));
        end

    end
end

%======================


%==========================
%% DistanceMatrix: prepare
%==========================

% double-loop through the data and compare. If the code becomes too slow
% here, remember all the old distances, and only update the necessary
% distances! 

nSets = length(parameters.distanceOrder);
nComparisons = (nSets-1)*nSets/2;
associatedInfo = zeros(nComparisons,8);
ct = 0;

if isInit || recalc
    % save the distanceIdx
    
    parameters.initialM = zeros(nComparisons,2);
    parameters.initialF = zeros(nComparisons,1);
    parametres.initialFIdx = zeros(nComparisons,2);
    
    
    % last recalc is what we'll use to make the recalc-groups
    
    
    if isInit
    parameters.initialProbIdx = zeros(nComparisons,2);
    parameters.initialProb = zeros(nComparisons,1);
    parameters.lastRecalc = 1;
    else
        parameters.lastRecalc = find(parameters.groupIdx(1,:),1,'last')+1;
    end
       
else
    % find the latest groupIdx
    lastIdx = find(parameters.groupIdx(1,:),1,'last');
    groupIdx = parameters.groupIdx(:,[1,parameters.lastRecalc,lastIdx]);
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
        
        if isInit || recalc
            % for every pair: get test statistic
            [data(iIdx).paramVec,data(jIdx).paramVec,...
                data(iIdx).varCovMat,data(jIdx).varCovMat] = ...
                armaxOrderMatch(data(iIdx),data(jIdx));
            
            fRatio = globalCoefTestStat(...
                data(iIdx).paramVec, data(jIdx).paramVec,...
                data(iIdx).varCovMat, data(jIdx).varCovMat);
          
            
            % fTest
            
            % get degrees of freedom. Use maximum for both to get strict
            % test
            n1 = max(sum(data(iIdx).orderLen),sum(data(jIdx).orderLen));
            n2 = max(data(iIdx).numObserve - sum(data(iIdx).orderLen),...
                data(jIdx).numObserve - sum(data(jIdx).orderLen));
            
            
            % directly calculate -log10(p)
            [logP,isExtrapolated] =...
                 fcdfExtrapolateLog(fRatio,n1,n2);
             
            % store associatedInfo
            associatedInfo(ct,[2:5,8]) = [logP, fRatio, n1, n2, isExtrapolated];
            
            
            % if init: store logP with associated idx. min/max is the same
            % as the actual probabilities
            if isInit
                
                % remember test 
            parameters.initialProbIdx(ct,:) = [i,j];
            parameters.initialProb(ct) = logP;
            
            % min/max probabilities are the same as the actual
            % probabilities if init
            associatedInfo(ct,2:3) = associatedInfo(ct,1);
            
            else
                % min/max have to be found from initialProb
                
                % find comparisons between elements in iSets and jSets. Take if
            % in a row of initialProbIdx there is an entery of iSets and
            % an entry of jSets
            comparisonIdx = ...
                any(ismember(parameters.initialProbIdx,iSets),2) &...
                any(ismember(parameters.initialProbIdx,jSets),2);
            
            % min/max
            associatedInfo(ct,6) = ...
                min(parameters.initialProb(comparisonIdx));
            associatedInfo(ct,7) = ...
                max(parameters.initialProb(comparisonIdx));
            
            end
            
            % store fRatios, n1, n2, corresponding index
            parameters.initialF(ct,1) = fRatio;
            parameters.initialM(ct,:) = [n1, n2];
            parameters.initialFIdx(ct,:) = [i,j];
            
        else
            
            % for every pair: check groupIdx, find the test statistics for
            % the comparison of every single test, then take the weighted
            % mean and calculate the probability from the combined
            % statistic
            
            % find sets (original data sets that belong to groups i and j
            iSets = groupIdx(groupIdx(:,3) == iIdx,2);
            jSets = groupIdx(groupIdx(:,3) == jIdx,2);
            
            % find comparisons between elements in iSets and jSets. Take if
            % in a row of initialProbIdx there is an entery of iSets and
            % an entry of jSets
            comparisonIdx = ...
                any(ismember(parameters.initialFIdx,iSets),2) &...
                any(ismember(parameters.initialFIdx,jSets),2);
            
            % weighted mean
            fRatio = sum(...
                parameters.initialF(comparisonIdx).*...
                parameters.initialM(comparisonIdx,1))/...
                sum(parameters.initialM(comparisonIdx,1));
            n1= max(parameters.initialM(comparisonIdx,1));
            n2 = max(parameters.initialM(comparisonIdx,2));
            % test. Use maximum (not sum) for degrees of freedom -
            % otherwise, we are favoring single-to-group linkage over
            % group-to-group
            [logP,isExtrapolated] =...
                fcdfExtrapolateLog(fRatio,n1,n2);
            
            % store associatedInfo
            associatedInfo(ct,[2:5,8]) = [logP, fRatio, n1, n2, isExtrapolated];
            
            % min/max
            % find sets (original data sets that belong to groups i and j
            iSets = groupIdx(groupIdx(:,3) == iIdx,1);
            jSets = groupIdx(groupIdx(:,3) == jIdx,1);
            comparisonIdx = ...
                any(ismember(parameters.initialProbIdx,iSets),2) &...
                any(ismember(parameters.initialProbIdx,jSets),2);
            associatedInfo(ct,6) = ...
                min(parameters.initialProb(comparisonIdx));
            associatedInfo(ct,7) = ...
                max(parameters.initialProb(comparisonIdx));
    
            
        end % if isInit
        
    end % loop j
end % loop i


%=======================
%% finish up
%=======================

%   * there used to be alternatives to -log(p) here. Go back to older
%   revisions to check them out if needed *


% read distance
distance = associatedInfo(:,2);

% check whether compatible pairs
if ~isempty(parameters.cutoff)
    % check whether we cutoff according to probability or according to
    % number of groups
    if parameters.cutoff == -abs(round(parameters.cutoff))
        % negative integer: cut number of groups
        if nSets <= -parameters.cutoff
            % all are different if there are only -cutoff sets around
            associatedInfo(:,1) = 1;
        end
    else
        associatedInfo(:,1) = associatedInfo(:,2) > -log10(parameters.cutoff);
    end
end

% % store initialProb
% if isInit
%     parameters.initialProb = associatedInfo(:,1);
% end

