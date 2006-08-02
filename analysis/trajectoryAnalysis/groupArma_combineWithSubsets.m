function [data,recalc] = groupArma_combineWithSubsets(data,parameters,ijk)
%GROUPARMA_COMBINEWITHSUBSETS creates a new set of lengthSeries by picking only a subset of the data
%
% SYNOPSIS: [data,recalc] = groupArma_combineWithSubsets(data,parameters,ijk);
%
% INPUT data: input data to groupArmaDescriptors (this function will add a
%               field lengthSeriesInfo, with [dataSetIdx, lengthSeriesIdx]
%               showing which lengthSeries have been chosen for this group
%		parameters: input to grouping costfunctions, ijk: groups that are
%               to be combined (1:2), new group index (3)
%
% OUTPUT data: data with data(ijk(3)).lengthSeries if recalc is true
%        recalc: true or false dependent on whether to call recalcArma
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 31-Jul-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check whether we combine at all
if parameters.mode(2) == 0 || parameters.associatedInfo(2) < -log10(parameters.mode(2))
    % combine (mode2==0 means that there is no cutoff)
    recalc = true;

    % find data sets that are about to get grouped
    lastIdx = find(parameters.groupIdx(1,:),1,'last');
    groupedSets = parameters.groupIdx(...
        parameters.groupIdx(:,lastIdx) == ijk(3));
    nGroupedSets = length(groupedSets);

    % find number of observations in each of the sets
    totalNumObserve = cat(1,data(groupedSets).numObserve);

    % attempted number of observations is the average number of
    % observations from each set (or something else - we can
    % add to parameters.mode)
    goalNumObserveTotal = mean(totalNumObserve);

    % initialize the new length series
    data(ijk(3)).lengthSeries = [];
    data(ijk(3)).lengthSeriesInfo = [];

    % for every set: choose random trajectories, until we have
    % something that is close to the number of points we look
    % for. Try 5 times and select the best fit.
    for iSet = 1:nGroupedSets
        % find how many tiempoints we need. mode3==1 means to
        % use (roughly) the same number for all, mode3==2 means
        % to use proportional to numObserve.
        % There is never the risk of asking for too many
        % timepoints - we anyway attempt to get the number that
        % comes closest to target.
        if length(parameters.mode > 2) && parameters.mode(3) == 2
            goalNumObserve = goalNumObserveTotal / nGroupedSets;
        else
            % use default mode3=1
            goalNumObserve = goalNumObserveTotal * ...
                totalNumObserve(iSet)/sum(totalNumObserve);
        end

        % numObserveIndividual only has one col normally, but
        % as I made a mistake debugging, I have to select the
        % first col here.
        numObserveIndividual = data(groupedSets(iSet)).numObserveIndividual(:,1);
        nTraj = length(data(groupedSets(iSet)).lengthSeries);
        randomIdx = zeros(length(numObserveIndividual),5);
        for i=1:5
            randomIdx(:,i) = randperm(nTraj);
        end
        nTimepoints = numObserveIndividual(randomIdx);
        % cumsum to add up the data points
        nTimepoints = cumsum(nTimepoints,1);

        % subtract the goal and find the minimum (more
        % convenient than going through sub2ind)
        nTimepoints = abs(nTimepoints - goalNumObserve);
        minDelta = min(nTimepoints(:));
        % minDelta(1) in case there are multiple options
        [row,col] = find(nTimepoints == minDelta(1));

        % now we have the sets we want from this group
        data(ijk(3)).lengthSeries = cat(2,...
            data(ijk(3)).lengthSeries,...
            data(groupedSets(iSet)).lengthSeries(randomIdx(1:row,col)));

        % also store which of the series we have used
        data(ijk(3)).lengthSeriesInfo = cat(1,...
            data(ijk(3)).lengthSeriesInfo,...
            [groupedSets(iSet)*ones(row,1),randomIdx(1:row,col)]);
    end % loop sets
else
    recalc = false;
end % if combine

