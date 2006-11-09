function [data] = groupArma_recalcArma(data,ijk,guess,weights,name)
%GROUPARMA_RECALCARMA is a utility to recalculate ARMA descriptors
%
% SYNOPSIS: [data] = groupArma_recalcArma(data,ijk)
%
% INPUT data: structure array with ARMA-fitResults
%		ijk : 1-by-3 array with data sets (ij) that are grouped into a new
%             set (k). If ijk is a 1-by-2 cell array, the list of data sets
%             in the first element will be used to calculate the parameters
%             of the second element.
%       guess : (opt) Type of initial guess to be used {1}
%               1: ARMA descriptors (and model order) of largest data set
%               2: Try range of orders, use BIC to select best
%       weights : (opt) vector with the same number of elements as there
%             are data sets to be grouped. Used to calculate weighted ARMA
%             descriptors
%       name : (opt) string with name for the new time series
%
% OUTPUT data: updated data
%
% REMARKS If the data in ijk(3) already contains a length series, the
%         series aren't being updated!
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: jdorn
% DATE: 19-Jun-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%====================
%% CHECK INPUT
%====================

if nargin < 2 || isempty(data) || isempty(ijk)
    error('not enough input arguments')
end
if ~iscell(ijk)
    % transform to cell array
    ijk = {ijk(1:2),ijk(3)};
end
if nargin < 3 || isempty(guess)
    guess = 1;
end
% check for weights - set 1 if there aren't any
if nargin < 4 || isempty(weights)
    % don't do anything
else
    % add field .weight to every element in data
    for i=length(weights)
        data(ijk{1}(i)).weight = weights(i);
    end
end

% check for name. If none, make either a combination of the names of the
% two input series, or list the number of movies that have been combined
% here
if nargin < 5 || isempty(name)
    nCombine = length(ijk{1});
    if nCombine == 2
        name = sprintf('%s & %s',...
    data(ijk{1}(1)).name,data(ijk{1}(2)).name);
    else
        name = sprintf('Combination of %i length series',nCombine);
    end
end

%====================


%====================
%% RECALC ARMA
%====================

% make new length series if necessary
if length(data)<ijk{2} || isempty(data(ijk{2}).lengthSeries)
    data(ijk{2}).lengthSeries = ...
        cat(2,data(ijk{1}).lengthSeries);
end

% set initial guesses for new data. Take from set with most
% observations
switch guess
    case 1
numObserve = cat(1,data(ijk{1}).numObserve);
[dummy,selectIdx] = max(numObserve);
initialIdx = ijk{1}(selectIdx);
% ar, ma Param have transformed values in the first col
initialGuess.arParamP0 = data(initialIdx).arParamK(2,:);
initialGuess.maParamP0 = data(initialIdx).maParamK(2,:);
initialGuess.xParam0 = data(initialIdx).xParamK;
    case 2
        initialGuess = [0,3;0,3;-1,-1];
end

% recalculate
fitResults = armaxFitKalman(data(ijk{2}).lengthSeries,[],initialGuess);

% assign to data
if guess == 2
    % there are multiple models in fitResults. Choose best one
    success = cat(1,fitResults.success);
    bic = catStruct(1,'fitResults.selectCrit.bic');
    % remove non-successes
    bic(success == 0) = NaN;
    [dummy,idx] = nanmin(bic);
    fitResults = fitResults(idx);
end

for fn=fieldnames(fitResults)'
    data(ijk{2}).(char(fn)) = fitResults.(char(fn));
end

% store orderLen
if guess == 1
data(ijk{2}).orderLen = data(ijk{1}(selectIdx)).orderLen;
else
    data(ijk{2}).orderLen = [length(fitResults.arParamK) - 1,...
        length(fitResults.maParamK) - 1];
end
% store name
data(ijk{2}).name = name;