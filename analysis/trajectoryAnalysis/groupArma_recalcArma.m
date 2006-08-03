function [data] = groupArma_recalcArma(data,ijk)
%GROUPARMA_RECALCARMA is a utility to recalculate ARMA descriptors
%
% SYNOPSIS: [data] = groupArma_recalcArma(data,ijk)
%
% INPUT data: structure array with ARMA-fitResults
%		ijk data sets (ij) that are grouped into a new set (k)
%
% OUTPUT data: updated data
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: jdorn
% DATE: 19-Jun-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make new length series if necessary
if length(data)<ijk(3) || isempty(data(ijk(3)).lengthSeries)
data(ijk(3)).lengthSeries = ...
    [data(ijk(1)).lengthSeries,data(ijk(2)).lengthSeries];
end

% set initial guesses for new data. Take from set with more
% observations (if ijk2 has more, idx into ijk will become 2
selectIdx = (data(ijk(2)).numObserve > data(ijk(1)).numObserve)+1;
initialIdx = ijk(selectIdx);
% ar, ma Param have transformed values in the first col
initialGuess.arParamP0 = data(initialIdx).arParamK(2,:);
initialGuess.maParamP0 = data(initialIdx).maParamK(2,:);
initialGuess.xParam0 = data(initialIdx).xParamK;

% recalculate
fitResults = armaxFitKalman(data(ijk(3)).lengthSeries,[],initialGuess);

% assign to data
for fn=fieldnames(fitResults)'
    data(ijk(3)).(char(fn)) = fitResults.(char(fn));
end

% store orderLen
data(ijk(3)).orderLen = data(ijk(selectIdx)).orderLen;
% store name
data(ijk(3)).name = sprintf('%s & %s',...
    data(ijk(1)).name,data(ijk(2)).name);