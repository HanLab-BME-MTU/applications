function centers = fuzzyGroupArma_centerFunction(data,membership,typicality,m,centerFunctionParameters)
%FUZZYGROUPARMA_CENTERFUNCTION calculates new cluster prototypes from data and memberships
%
% SYNOPSIS: fuzzyGroupArma_centerFunction
%
% INPUT   data: data in the same form as output from armaxFitKalman, plus
%               the additional field orderLen
%         membership : nData-by-nCenters array of memberships
%         typicality : nData-by-nCenters array of typicalities
%         m : fuzzyness-exponent
%         centerFunctionParameters.names : center-names
%                                 .recalc : ARMA recalc-option
%
% OUTPUT  centers : ARMA descriptors of cluster prototypes
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 09-Nov-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% find number of data points, number of centers
[nData, nCenters] = size(membership);
% init centers with data. Initialize field fitResults
data(1).fitResults = [];
[centers(1:nCenters)] = deal(data(1));

% assign empty names if no names given
if isempty(centerFunctionParameters.names)
    centerFunctionParameters.names = cell(nCenters,1);
end

% m=1

% to remove the influence of the trajectory length on the estimation,
% divide the weights by the length of the trajectories
numObserve = cat(1,data.numObserve);
 
% for every center: recalc ARMA
for iCenter = 1:nCenters
    % -> check whether we may have to put membership to the mth power!
    [data, fitResults] = groupArma_recalcArma(data,{1:nData,nData+1},...
        centerFunctionParameters.recalc,...
        (membership(:,iCenter).*typicality(:,iCenter)).^m./numObserve,centerFunctionParameters.names{iCenter});

    % data(nData+1) is the new center. Store and remove
    centers(iCenter) = data(end);
    centers(iCenter).fitResults = fitResults;
    data(end) = [];

end

