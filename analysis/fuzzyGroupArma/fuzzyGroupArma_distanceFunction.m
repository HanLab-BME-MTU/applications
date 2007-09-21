function distances = fuzzyGroupArma_distanceFunction(data,centers,distanceFunctionParameters)
%FUZZYGROUPARMA_DISTANCEFUNCTION is the distance function for possibilistic clustering in fuzzyGroupArma
%
% SYNOPSIS: distances = fuzzyGroupArma_distanceFunction(data,centers,distanceFunctionParameters)
%
% INPUT   data: data in the same form as output from armaxFitKalman, plus
%               the additional field orderLen
%         centers : ARMA descriptors of the cluster prototype
%         distanceFunctionParameters : placeholder
%
% OUTPUT  distances : nData-by-nCenters array of -log10(p) for the
%                     comparison between individual data items and the
%                     centers.
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 08-Nov-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% count input
nData = length(data);
nCenters = length(centers);

% ready output
distances = zeros(nData,nCenters);

% double-loop
for iCenter = 1:nCenters
    for iData = 1:nData
        % use armaxModelComp for p-values to ensure consistency
        distances(iData,iCenter) = armaxModelComp(data(iData),...
            centers(iCenter));
    end
end

