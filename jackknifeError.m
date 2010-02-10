function [sigmaJN] = jackknifeError(datavector)
% jackknife error computation from data
% INPUT:    data =  assumed to be samples of parameter estimation with
%                   succesive data sets deleted
% OUTPUT:   sigmaJN = jackknife error


% number of samples
ns = length(datavector);

% jackknife variance
varJN = ((ns-1)/ns)*sum( (datavector-nanmean(datavector)).^2 );

% jackknife error
sigmaJN = sqrt(varJN);