function [sigmaJN]=jackknifeError(datavector);
% jackknife error computation from data
% INPUT:    data =  assumed to be samples of parameter estimation with
%                   succesive data sets deleted
% OUTPUT:   sigmaJN = jackknife error


% number of samples
ns = length(datavector);

% average parameter value from all jackknife estimations
pav = nanmean(datavector);

% sum of variances
vsum = sum( (datavector-pav).^2 );

% jackknife variance
varJN = ((ns-1)/ns)*vsum;

% jackknife error
sigmaJN = sqrt(varJN);



end % of function


