function [compactRes, data] = lifetimeCompactFitData(data, restrict, shape)
% perform multiple-population fit on lifetime data, and write all relevant
% results, including cell-to-cell variation and statistics, into a copmact
% results file
% INPUT:    data    = experiment data, has to contain fields
%           restrict    = time restriction in seconds, e.g. 300
%           shape   = shape for populations, e.g. [2 2 1], where 1
%                   indicates an exponential distribution, and 2 a Rayleigh
%                   distribution
% OUTPUT:   compactRes    = compact results
%           data (optional) = with added field .lftHist if not already
%           present
%
% Dinah Loerke, last modified Mar 13, 2008
% Francois Aguet, last modified Jun 8 2010

data = fillStructLifetimeHist(data);

% combine lifetime histograms of the slow and fast movies separately
res = stageFitLifetimes(data);

if length(unique([data.framerate])) == 1
    mer_E2E = fitWeibullToHist(res, restrict, shape);
else
    mer_E2E = mergeHistogramsFitWeibull(res, restrict, shape);
end

% fit results
resmat1 = mer_E2E(1).compactFitRes;

% determine cell-to-cell error with jackknife
ns = length(mer_E2E);
E2Emat = cat(3, mer_E2E.compactFitRes);
E2Eval = sqrt((ns-1)/ns*sum((E2Emat - repmat(nanmean(E2Emat, 3), [1 1 ns])).^2, 3));



% enter output results
compactRes.numcells = length(data);
compactRes.numtraj = mer_E2E(1).numcells;
compactRes.contr = resmat1(:,1);
compactRes.contrError = round(100*E2Eval(:,1))/100;
compactRes.tau = resmat1(:,2);
compactRes.tauError = round(100*E2Eval(:,2))/100;
compactRes.tau50 = resmat1(:,3);

Dens = arrayfun(@(x) nanmean(x.pitDensity), data);

compactRes.contr;
compactRes.tau;

matrix(:,1) = compactRes.contr;
matrix(:,2) = compactRes.contrError;
matrix(:,3) = compactRes.tau;
matrix(:,4) = compactRes.tauError;
matrix(:,5) = compactRes.tau50;

matrix(1,6) = compactRes.numcells;
matrix(2,6) = compactRes.numtraj;
matrix(3,6) = 1000*nanmean(Dens);
matrix(4,6) = 1000*nanstd(Dens);

compactRes.matrix = matrix;