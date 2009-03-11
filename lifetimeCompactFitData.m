function [compactRes, data] = lifetimeCompactFitData(data, restrict, shape);
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

[data] = fillStructLifetimeHist(data);

% k= number of movies/cells
k = length(data);

% combine lifetime histograms of the slow and fast movies separately
res     = stageFitLifetimesPlat(data);
res_E2E = stageFitLifetimes_E2E(data);
mer     = mergeFastSlowHistogramsPlat(res, restrict, shape);
mer_E2E = mergeFastSlowHistogramsPlat_E2E(res_E2E, restrict, shape);

% n = number of trajectories
n = mer.numcells;

% fit results
resmat1 = mer.compactFitRes;

% determine cell-to-cell error with jackknife
for s = 1:length(mer_E2E)
    E2Emat(:,:,s) = mer_E2E(s).compactFitRes;
end
[mx,my,mz] = size(E2Emat);
for inx=1:mx
    for iny=1:my
        E2Eval(inx,iny) = jackknifeError(E2Emat(inx,iny,:));
    end
end

% enter output results
compactRes.numcells = k;
compactRes.numtraj = n;
compactRes.contr = resmat1(:,1);
compactRes.contrError = round(10*E2Eval(:,1))/10;
compactRes.tau = resmat1(:,2);
compactRes.tauError = round(10*E2Eval(:,2))/10;
compactRes.tau50 = resmat1(:,3);

 [data] = determinePitDensities(data);
 for i=1:length(data)
     Dens(i) = nanmean(data(i).pitDensity);
 end
 densAV = nanmean(Dens);
 densSTD = nanstd(Dens);


matrix(:,1) = compactRes.contr;
matrix(:,2) = compactRes.contrError;
matrix(:,3) = compactRes.tau;
matrix(:,4) = compactRes.tauError;
matrix(:,5) = compactRes.tau50;

matrix(1,6) = compactRes.numcells;
matrix(2,6) = compactRes.numtraj;
matrix(3,6) = 1000*densAV;
matrix(4,6) = 1000*densSTD;





compactRes.matrix = matrix;


end
