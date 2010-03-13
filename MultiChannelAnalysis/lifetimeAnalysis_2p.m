function output = lifetimeAnalysis_2p(data, Results, kWeibull, cutoffIdx, name)
%(Results, restrict, kWeibull, estartvec, efixvec)

% data = fillStructLifetimeHist(data);
% data = determinePitDensities(data);
% data = determineLifetimeInfo(data);
% 
% Results = stageFitLifetimes(data);


if (nargin < 3)
    cutoffIdx = 1;
end
if (nargin < 4)
    name = [];
else
    name = [' (' name ')'];
end

t = Results.hist_slow(1,:);
histVect = Results.hist_slow(2,:);
% cut-off for 2-population fit:
%idx = find(histVect==max(histVect));
offset = sum(histVect(1:cutoffIdx-1));
histVect(1:cutoffIdx-1) = [];

tc = t(cutoffIdx:end);

% optional input restrict: restrict final fitting analysis to a stretch of
% the data, e.g. to the first 300s of the histogram, since measured
% frequncies for lifetimes >300s are more than 50% speculative due to the
% correction for movie length
% if nargin>1
%     resT = restrict;
% else
%     resT = max(tvecslow);
% end

if nargin < 2
    kWeibull = [1 2 2];
end


opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e6, ...
    'MaxIter', 1e6, ...
    'Display', 'off', ...
    'TolX', 1e-12, ...
    'Tolfun', 1e-12);


%========================================================================
% 1. Fit histogram
%========================================================================
N = length(kWeibull);
aVect = ones(1,N)/N;
lambdaVect = (1:N)/(N+1)*t(floor(end/10));

% fit
initV = [lambdaVect aVect];
[p, resnorm] = lsqnonlin(@nWeibullCost, initV, [], [], opts, tc, histVect, kWeibull);
lambdaVect = abs(p(1:N));
aVect = abs(p(N+1:2*N));

% Schwarz criterion (assumption: normal errors)
n = length(histVect);
k = length(initV);
BIC = n*log(resnorm/n) + k*log(n);
fprintf('BIC = %.2f\n', BIC);

% Fitted curve
[w W] = nWeibull(t, kWeibull, lambdaVect, aVect);

% Plot histogram & fit
figure;
plot(tc, histVect, 'k-', 'LineWidth', 1.5);
hold on;
plot(tc, histVect, 'k.', 'MarkerSize', 10);
for k = 1:N
    plot(t, W(k,:), 'b', 'LineWidth', 2);
end
plot(t, w, 'r-', 'LineWidth', 3);
axis([0 t(end) 0 max(histVect)]);
set(gca, 'FontName', 'Helvetica', 'FontSize', 14, 'LineWidth', 2);
xlabel('t [s]', 'FontName', 'Helvetica', 'FontSize', 14);
ylabel('Relative frequency', 'FontName', 'Helvetica', 'FontSize', 14);
title(['Lifetime histogram' name], 'FontName', 'Helvetica', 'FontSize', 14);

%print('-depsc2', '-r300', 'histogramFit.eps');
output.populationContributions = aVect/sum(aVect);
output.tau = lambdaVect;
output.nCell = length(data);
output.nCCP = Results.numcells_slow;
%output.tau50 = jackknifed value 


% ========================================================================
% 2. Fit cumulative histogram
% =========================================================================

cumulativeHist = cumsum(histVect) + offset;

% Fit (retaining 'lambda' values from histogram fit)
p = lsqnonlin(@nWeibullCDFCost, [aVect 0], [], [], opts, tc, cumulativeHist, kWeibull, lambdaVect);
aVect = abs(p(1:N));
dt = p(N+1);

% Fitted curve
[w W] = nWeibullCDF(t-dt, kWeibull, lambdaVect, aVect);

% Plot cumulative histogram w/ fit
figure;
plot(tc, cumulativeHist, 'k-', 'LineWidth', 1);
hold on;
plot(tc, cumulativeHist, 'k.', 'Markersize', 10);
for k = 1:N
    plot(t, W(k,:), 'b', 'LineWidth', 1.5);
end
plot(t, w, 'r-', 'LineWidth', 2);
axis([0 t(end) 0 1]);
set(gca, 'YTick', 0:0.1:1, 'FontName', 'Helvetica', 'FontSize', 14, 'LineWidth', 1.5);
xlabel('t [s]', 'FontName', 'Helvetica', 'FontSize', 14);
ylabel('Relative frequency', 'FontName', 'Helvetica', 'FontSize', 14);
title(['Cumulative histogram' name], 'FontName', 'Helvetica', 'FontSize', 14);
%print('-depsc2', '-r300', 'cumulativeHistogramFit.eps');




function v = nWeibullCost(p, t, data, kVect)
N = length(kVect);
lambdaVect = abs(p(1:N));
aVect = abs(p(N+1:2*N));

v = data - nWeibull(t, kVect, lambdaVect, aVect);


function v = nWeibullCDFCost(p, t, data, kVect, lambdaVect)
N = length(kVect);
%lambdaVect = abs(p(1:N));
%aVect = abs(p(N+1:2*N));
aVect = abs(p(1:N));
dt = p(N+1);

v = data - nWeibullCDF(t-dt, kVect, lambdaVect, aVect);


function [w W] = nWeibull(t, k, lambda, A)
N = length(k);
nt = length(t);

K = repmat(k', [1 nt]);
L = repmat(lambda', [1 nt]);
t = repmat(t, [N 1])./L;

W = repmat(A', [1 nt]).*K./L.*t.^(K-1).*exp(-t.^K);
w = sum(W, 1);


function [w W] = nWeibullCDF(t, k, lambda, A)
N = length(k);
nt = length(t);

K = repmat(k', [1 nt]);
t = repmat(t, [N 1]) ./ repmat(lambda', [1 nt]);

W = repmat(A', [1 nt]) .* (1 - exp(-t.^K));
w = sum(W, 1);