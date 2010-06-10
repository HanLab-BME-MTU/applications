function output = lifetimeAnalysis(data, Results, kWeibull, cutoffIdx, name)
%(Results, restrict, kWeibull, estartvec, efixvec)

% data = fillStructLifetimeHist(data);
% data = determinePitDensities(data);
% data = determineLifetimeInfo(data);
% 
% Results = stageFitLifetimes(data);

if nargin < 2
    kWeibull = [1 2 2];
end
if (nargin < 3)
    cutoffIdx = 1;
end
if (nargin < 4)
    name = [];
else
    name = [' (' name ')'];
end

histName = 'hist_2s';


N = length(kWeibull);

t = Results.(histName)(1,:);
tc = t(cutoffIdx:end);
t = (0:size(Results.(histName),2))*data(1).framerate;

% optional input restrict: restrict final fitting analysis to a stretch of
% the data, e.g. to the first 300s of the histogram, since measured
% frequncies for lifetimes >300s are more than 50% speculative due to the
% correction for movie length
% if nargin>1
%     resT = restrict;
% else
%     resT = max(tvecslow);
% end

% (1:N)/(N+1)*t(floor(end/10))

histVect = Results.(histName)(2,:);
histVect = histVect / sum(histVect);
histVect(1:cutoffIdx-1) = [];
%histMean = sum(histVect.*tc);

% initial values
initVect = zeros(1,1+3*N);
initVect(2:3:end) = ones(1,N)/N; % A
% initVect(3:3:end) = 2.^(1:N)/2^N * histMean; % lambda
initVect(3:3:end) = 25 * (1:N);

initVect(4:3:end) = kWeibull;
estVect = [1 repmat([1 1 0], [1 N])];


% loop though individual histograms
for k = 1:size(Results.histMatrix_2s, 1)
    histVect = Results.histMatrix_2s(k,:);
    histVect = histVect/sum(histVect);
    offset = sum(histVect(1:cutoffIdx-1));
    histVect(1:cutoffIdx-1) = [];
    [prmVect] = fitHistogram(tc, histVect, offset, initVect, estVect);
    output.indPrm(k,:) = prmVect;
end

% loop through averaged histograms in input structure
nHist = size(Results.(histName),1)-1;
for k = 1:nHist

    % re-normalize histogram
    histVect = Results.(histName)(k+1,:);
    histVect = histVect / sum(histVect);
    offset = sum(histVect(1:cutoffIdx-1));
    histVect(1:cutoffIdx-1) = [];  

    [prmVect, histVect, cumulativeHist, ~, BIC] = fitHistogram(tc, histVect, offset, initVect, estVect);

    if k==1
        fprintf('BIC = %.2f %s\n', BIC, name);
        output.histVect = histVect;
        output.cumulativeHist = cumulativeHist;
    end

    output.prmVect(k,:) = prmVect;
end
prmVect = output.prmVect(1,:);

% Jackknife standard deviation
%sqrt((nHist-2)/(nHist-1)*sum((output.prmVect(2:end,:) - repmat(mean(output.prmVect(2:end,:)), [nHist-1 1])).^2,1))
JK = sqrt((nHist-2)*var(output.prmVect(2:end,:), 1, 1));

% Store results
output.populationContributions = prmVect(2:3:end) / sum(prmVect(2:3:end));
output.A_JKerror = JK(2:3:end);
output.A_c2c = mean(output.indPrm(:,2:3:end),1);
output.A_std_c2c = std(output.indPrm(:,2:3:end),1);


output.tau = prmVect(3:3:end);
output.tau_JKerror = JK(3:3:end);
output.tau_c2c = mean(output.indPrm(:,3:3:end),1);
output.tau_std_c2c = std(output.indPrm(:,3:3:end),1);

output.median = output.tau .* nthroot(-log(0.5), kWeibull);
output.percentile25 = output.tau .* nthroot(-log(0.75), kWeibull);
output.percentile75 = output.tau .* nthroot(-log(0.25), kWeibull);
output.range50 = arrayfun(@(x) boxwhiskerPerRange(output.tau(x), kWeibull(x), 0.5), 1:N);

output.nCell = length(data);
output.nCCP = Results.numtracks_2s;


% Plot results

% Initial histogram fit
% figure;
% plot(tc, histVect, 'k.-', 'LineWidth', 1, 'MarkerSize', 10);
% hold on;
% plot(t, W, 'b', 'LineWidth', 1.5);
% plot(t, w, 'r-', 'LineWidth', 2);
% axis([0 t(end) 0 1.1*max(histVect)]);
% set(gca, 'FontName', 'Helvetica', 'FontSize', 14, 'LineWidth', 1.5);
% xlabel('t [s]', 'FontName', 'Helvetica', 'FontSize', 14);
% ylabel('Relative frequency', 'FontName', 'Helvetica', 'FontSize', 14);
% title(['Lifetime histogram' name], 'FontName', 'Helvetica', 'FontSize', 14);


% Plot cumulative histogram w/ fit
[w W] = nWeibull(t, prmVect, 'CDF');

figure;
plot(tc, output.cumulativeHist, 'k.-', 'LineWidth', 1, 'MarkerSize', 15);
hold on;
plot(t, w, 'r', 'LineWidth', 2);
plot(t, W, 'b', 'LineWidth', 1.5);
axis([0 t(end) 0 1]);
set(gca, 'YTick', 0:0.1:1, 'FontName', 'Helvetica', 'FontSize', 14, 'LineWidth', 1.5);
xlabel('t [s]', 'FontName', 'Helvetica', 'FontSize', 14);
ylabel('Relative frequency', 'FontName', 'Helvetica', 'FontSize', 14);
title(['Cumulative histogram' name], 'FontName', 'Helvetica', 'FontSize', 14);
%print('-depsc2', '-r300', 'CumulativeHistogramFit.eps');


% Plot CDF fit parameters on PDF
[w W] = nWeibull(t, prmVect, 'PDF');

figure;
plot(tc, output.histVect, 'k.-', 'LineWidth', 1, 'MarkerSize', 10);
hold on;
plot(t, w*data(1).framerate, 'r', 'LineWidth', 2);
plot(t, W*data(1).framerate, 'b', 'LineWidth', 1.5);
axis([0 t(end) 0 1.1*max(output.histVect)]);
set(gca, 'FontName', 'Helvetica', 'FontSize', 14, 'LineWidth', 1.5);
xlabel('t [s]', 'FontName', 'Helvetica', 'FontSize', 14);
ylabel('Relative frequency', 'FontName', 'Helvetica', 'FontSize', 14);
title(['Histogram' name], 'FontName', 'Helvetica', 'FontSize', 14);
%print('-depsc2', '-r300', 'HistogramFit.eps');



function [prmVect, histVect, cumulativeHist, offset, BIC] = fitHistogram(tc, histVect, offset, prmVect, estVect)
N = (length(prmVect)-1)/3;

% fit to histogram to refine parameter estimates
prmVect = fitNWeibull(tc, histVect, prmVect, estVect, 'PDF');

cumulativeHist = cumsum(histVect) + offset;

% fit to cumumative histogram
[prmVect, ~, ~, BIC] = fitNWeibull(tc, cumulativeHist, prmVect, estVect, 'CDF');

% renormalize histogram and correct offset
histVect = histVect / (1-prmVect(1));
offset = (offset-prmVect(1)) / (1-prmVect(1));
cumulativeHist = cumsum(histVect) + offset;

% renormalize amplitudes
prmVect(2:3:end) = prmVect(2:3:end)/(1-prmVect(1));
% correct offset
prmVect(1) = 0;

% sort parameter vector according to population mean
[~, order] = sort(prmVect(3:3:end));
idx = reshape(2:3*N+1, [3 N]);
idx = reshape(idx(:,order), [1 3*N]);
prmVect(2:end) = prmVect(idx);
