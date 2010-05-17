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
N = length(kWeibull);

t = Results.hist_slow(1,:);
tc = t(cutoffIdx:end);
t = (0:size(Results.hist_slow,2))*data(1).framerate;

% optional input restrict: restrict final fitting analysis to a stretch of
% the data, e.g. to the first 300s of the histogram, since measured
% frequncies for lifetimes >300s are more than 50% speculative due to the
% correction for movie length
% if nargin>1
%     resT = restrict;
% else
%     resT = max(tvecslow);
% end



% loop through all histograms in input structure
nHist = size(Results.hist_slow,1)-1;

for k = 1:nHist

    % re-normalize histogram
    histVect = Results.hist_slow(k+1,:);
    histVect = histVect / sum(histVect);
    
    % cut-off for 2-population fit:
    %idx = find(histVect==max(histVect));
    offset = sum(histVect(1:cutoffIdx-1));
    histVect(1:cutoffIdx-1) = [];
    
    
    % initial values
    prmVect = zeros(1,1+3*N);
    prmVect(2:3:end) = ones(1,N)/N; % A
    prmVect(3:3:end) = (1:N)/(N+1)*t(floor(end/10)); % lambda
    prmVect(4:3:end) = kWeibull;
    estVect = [1 repmat([1 1 0], [1 N])];
    
    % fit to histogram to refine parameter estimates
    [prmVect, residuals, estimatesSigma, BIC] = fitNWeibull(tc, histVect, prmVect, estVect, 'PDF');
    [w W] = nWeibull(t, prmVect, 'PDF');
    
    
    cumulativeHist = cumsum(histVect) + offset;
    
    % fix lambda values from histogram fit
    % estVect(3:3:end) = 0;
    
    % fit to cumumative histogram
    [prmVect, residuals, estimatesSigma, BIC] = fitNWeibull(tc, cumulativeHist, prmVect, estVect, 'CDF');
    if k==1
        fprintf('BIC = %.2f %s\n', BIC, name);
    end
    
    % renormalize histogram and correct offset
    histVect = histVect / (1-prmVect(1));
    offset = (offset-prmVect(1)) / (1-prmVect(1));
    cumulativeHist = cumsum(histVect) + offset;
    
    % renormalize amplitudes
    prmVect(2:3:end) = prmVect(2:3:end)/(1-prmVect(1));
    % correct offset
    prmVect(1) = 0;
    
    % sort parameter vector according to population mean
    [dummy order] = sort(prmVect(3:3:end));
    idx = reshape(2:3*N+1, [3 N]);
    idx = reshape(idx(:,order), [1 3*N]);
    prmVect(2:end) = prmVect(idx);
    output.prmVect(k,:) = prmVect;
end


% Jackknife standard deviation
%sqrt((nHist-2)/(nHist-1)*sum((output.prmVect(2:end,:) - repmat(mean(output.prmVect(2:end,:)), [nHist-1 1])).^2,1))
JK = sqrt((nHist-2)*var(output.prmVect(2:end,:), 1, 1));

% Store results
output.populationContributions = prmVect(2:3:end) / sum(prmVect(2:3:end));
output.A_JKerror = JK(2:3:end);

output.tau = prmVect(3:3:end);
output.tau_JKerror = JK(3:3:end);

output.median = output.tau .* nthroot(-log(0.5), kWeibull);
output.percentile25 = output.tau .* nthroot(-log(0.75), kWeibull);
output.percentile75 = output.tau .* nthroot(-log(0.25), kWeibull);
output.range50 = arrayfun(@(x) boxwhiskerPerRange(output.tau(x), kWeibull(x), 0.5), 1:N);

output.nCell = length(data);
output.nCCP = Results.numcells_slow;


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
plot(tc, cumulativeHist, 'k.-', 'LineWidth', 1, 'MarkerSize', 15);
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
plot(tc, histVect, 'k.-', 'LineWidth', 1, 'MarkerSize', 10);
hold on;
plot(t, w*data(1).framerate, 'r', 'LineWidth', 2);
plot(t, W*data(1).framerate, 'b', 'LineWidth', 1.5);
axis([0 t(end) 0 1.1*max(histVect)]);
set(gca, 'FontName', 'Helvetica', 'FontSize', 14, 'LineWidth', 1.5);
xlabel('t [s]', 'FontName', 'Helvetica', 'FontSize', 14);
ylabel('Relative frequency', 'FontName', 'Helvetica', 'FontSize', 14);
title(['Histogram' name], 'FontName', 'Helvetica', 'FontSize', 14);
%print('-depsc2', '-r300', 'HistogramFit.eps');

