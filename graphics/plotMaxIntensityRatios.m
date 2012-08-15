function plotMaxIntensityRatios(dataSets, varargin)

nSet = numel(dataSets);

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('CohortLB', [5  11 16 21 41 61]);
ip.addParamValue('CohortUB', [10 15 20 40 60 120]);
ip.addParamValue('Cutoff_f', 5);
ip.addParamValue('DisplayMode', '');
ip.addParamValue('FaceColor', []);
ip.addParamValue('EdgeColor', []);
ip.parse(varargin{:});
cfM = ip.Results.FaceColor;
cfE = ip.Results.EdgeColor;
if isempty(cfM)
    hv = (0:nSet-1)/max(nSet,3);
    cfM = hsv2rgb([hv(end:-1:1)' 0.3*ones(nSet,1) ones(nSet,1)]);
end
if isempty(cfE)
    hv = rgb2hsv(cfM);
    ceM = hsv2rgb([hv(:,1) ones(nSet,2)]);
end

lb = ip.Results.CohortLB;
ub = ip.Results.CohortUB;
nc = numel(lb);

% for each group of data, generate intensity ratio cohorts


pct = cell(1,nSet);
for s = 1:nSet
    data = dataSets{s};
    nd = numel(data);
    lftData = getLifetimeData(data);
    A = arrayfun(@(i) i.A(i.lifetime_s(i.catIdx==1)>=ip.Results.Cutoff_f,:), lftData, 'UniformOutput', false);
    maxA_all = cellfun(@(i) nanmax(i,[],2)', A, 'UniformOutput', false);
    a = rescaleEDFs(maxA_all, 'Display', false);
    % apply scaling
    for i = 1:nd
        maxA_all{i} = a(i) * maxA_all{i};
        A{i} = a(i)*A{i};
    end
    A = vertcat(A{:});
    lifetime_s = [lftData.lifetime_s];
    lifetime_s = lifetime_s([lftData.trackLengths]>=ip.Results.Cutoff_f & [lftData.catIdx]==1);
    
    ratio = nanmax(A,[],2) ./ A(:,1);
    
    pct{s} = zeros(5,nc);
    for c = 1:nc
        cidx = lb(c)<=lifetime_s & lifetime_s<ub(c);
        pct{s}(:,c) = prctile(ratio(cidx), [50 25 75 5 95]);
    end
end

%%
fset = loadFigureSettings(ip.Results.DisplayMode);

cohortLabels = arrayfun(@(i) [' ' num2str(lb(i)) '-' num2str(ub(i)) ' s'], 1:nc, 'UniformOutput', false);

figure(fset.fOpts{:});
axes(fset.axOpts{:});
hold on;
xlabel('Lifetime cohort', fset.lfont{:});
boxplot2(pct, ...
    'FaceColor', cfM, 'EdgeColor', ceM, 'ErrorbarColor', ceM,...
    'GroupDistance', 1, 'BarWidth', 0.8, 'LineWidth', 1,...
    'XLabels', cohortLabels,...
    'YLim', [0 6]);

ylabel('Max. / 1^{st} frame intensity', fset.lfont{:});



