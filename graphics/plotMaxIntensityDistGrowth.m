% plot max. intensity distribution for lifetime cohorts @ frame rate

function plotMaxIntensityDistGrowth(data)

lftData = getLifetimeData(data);
% cohorts
% lb = 5:120;
% ub = lb; 
lb = 5:5:120;
ub = lb+4;
nc = numel(lb);

maxA_all = arrayfun(@(i) nanmax(i.intMat_Ia,[],2), lftData, 'UniformOutput', false);

% Rescale EDFs (correction for FP-fusion expression level)
a = rescaleEDFs(maxA_all, 'Display', false);

% apply scaling
for i = 1:numel(maxA_all)
    lftData(i).intMat_Ia = a(i) * lftData(i).intMat_Ia;
    maxA_all{i} = a(i) * maxA_all{i};
end

% Concatenate data
maxA = vertcat(maxA_all{:});
lifetime_s = arrayfun(@(i) i.lifetime_s([i.catIdx]==1), lftData, 'UniformOutput', false);
lifetime_s = [lifetime_s{:}];

% generate cohorts
maxAcohorts = cell(1,nc);
for c = 1:nc
    idx = lb(c)<=lifetime_s & lifetime_s<=ub(c);
    maxAcohorts{c} = maxA(idx);
end

cmap = jet(nc);
dx = 40;
xi = 0:1:450;

figure;
hold on;

for c = 1:nc
    %ni = hist(maxAcohorts{c}, xi);
    %ni = ni/sum(ni)/dx;
    %plot(xi, ni, 'Color', cmap(c,:));
    [f, xi] = ksdensity(maxAcohorts{c}, xi);
    plot(xi, f, 'Color', cmap(c,:));
end

% compute convolution parameters
mu = zeros(1,nc-1);
sigma = zeros(1,nc-1);
f0 = ksdensity(maxAcohorts{1}, xi);
for c = 1:nc-1
    f = ksdensity(maxAcohorts{c+1}, xi);
    p = getGaussianConvPrms(xi, f0, f);
    mu(c) = p(1);
    sigma(c) = p(2);
end


figure;
hold on;
plot(lb(2:end), mu, 'r')
plot(lb(2:end), sigma)
legend('\mu', '\sigma', 'Location', 'NorthWest');
