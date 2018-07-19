function lftData = simLifetimeData(k, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('k');
ip.addOptional('ns', 2e4);
ip.addParamValue('N', 10);
ip.addParamValue('cutoff', 4);
ip.addParamValue('Display', false, @islogical);
ip.parse(k, varargin{:});

% # states
ns = numel(k)+1;
np = ns/2;

% time vector
dt = 0.01;
t = 0:dt:400;

% Intializations & bounds
S0 = [1 zeros(1,ns-1)];


% ground truth
sol = ode45(@(t,y) getStateMatrix(np, k)*y, [0 t(end)], S0);
Y = deval(sol, t);

popPDF = zeros(np,numel(t));
for k = 1:np
    p = Y(2*k-1,:);
    popPDF(k,:) = p/sum(p)/dt * Y(2*k,end);
end
pdf = sum(popPDF, 1);

n = sum(pdf)*dt;
pdf = pdf/n;
popPDF = popPDF/n;
cdf = sum(Y(2:2:end,:),1);
popCDF = Y(2:2:end,:);

[ucdf, uidx] = unique(cdf);



dti = 1;
ti = ip.Results.cutoff:dti:t(end);

lftData.t = ti;


for n = 1:ip.Results.N

    samples = interp1(ucdf, t(uidx), rand(1, ip.Results.ns));
    
    %[f_ecdf, t_ecdf] = ecdf(samples);
    
    % histogram & cdf
    lftData.nSamples(n) = sum(samples>ip.Results.cutoff-dti/2);
    ni = hist(samples(samples>ip.Results.cutoff-dti/2), ti);
    ni = ni/sum(ni)/dti;
    cdfi = cumsum(ni)*dti;
    a = interp1(t, cdf, ti(1)-dti/2);
    
    % least-squares fit
    % sum((1-cdfi).*(interp1(t, CDF, ti)-cdfi)) / sum ((1-cdfi).^2)

    lftData.lftHist{n} = ni;
    
end


M = vertcat(lftData.lftHist{:});
meanHist = mean(M,1);
% histSEM = std(M,[],1) / sqrt(length(data));
% res.t = t_hist;
lftData.meanHist = meanHist;

lftData.meanECDF = cumsum(meanHist)*dti;

% [ti, ni, t, pdf, popPDF, cdf, popCDF]





if ip.Results.Display
    fset = loadFigureSettings();
    
    %------------------------------------
    % Display result of best fit
    %------------------------------------
    switch np
        case 1
            colorOrder = fset.ceG;
        case 2
            colorOrder = [fset.ceR; fset.ceG];
        case 3
            colorOrder = [fset.ceR; fset.ceB2; fset.ceG];
        case 4
            colorOrder = [fset.ceR; fset.ceB2; fset.ceB; fset.ceG];
    end
    colorOrderFill = rgb2hsv(colorOrder);
    colorOrderFill(:,2) = colorOrderFill(:,2)*0.3;
    colorOrderFill = hsv2rgb(colorOrderFill);
    
    figure('Position', [240 378 1150 360], 'PaperPositionMode', 'auto', 'Color', 'w', 'InvertHardcopy', 'off');
    %---------------------------------
    % Lifetime histogram
    %---------------------------------
    axes('Units', 'Pixels', 'Position', [95 65 450 270]);
    set(gca, 'ColorOrder', colorOrder);
    hold on;
    hp(1) = plot(ti, meanHist*(1-a), '.', 'MarkerSize', 20, 'Color', [0 0 0]);
    hi = plot(t, popPDF, 'LineWidth', 2);
    hp(2) = hi(1);
    hp(3) = plot(t, pdf, '--', 'Color', fset.ceB, 'LineWidth', 4);
    
    axis([0 100 0 0.05]);
    set(gca, 'LineWidth', 2, 'Layer', 'top', fset.sfont{:});
    xlabel('Lifetime (s)', fset.lfont{:});
    ylabel('Frequency', fset.lfont{:});
    
    %---------------------------------
    % Inset with amplitudes
    %---------------------------------
    
    ha = axes('Units', 'Pixels', 'Position', [95+450-32*np-10 260 32*np 65]);
    xlabels = arrayfun(@(i) ['P' num2str(i)], 1:np, 'UniformOutput', false);
    
    barplot2(popCDF(:,end)', 'AdjustFigure', false, 'XLabels', xlabels,...
        'FaceColor', colorOrderFill, 'EdgeColor', colorOrder,...
        'BarWidth', 0.6, 'GroupDistance', 0.5, 'Angle', 45);
    set(ha, fset.tfont{:}, 'YTick', 0:0.2:0.8, 'YLim', [0 0.8]);
    ylabel('Contrib.', fset.sfont{:})
  
    
    axes('Units', 'Pixels', 'Position', [95+95+450 65 450 270]);
    
    set(gca, 'ColorOrder', colorOrder);
    hold on;
    plot(ti+dti/2, cdfi*(1-a)+a, 'k.', 'MarkerSize', 20)
    plot(t, popCDF, 'LineWidth', 2);
    plot(t, cdf, '--', 'Color', fset.ceB, 'LineWidth', 4);
    
    axis([0 100 0 1]);
    set(gca, 'LineWidth', 2, 'Layer', 'top', fset.sfont{:});
    xlabel('Lifetime (s)', fset.lfont{:});
    ylabel('Frequency', fset.lfont{:});
    
%     plot(ti+dti/2, cdfi*(1-a)+a, 'k.');
%     plot(t, cdf, 'r');
%     plot(t_ecdf, f_ecdf, 'b--');
%     axis([0 100 0 1]);
%     legend('Cumulative hist.', 'CDF (model)', 'ECDF (samples)', 'Location', 'SouthEast')
end


% less precise
% a = 1 - interp1(t, CDF, ti(1)-dt/2);
% plot(ti, ni*a, 'k.');
% sum((x-a*ni).^2)

% a = sum(interp1(t, PDF, ti))*dti;
% hp(2) = plot(ti, ni*a, 'k.');
% sum((x-a*ni).^2)


% hp(1) = plot(t, PDF, 'r');

% set(gca, 'LineWidth', 1.5, 'TickDir', 'Out', 'Layer', 'top');

% legend(hp, 'Ground truth', 'Histogram of samples');
% axis([0 200 0 0.03]);

