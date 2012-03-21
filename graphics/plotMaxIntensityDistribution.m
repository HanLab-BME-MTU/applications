function plotMaxIntensityDistribution(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Mode', 'pdf', @(x) any(strcmpi(x, {'pdf', 'cdf'})));
ip.parse(varargin{:});
mode = ip.Results.Mode;

[~, figPath] = getCellDir(data(1));
figPath = [figPath 'Figures' filesep];
[~,~] = mkdir(figPath);


% cohorts: 4, 5, 6, 7, 8, 9, 10, 11-20, 21-40, 41-60 etc.
lb = [3:10 11 16 21 41 61 81 101 141];
ub = [3:10 15 20 40 60 80 100 140 200]; 

dxi = 5;
xi = 0:dxi:400;

%%
% res = getMaxIntensityDistributions(data, lb, ub);
% load('dataOXmaxInt_Ia.mat');
load('dataOXmaxInt_allTracks.mat');
% res = res([1:3 5:6 8:10]); % remove outliers

[a medIdx] = rescaleEDFs({res.maxA_all}, 'Display', false);


for i=1:numel(res)
    res(i).maxA = cellfun(@(x) x*a(i), res(i).maxA, 'UniformOutput', false);
    res(i).maxA4 = cellfun(@(x) x*a(i), res(i).maxA4, 'UniformOutput', false);
    %res(i).maxA_all = res(i).maxA_all*a(i);
end

%%
pres = struct([]);
for k = 1:numel(lb)
    tmp = vertcat(res.maxA);
    pres(k).maxA = horzcat(tmp{:,k});
    
    tmp = vertcat(res.maxA4);
    pres(k).maxA4 = horzcat(tmp{:,k});
    
    ni = hist(pres(k).maxA, xi);
    pres(k).ni = ni/sum(ni)/dxi;
    ni4 = hist(pres(k).maxA4, xi);
    pres(k).ni4 = ni4/sum(ni4)/dxi;
    
    tmp = vertcat(res.lft);
    pres(k).lft = horzcat(tmp{:,k});
    
    %mu = mean(pres(k).maxA);
    %mu3 = mean((maxA-mu).^3);
    %mu4 = mean((maxA-mu).^4);
    %sigma = std(pres(k).maxA);
    %pres(k).skew = mu3/sigma^3;
    %pres(k).kurt = mu4/sigma^4-3;
    
    if lb(k)==ub(k)
        pres(k).cohortLabel = num2str(lb(k));
    else
        pres(k).cohortLabel = [num2str(lb(k)) ' - ' num2str(ub(k))];
    end
end

% [k0, nVec, xVec, fVec, FVec, aVec] = fitGammaDistN({pres(:).maxA});
% aVec = [ones(1,7) aVec];


%%
% verify whether scaling of A_pstd is correct
% figure;
% hold on;
% 
% dx = 1;
% xi = 0:dx:max(res(medIdx).maxA_all_pstd);
% ni = hist(res(medIdx).maxA_all_pstd, xi);
% ni = ni/sum(ni)/dx;
% %[ni,xi] = ksdensity(samples{medIdx}, 'npoints', 1000);
% plot(xi, ni, 'r.-', 'LineWidth', 2);
% idx = setdiff(1:nd, medIdx);
% for i = 1:nd-1
%     ni = hist(res(idx(i)).maxA_all_pstd*a(idx(i)), xi);
%     ni = ni/sum(ni)/dx;
%     %[ni,xi] = ksdensity(samples{idx(i)}, 'npoints', 1000);
%     plot(xi, ni, 'k.-', 'LineWidth', 2);
% end

%%

% estimate Gaussian convolution that mediates transition from one cohort to the next

% for k = 1:numel(lb)-1
%     p = getGaussianConvPrms(xi, pres(k).ni, pres(k+1).ni);
%     mu(k) = p(1);
%     sigma(k) = p(2);
% end
% % meanLft = mean([lb; ub],1);
% dt = ub-lb+1;
% figure; plot(lb(1:end-1), sigma, '.-');
% figure; plot(lb(1:end-1), sigma./dt(1:end-1), '.-');
% % figure; plot(lb(1:end-1), mu, '.-');


%%
cb = [0 0.8 0];
cf3 = [1 1 1]*0.6;
ce3 = [1 1 1]*0.3;
tfont = {'FontName', 'Helvetica', 'FontSize', 14};
sfont = {'FontName', 'Helvetica', 'FontSize', 14};
lfont = {'FontName', 'Helvetica', 'FontSize', 18};

ny = 7; % # axes in y

pos = get(0, 'DefaultFigurePosition');
pos(3:4) = [920 800];
figure('Position', pos, 'Color', [1 1 1], 'PaperPositionMode', 'auto');
hbg = axes('Units', 'pixels', 'Position', [0 0 pos(3:4)]);
% axis(hbg, [0 0 pos(3:4)])
hold(hbg, 'on');
for k = 1:numel(lb)
    aw = 6*40;
    hi = axes('Units', 'pixels', 'Position', [80+floor((k-1)/ny)*(170+115) (ny-mod(k-1,ny)-1)*100+70 aw 80]);
    hold on;
    box off;
    set(hi, 'XGrid', 'on', 'GridLineStyle', ':');
    xa = 0:40:360;
    pct = prctile(pres(k).maxA, [5 50 95]);
    
    %[kGamma(k), nGamma(k), x, f, a, kappa(k)] = fitGammaDist(pres(k).maxA);
    if strcmpi(mode, 'pdf')
        bar(xi, pres(k).ni, 'BarWidth', 1, 'FaceColor', cf3, 'EdgeColor', ce3, 'LineWidth', 1);
        %plot(xVec, fVec{k-ny}, 'r', 'LineWidth', 1.5);
        
        if ub(k)<10
            [mu_g(k) sigma_g(k) xg g] = fitGaussianModeToHist(xi, pres(k).ni);
            plot(xg, g, 'g', 'LineWidth', 1.5);
            plot(norminv(0.95, mu_g(k), sigma_g(k))*[1 1], [0 0.03], 'g--', 'LineWidth', 1.5);
        end
        
        ya = 0:0.01:0.05;
        plot(repmat(pct, [2 1]), repmat([0 0.03]', [1 numel(pct)]), 'r--', 'LineWidth', 1.5);
        axis([xa(1) xa(end) ya(1) ya(end)]);
        
        if lb(k)==ub(k)
            text(xa(end), ya(end), ['t_L = ' num2str(lb(k))],...
                'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', sfont{:});
        else
            text(xa(end), ya(end), ['t_L \in [' num2str(lb(k)) '...' num2str(ub(k)) ']'],...
                'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', sfont{:});
        end
        
    else
        [f_ecdf, t_ecdf] = ecdf(pres(k).maxA);
        %plot(t_ecdf, f_ecdf*aVec(k)+1-aVec(k), 'k', 'LineWidth', 1.5);
        plot(t_ecdf, f_ecdf, 'k', 'LineWidth', 1.5);
        ya = 0:0.2:1;
        axis([xa(1) xa(end) 0 1.05]);
        plot(xVec, FVec{k}, 'r');
        
        if lb(k)==ub(k)
            text(xa(end), ya(1), ['t_L = ' num2str(lb(k))],...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', sfont{:});
        else
            text(xa(end), ya(1), ['t_L \in [' num2str(lb(k)) '...' num2str(ub(k)) ']'],...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', sfont{:});
        end
        
    end
    set(hi, 'YTick', ya , 'LineWidth', 2, tfont{:}, 'Layer', 'bottom', 'TickDir', 'out',...
        'TickLength', [0.02 0], 'XTick', xa, 'XTickLabel', []);
    
        
    
    if mod(k,ny)==0 || k == numel(lb)
        xlabel('Max. fluo. intensity (A.U.)', lfont{:});
        set(hi, 'TickLength', [0.02 0], 'XTickLabel', xa);
    end
    if k>ny
        set(hi, 'YTickLabel', []);
    end
    
    if k==(ny+1)/2
       ylabel('Frequency', lfont{:}); 
    end
    
    if k>1
        [pval hval] = ranksum(pres(k-1).maxA, pres(k).maxA);
        %[hval pval] = kstest2(pres(k-1).maxA, pres(k).maxA);
        if hval==0 % indicate that the distributions are the same
            
            %80+floor((k-1)/7)*(125+115) (7-mod(k-1,7)-1)*100+70 aw 80
            % x: after box:
            x0 = 80+floor((k-1)/7)*(125+115) + aw + 10;
            y0 = (7-mod(k-1,7)-1)*100+70 + 50;
            %y0 = 70;
            %plot(hbg, [x0 x0], [0 200], 'm');
            plot(hbg, [x0 x0], [y0 y0+80], 'Color', cb, 'LineWidth', 2); % height: 80, spacer: 20
            plot(hbg, [x0-3 x0], [y0 y0]+80-1, 'Color', cb, 'LineWidth', 2);
            plot(hbg, [x0-3 x0], [y0 y0]+1, 'Color', cb, 'LineWidth', 2);
            text(x0+3, y0+36, '*', lfont{:}, 'Parent', hbg, 'VerticalAlignment', 'middle')
            
        end
    end
    
end

% bring background axes to front
ch = get(gcf, 'Children');
set(gcf, 'Children', [hbg; setdiff(ch, hbg)]);


% plot(hbg, rand(1,10))
set(hbg, 'Color', 'none', 'XLim', [0 pos(3)], 'YLim', [0 pos(4)]);
axis(hbg, 'off');


if numel(data)==1
c    filename = ['_' getCellDir(data)];
else
    filename = '_pooled';
end
% print('-depsc2', '-loose', [figPath 'maxIntensityDist' filename '.eps']);
%%
% set threshold, plot lifetime distributions for the two resulting classes

T_int = norminv(0.95, mu_g(2), sigma_g(2));
% T = mu_g(2);
xh = 1:data(1).movieLength;

nr = numel(res);
nx = data(1).movieLength;
histIn = zeros(nr, nx);
histOut = zeros(nr, nx);

histAll = zeros(nr, nx);

pctIn = zeros(1,nr);
pctOut = zeros(1,nr);

lftIn = cell(1,nr);
lftOut = cell(1,nr);

for i = 1:nr
    idx = a(i)*res(i).maxA_all>=T_int;
    
    lftIn{i} = res(i).lft_all(idx);
    lftOut{i} = res(i).lft_all(~idx);
    lftIn{i}(lftIn{i}<3) = [];
    lftOut{i}(lftOut{i}<3) = [];
    
    pctIn(i) = sum(idx)/numel(idx);
    pctOut(i) = sum(~idx)/numel(idx);
    
    histIn(i,:) = hist(lftIn{i}, xh);
    histIn(i,:) = histIn(i,:)/sum(histIn(i,:));
    histOut(i,:) = hist(lftOut{i}, xh);
    histOut(i,:) = histOut(i,:)/sum(histOut(i,:));
end


lftIn = [lftIn{:}];
lftOut = [lftOut{:}];
nhOut_mean = mean(histOut,1);
nhIn_mean = mean(histIn,1);

nhIn_mean(1:2) = NaN; % only use lifetimes >= 3 frames
nhOut_mean(1:2) = NaN;


%%
[mu, mu_std, a_exp] = fitExpPDF(lftOut(lftOut>2));

figure;
hold on;
plot(xh, nhOut_mean, 'k.-', 'LineWidth', 2, 'MarkerSize', 20);
plot(xh, nhIn_mean, 'r.-', 'LineWidth', 2, 'MarkerSize', 20);

plot(xh, 1/mu*exp(-1/mu*xh) / a_exp, 'b');

[mu, mu_std, a_exp] = fitExpToHist(xh, nhOut_mean);
plot(xh, 1/mu * exp(-1/mu*xh) / a_exp, 'm');

axis([0 100 0 0.15]);
set(gca, 'LineWidth', 1.5, sfont{:});
xlabel('Lifetime (s)', lfont{:});
ylabel('Frequency', lfont{:});
legend(['% tracks: ' num2str(mean(pctOut)*100, '%.1f') ' ± ' num2str(std(pctOut)*100, '%.1f')],...
    ['% tracks: ' num2str(mean(pctIn)*100, '%.1f') ' ± ' num2str(std(pctIn)*100, '%.1f')], 'Location', 'NorthEast');

% print('-depsc2', '-loose', [figPath 'lftSep_T=' num2str(T, '%.1f') '.eps']);

%%
% Figure with split for a range of thresholds, shown in color gradient plots

T_lft = 4;

% tvec = 40:10:200; % intensity threshold
% tvec = 60;
tvec = [50 60 70 80];
% tvec = norminv(0.95, mu_g(2), sigma_g(2));


nt = numel(tvec);
meanHistIn = zeros(nt,nx);
meanHistOut = zeros(nt,nx);
meanHistAll = zeros(nt,nx);


for ti = 1:numel(tvec);
    T_int = tvec(ti);

    % pool all data sets
    lftIn = cell(1,nr);
    lftOut = cell(1,nr);
    for i = 1:nr
        idx = a(i)*res(i).maxA_all>=T_int;
        
        lftIn{i} = res(i).lft_all(idx);
        lftOut{i} = res(i).lft_all(~idx);
        lftIn{i}(lftIn{i}<T_lft) = [];
        lftOut{i}(lftOut{i}<T_lft) = [];
        
        %pctIn(i) = sum(idx)/numel(idx);
        %pctOut(i) = sum(~idx)/numel(idx);
        
        histAll(i,:) = hist([lftIn{i} lftOut{i}], xh);
        
        % normalization factor
        normf = sum(histAll(i,:));        
        histAll(i,:) = histAll(i,:)/normf;
        
        
        histIn(i,:) = hist(lftIn{i}, xh);
        histIn(i,:) = histIn(i,:)/normf;
        histOut(i,:) = hist(lftOut{i}, xh);
        histOut(i,:) = histOut(i,:)/normf;
    end
    
    meanHistIn(ti,:) = mean(histIn,1);
    meanHistOut(ti,:) = mean(histOut,1);
    meanHistAll(ti,:) = mean(histAll,1);
    meanHistIn(ti,1:T_lft-1) = NaN;
    meanHistOut(ti,1:T_lft-1) = NaN;
    meanHistAll(ti,1:T_lft-1) = NaN;
    
    [mu, ~, A_exp, xpdf] = fitExpToHist(xh, meanHistOut(ti,:));
    expfitOut(ti,:) = xpdf;%A_exp * 1/mu * exp(-1/mu*xh) ;
%     disp('');
%     [mu, mu_std, a_exp] = fitExpPDF([lftOut{:}]);
%     expFit(ti,:) = 1/mu*exp(-1/mu*xh) / a_exp;
end

figure;
hold on;
plot(xh, meanHistAll, 'k');
plot(xh, meanHistOut, 'r');

plot(xh, expfitOut, 'c');

plot(xh, meanHistIn, 'g');
plot(xh, meanHistAll, 'k');
% plot(xh, expFit, 'b');

% plot(xh, meanHistIn+meanHistOut, 'm--');
set(gca, 'XLim', [0 40]);
set(gca, 'YLim', [0 0.1]);
% axis([0 80 0 0.4]);
%%

% figure;
% hold on;
% plot3(tvec(1)*ones(size(xh)), xh, meanHistOut(1,:));
% 
% % figure; waterfall(meanHistOut)

%%
figure('Position', pos, 'Color', [1 1 1], 'PaperPositionMode', 'auto');

for k = 1:numel(lb)
    
    % for each lifetime cohort, get cumulative intensity of the first 4 frames
    M = arrayfun(@(i) vertcat(i.int4{k}{:}), res(1), 'UniformOutput', false);
    M = vertcat(M{:});   
    
    Tindex = arrayfun(@(i) vertcat(i.maxA{k}), res(1), 'UniformOutput', false);
    Tindex = horzcat(Tindex{:});
    
    aw = 3*40;
    hi = axes('Units', 'pixels', 'Position', [80+floor((k-1)/ny)*(170+115) (ny-mod(k-1,ny)-1)*100+70 aw 80]);
    hold on;
    box off;

    xa = 1:4;
    ya = 0:200:1000;
    
    %plot(M', 'k');
    
    plot(M(Tindex>=T_int,:)', 'g-');
    plot(M(Tindex<T_int,:)', 'r--');

    
    plot(cumsum(T_int*ones(1,4)), 'k', 'LineWidth', 2);
   
    %plot(repmat(pct, [2 1]), repmat([0 0.03]', [1 numel(pct)]), 'r--', 'LineWidth', 1.5);
    axis([0.5 4.5 ya(1) ya(end)]);
    %set(gca, 'XLim', [0.5 4.5]);
    
    if lb(k)==ub(k)
        text(5, ya(end), ['t_L = ' num2str(lb(k))],...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', sfont{:});
    else
        text(5, ya(end), ['t_L \in [' num2str(lb(k)) '...' num2str(ub(k)) ']'],...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', sfont{:});
    end
    
    
    set(hi, 'YTick', ya , 'LineWidth', 2, tfont{:}, 'Layer', 'bottom', 'TickDir', 'out',...
       'TickLength', [0.02 0], 'XTick', xa, 'XTickLabel', []);
    
        
    
    if mod(k,ny)==0 || k == numel(lb)
        xlabel('Time (s)', lfont{:});
        set(hi, 'TickLength', [0.02 0], 'XTickLabel', xa);
    end
    if k>ny
        set(hi, 'YTickLabel', []);
    end
    
    if k==(ny+1)/2
       ylabel('Cumulative intensity', lfont{:}); 
    end
    
   
    
end








%%

% meanLft = mean([lb; ub],1);
% 
% pos = get(0, 'DefaultFigurePosition');
% pos(3) = 400;
% pos(4) = 300;
% 
% figure('Units', 'pixels', 'Color', [1 1 1], 'PaperPositionMode', 'auto');
% % axes('Units', 'Pixels', 'Position', [80 80 300 200]);
% plot(meanLft, nVec, 'k.-', 'LineWidth', 2, 'MarkerSize', 16);
% set(gca, 'Box', 'off', sfont{:}, 'LineWidth', 2, 'TickDir', 'out');
% % set(gca, 'XTick', meanLft, 'XTickLabel', {pres.cohortLabel});
% xlabel('Lifetime (s)', tfont{:});
% ylabel('n', tfont{:});
% 
% print('-depsc2', '-loose', [figPath 'maxIntensityDist' filename '_gammaN.eps']);
% 

% pos(3) = 800;
% 
% figure('Units', 'pixels', 'Position', pos, 'Color', [1 1 1], 'PaperPositionMode', 'auto');
% axes('Units', 'Pixels', 'Position', [80 80 300 200]);
% plot(meanLft, 1./kGamma, 'k.-', 'LineWidth', 2, 'MarkerSize', 16);
% set(gca, 'Box', 'off', sfont{:}, 'LineWidth', 2, 'TickDir', 'out');
% % set(gca, 'XTick', meanLft, 'XTickLabel', {pres.cohortLabel});
% xlabel('Lifetime (s)', tfont{:});
% ylabel('k (s^{-1})', tfont{:});
% 
% axes('Units', 'Pixels', 'Position', [440 80 300 200]);
% plot(meanLft, nGamma, 'k.-', 'LineWidth', 2, 'MarkerSize', 16);
% set(gca, 'Box', 'off', sfont{:}, 'LineWidth', 2, 'TickDir', 'out');
% xlabel('Lifetime (s)', tfont{:});
% ylabel('n', tfont{:});
% 
% print('-depsc2', '-loose', [figPath 'maxIntensityDist' filename '_gParam.eps']);


