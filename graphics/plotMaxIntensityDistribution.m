function plotMaxIntensityDistribution(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Mode', 'pdf', @(x) any(strcmpi(x, {'pdf', 'cdf'})));
ip.addParamValue('XTicks', 0:40:360);
ip.addParamValue('FirstNFrames', []);
ip.parse(varargin{:});

mode = ip.Results.Mode;
xa = ip.Results.XTicks;

[~, figPath] = getCellDir(data(1));
figPath = [figPath 'Figures' filesep];
[~,~] = mkdir(figPath);

% cohorts: 3, 4, 5, 6, 7, 8, 9, 10, 11-20, 21-40, 41-60 etc.
lb = [3:10 11 16 21 41 61 81 101 141];
ub = [3:10 15 20 40 60 80 100 140 200]; 

da = xa(2)-xa(1);

% dxi = 5;
dxi = da/8;
% xi = 0:dxi:400;
xi = 0:dxi:xa(end)+da;

ya = 0:0.01:0.05;
% ya = 0:0.05:0.25;
% ya = 0:0.1:0.5;


lftData = getLifetimeData(data);
maxA_all = arrayfun(@(i) nanmax(i.intMat_Ia,[],2), lftData, 'UniformOutput', false);

% Rescale EDFs (correction for FP-fusion expression level)
a = rescaleEDFs(maxA_all, 'Display', false);

% apply scaling
for i = 1:numel(data)
    lftData(i).intMat_Ia = a(i) * lftData(i).intMat_Ia;
    maxA_all{i} = a(i) * maxA_all{i};
end

% Concatenate data
if isempty(ip.Results.FirstNFrames) % full distribution
    maxA = vertcat(maxA_all{:});
else
    f = ip.Results.FirstNFrames;
    maxA = arrayfun(@(i) nanmax(i.intMat(:,1:f),[],2), lftData, 'UniformOutput', false);
    maxA = vertcat(maxA{:});
end
lifetime_s = arrayfun(@(i) i.lifetime_s([i.catIdx]==1), lftData, 'UniformOutput', false);
lifetime_s = [lifetime_s{:}];


% [k0, nVec, xVec, fVec, FVec, aVec] = fitGammaDistN({maxIntDistCat(:).maxA});
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
% cb = [0 0.8 0];
cb = [0 0 0];
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


nc = numel(lb);
ni = cell(1,nc);
maxAcohort = cell(1,nc);
for k = 1:numel(lb)
    aw = 6*40;
    hi = axes('Units', 'pixels', 'Position', [80+floor((k-1)/ny)*(170+115) (ny-mod(k-1,ny)-1)*100+70 aw 80]);
    hold on;
    box off;
    set(hi, 'XGrid', 'on', 'GridLineStyle', ':');
    maxAcohort{k} = maxA(lb(k)<=lifetime_s & lifetime_s<=ub(k));
    pct = prctile(maxAcohort{k}, [5 50 95]);
    
    %[kGamma(k), nGamma(k), x, f, a, kappa(k)] = fitGammaDist(maxIntDistCat(k).maxA);
    if strcmpi(mode, 'pdf')
        ni0 = hist(maxAcohort{k}, xi);
        ni{k} = ni0/sum(ni0)/dxi;
        
        bar(xi, ni{k}, 'BarWidth', 1, 'FaceColor', cf3, 'EdgeColor', ce3, 'LineWidth', 1);
        %stairsXT(xi, ni{k}, 'FaceColor', cf3, 'EdgeColor', ce3, 'LineWidth', 1);
        
        if ub(k)<10
            [mu_g(k) sigma_g(k) xg g] = fitGaussianModeToHist(xi, ni{k});
            plot(xg, g, 'g', 'LineWidth', 1.5);
            plot(norminv(0.95, mu_g(k), sigma_g(k))*[1 1], [0 3/5*ya(end)], 'g--', 'LineWidth', 1.5);
        end
        
        plot(repmat(pct, [2 1]), repmat([0 3/5*ya(end)]', [1 numel(pct)]), 'r--', 'LineWidth', 1.5);
        axis([xa(1) xa(end) ya(1) ya(end)]);
        
        if lb(k)==ub(k)
            text(xa(end), ya(end), ['t_L = ' num2str(lb(k))],...
                'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', sfont{:});
        else
            text(xa(end), ya(end), ['t_L \in [' num2str(lb(k)) '...' num2str(ub(k)) ']'],...
                'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', sfont{:});
        end
        
%         if isempty(ip.Results.FirstNFrames) && k>1
%             [p, y] = getGaussianConvPrms(xi, ni{k-1}, ni{k});
%             mu(k-1) = p(1);
%             sigma(k-1) = p(2);
%             stairsXT(xi, ni{k-1}, 'EdgeColor', 0.2*[1 1 1]);
%             plot(xi, y, 'b');
%         end
        
    else
        [f_ecdf, t_ecdf] = ecdf(maxIntDistCat(k).maxA);
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
        %[pval hval] = ranksum(maxIntDistCat(k-1).maxA, maxIntDistCat(k).maxA);
        [hval pval] = kstest2(maxAcohort{k-1}, maxAcohort{k});
        %pval
        if hval==0 % indicate that the distributions are the same
            
            %pos = [80+floor((k-1)/ny)*(170+115) (ny-mod(k-1,ny)-1)*100+70 aw 80]
            % x: after box:
            x0 = 80+floor((k-1)/ny)*(170+115) + aw + 10;
            y0 = (ny-mod(k-1,ny)-1)*100+70 + 50;

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
    filename = ['_' getCellDir(data)];
else
    filename = '_pooled';
end
% print('-depsc2', '-loose', [figPath 'maxIntensityDist' filename '.eps']);

%%
% plot result of Gaussian convolutions
% estimate Gaussian convolution that mediates transition from one cohort to the next

% meanLft = mean([lb; ub],1);
% cohortMean = (lb+ub)/2;
% cohortLength = ub-lb+1;
% figure; plot(lb(1:end-1), sigma, '.-'); xlabel('Previous cohort lower bound'); ylabel('sigma');
% figure; plot(lb(1:end-1), sigma./cohortLength(1:end-1), '.-'); ylabel('sigma/\Delta');
% figure; plot(lb(1:end-1), mu./cohortLength(1:end-1), '.-'); ylabel('mu');



%%
return
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
% Figure with split for a range of thresholds

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
%%

figure;
hold on;
% plot(xh, meanHistAll, 'k', 'LineWidth', 2);
tt = plot(xh, meanHistOut, 'r', 'LineWidth', 2);
hp(2) = tt(1);

tt = plot(xh, expfitOut, 'b', 'LineWidth', 1.5);
hp(3) = tt(1);

tt = plot(xh, meanHistIn, 'g', 'LineWidth', 2);
hp(1) = tt(1);

tt = plot(xh, meanHistAll, 'k--', 'LineWidth', 2);
hp(4) = tt(1);

% plot(xh, meanHistIn+meanHistOut, 'm--');
set(gca, 'XLim', [0 40]);
set(gca, 'YLim', [0 0.1]);
% axis([0 80 0 0.4]);
set(gca, 'LineWidth', 2, 'TickDir', 'out', 'Layer', 'top', 'FontSize', 16);
hl = legend(hp, 'Above threshold', 'Below threshold', 'Exp fit', 'All tracks');
set(hl, 'Box', 'off');
xlabel('Lifetime (s)', lfont{:})
ylabel('Frequency', lfont{:})

%%

% fit exponential to determine threshold






return

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
% % set(gca, 'XTick', meanLft, 'XTickLabel', {maxIntDistCat.cohortLabel});
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
% % set(gca, 'XTick', meanLft, 'XTickLabel', {maxIntDistCat.cohortLabel});
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


