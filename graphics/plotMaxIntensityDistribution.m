function plotMaxIntensityDistribution(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Mode', 'pdf', @(x) any(strcmpi(x, {'pdf', 'cdf'})));
ip.addParamValue('XTick', []);
ip.addParamValue('FirstNFrames', 5);
ip.addParamValue('CohortLB', [1  11 21 31 41 61 81 101]);
ip.addParamValue('CohortUB', [10 20 30 40 60 80 100 120]);
ip.parse(varargin{:});

mode = ip.Results.Mode;

[~, figPath] = getCellDir(data(1));
figPath = [figPath 'Figures' filesep];
[~,~] = mkdir(figPath);

% cohorts: 3, 4, 5, 6, 7, 8, 9, 10, 11-20, 21-40, 41-60 etc.
% lb = [3:10 11 16 21 41 61 81 101 141];
% ub = [3:10 15 20 40 60 80 100 140 200]; 

lb = ip.Results.CohortLB;
ub = ip.Results.CohortUB;
nc = numel(lb);

ny = min(4, ceil(nc/2)); % # axes in y

lftData = getLifetimeData(data);
maxA_all = arrayfun(@(i) nanmax(i.intMat_Ia,[],2), lftData, 'UniformOutput', false);

% Rescale EDFs (correction for FP-fusion expression level)
a = rescaleEDFs(maxA_all, 'Display', false);

% apply scaling
for i = 1:numel(data)
    lftData(i).intMat_Ia = a(i) * lftData(i).intMat_Ia;
    maxA_all{i} = a(i) * maxA_all{i};
end

% Concatenate maximum intensity and lifetime data
maxA = vertcat(maxA_all{:});

f = ip.Results.FirstNFrames;
maxAFirstN = arrayfun(@(i) nanmax(i.intMat_Ia(:,1:f),[],2), lftData, 'UniformOutput', false);
maxAFirstN = vertcat(maxAFirstN{:});

lifetime_s = arrayfun(@(i) i.lifetime_s([i.catIdx]==1), lftData, 'UniformOutput', false);
lifetime_s = [lifetime_s{:}];

% x-axis
xa = ip.Results.XTick;
if isempty(xa)
    db = prctile(maxA, [0.1 99.9]);
    mag = 10^floor(log10(db(2)));
    da = ceil(db(2)/mag)*mag/10;
    xa = (0:ceil(db(2)/da))*da;
else
    da = xa(2)-xa(1);
end
dxi = da/8;
xi = 0:dxi:xa(end)+da;

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
%========================================
% Generate cohorts
%========================================
maxAcohort = cell(1,nc);
maxAcohortFirstN = cell(1,nc);
ni = cell(1,nc);
niFirstN = cell(1,nc);
for k = 1:numel(lb)
    maxAcohort{k} = maxA(lb(k)<=lifetime_s & lifetime_s<=ub(k));
    maxAcohortFirstN{k} = maxAFirstN(lb(k)<=lifetime_s & lifetime_s<=ub(k));

    ni0 = hist(maxAcohort{k}, xi);
    ni{k} = ni0/sum(ni0);
    ni0 = hist(maxAcohortFirstN{k}, xi);
    niFirstN{k} = ni0/sum(ni0);
end

ymax = max(cellfun(@(i) max(i), niFirstN));
mag = 10^floor(log10(ymax/5));
dy = ceil(ymax/mag/5)*mag;
ya = (0:5)*dy;

%%
cb = [0 0 0];
cf0 = [1 1 1]*0.6;
ce0 = [1 1 1]*0.3;
fset = loadFigureSettings();


pos = get(0, 'DefaultFigurePosition');
pos(3:4) = [920 800];
figure('Position', pos, 'Color', [1 1 1], 'PaperPositionMode', 'auto');
hbg = axes('Units', 'pixels', 'Position', [0 0 pos(3:4)]);
hold(hbg, 'on');


for k = 1:numel(lb)
    aw = 6*40;
    
    hi = axes('Units', 'pixels', 'Position', [80+floor((k-1)/ny)*(170+115) (ny-mod(k-1,ny)-1)*100+70 aw 80],...
        'XLim', [xa(1) xa(end)], 'XTick', xa, 'XTickLabel', [], 'YLim', [ya(1) ya(end)], 'YTickLabel', [],...
        'TickLength', [0 0], 'Color', 'none');
    set(hi, 'XGrid', 'on', 'GridLineStyle', ':');
    
    hi = axes('Units', 'pixels', 'Position', [80+floor((k-1)/ny)*(170+115) (ny-mod(k-1,ny)-1)*100+70 aw 80],...
        'Color', 'none');
    hold on;
    box off;
    
    pct = prctile(maxAcohort{k}, [5 50 95]);
    
    if strcmpi(mode, 'pdf')
        
        bar(xi, niFirstN{k}, 'BarWidth', 1, 'FaceColor', cf0, 'EdgeColor', ce0, 'LineWidth', 1);
        bar(xi, ni{k}, 'BarWidth', 1, 'FaceColor', fset.cfB, 'EdgeColor', fset.ceB, 'LineWidth', 1);
        stairsXT(xi, niFirstN{k}, 'EdgeColor', ce0, 'LineWidth', 1);
        
        if ub(k)<10
            [mu_g(k) sigma_g(k) xg g] = fitGaussianModeToHist(xi, ni{k});
            plot(xg, g, 'g', 'LineWidth', 1.5);
            plot(norminv(0.95, mu_g(k), sigma_g(k))*[1 1], [0 3/5*ya(end)], 'g--', 'LineWidth', 1.5);
        end
        
        %plot(repmat(pct, [2 1]), repmat([0 3/5*ya(end)]', [1 numel(pct)]), 'r--', 'LineWidth', 1.5);
        axis([xa(1) xa(end) ya(1) ya(end)]);
        
        if lb(k)==ub(k)
            text(xa(end), ya(end), ['t_L = ' num2str(lb(k)) ' s'], 'BackgroundColor', [1 1 1],...
                'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', fset.ifont{:});
        else
            text(xa(end), ya(end), ['t_L \in [' num2str(lb(k)) '...' num2str(ub(k)) '] s'],...
                'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', fset.ifont{:},...
                'BackgroundColor', [1 1 1]);
        end
        
    else % CDF plot
        [f_ecdf, t_ecdf] = ecdf(maxIntDistCat(k).maxA);
        %plot(t_ecdf, f_ecdf*aVec(k)+1-aVec(k), 'k', 'LineWidth', 1.5);
        plot(t_ecdf, f_ecdf, 'k', 'LineWidth', 1.5);
        ya = 0:0.2:1;
        axis([xa(1) xa(end) 0 1.05]);
        plot(xVec, FVec{k}, 'r');
        
        if lb(k)==ub(k)
            text(xa(end), ya(1), ['t_L = ' num2str(lb(k))],...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', fset.ifont{:});
        else
            text(xa(end), ya(1), ['t_L \in [' num2str(lb(k)) '...' num2str(ub(k)) ']'],...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', fset.ifont{:});
        end
        
    end
    set(hi, 'YTick', ya, 'LineWidth', 1.5, fset.ifont{:}, 'Layer', 'top', 'TickDir', 'out',...
        'TickLength', [0.02 0], 'XTick', xa, 'XTickLabel', []);
    
    
    ncol = ceil(numel(lb)/ny);
    if mod(k,ny)==0
        set(hi, 'TickLength', [0.02 0], 'XTickLabel', xa);
        if floor(k/ny)==floor(ncol/2)+1
            hx = xlabel('Max. fluo. intensity (A.U.)', fset.tfont{:});
            if mod(ncol,2)==0
                xpos = get(hx, 'Position');
                xpos(1) = xa(1);
                xpos(2) = 1.2*xpos(2);
                set(hx, 'Position', xpos);
            end
        end
    end
    if k>ny
        set(hi, 'YTickLabel', []);
    end
    
    if k==ceil(ny/2)
        hy = ylabel('Frequency', fset.tfont{:});
        if mod(ny,2)==0
            ypos = get(hy, 'Position');
            ypos(2) = 0;
            ypos(1) = 1.2*ypos(1);
            set(hy, 'Position', ypos);
        end
    end
    
%     if k>1
%         %[pval hval] = ranksum(maxIntDistCat(k-1).maxA, maxIntDistCat(k).maxA);
%         [hval pval] = kstest2(maxAcohortFirstN{k-1}, maxAcohortFirstN{k});
%         %pval
%         if hval==0 % indicate that the distributions are the same
%             
%             %pos = [80+floor((k-1)/ny)*(170+115) (ny-mod(k-1,ny)-1)*100+70 aw 80]
%             % x: after box:
%             x0 = 80+floor((k-1)/ny)*(170+115) + aw + 10;
%             y0 = (ny-mod(k-1,ny)-1)*100+70 + 50;
% 
%             plot(hbg, [x0 x0], [y0 y0+80], 'Color', cb, 'LineWidth', 2); % height: 80, spacer: 20
%             plot(hbg, [x0-3 x0], [y0 y0]+80-1, 'Color', cb, 'LineWidth', 2);
%             plot(hbg, [x0-3 x0], [y0 y0]+1, 'Color', cb, 'LineWidth', 2);
%             text(x0+3, y0+36, '*', lfont{:}, 'Parent', hbg, 'VerticalAlignment', 'middle')            
%         end
%     end
end

% bring background axes to front
% ch = get(gcf, 'Children');
% set(gcf, 'Children', [hbg; setdiff(ch, hbg)]);


% plot(hbg, rand(1,10))
set(hbg, 'Color', 'none', 'XLim', [0 pos(3)], 'YLim', [0 pos(4)]);
axis(hbg, 'off');


if numel(data)==1
    filename = ['_' getCellDir(data)];
else
    filename = '_pooled';
end
% print('-depsc2', '-loose', [figPath 'maxIntensityDist' filename '.eps']);


return

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

