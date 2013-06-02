function plotMaxIntensityDistribution(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Mode', 'pdf', @(x) any(strcmpi(x, {'pdf', 'cdf'})));
ip.addParamValue('XTick', []);
ip.addParamValue('YTick', []);
ip.addParamValue('FirstNFrames', 5);
ip.addParamValue('CohortLB', [5  11 16 21 41 61 81]);
ip.addParamValue('CohortUB', [10 15 20 40 60 80 120]);
ip.addParamValue('DisplayMode', 'screen', @(x) any(strcmpi(x, {'screen', 'print'})));
ip.addParamValue('ShowSignificance', false, @islogical);
ip.addParamValue('ShowGaussians', false, @islogical);
ip.addParamValue('ShowFirstFrame', false, @islogical);
ip.addParamValue('LifetimeData', 'lifetimeData.mat');
ip.addParamValue('Cutoff_f', 5);
ip.parse(varargin{:});

[~, figPath] = getCellDir(data(1));
figPath = [figPath 'Figures' filesep];
[~,~] = mkdir(figPath);

mCh = strcmp(data(1).source, data(1).channels);

lb = ip.Results.CohortLB;
ub = ip.Results.CohortUB;
nc = numel(lb);

ny = nc;

lftData = getLifetimeData(data, 'ReturnValidOnly', true, 'Scale', true,...
    'Cutoff_f', 5, 'ExcludeVisitors', false, 'LifetimeData', ip.Results.LifetimeData);
A = arrayfun(@(i) i.A(:,:,mCh), lftData, 'unif', 0);
A = vertcat(A{:});
% maxA = vertcat(lftData.maxA);
maxA = nanmax(A,[],2);
f = ip.Results.FirstNFrames;

maxAFirstN = nanmax(A(:,1:f), [], 2);

lifetime_s = vertcat(lftData.lifetime_s);

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

%========================================
% Generate cohorts
%========================================
maxAcohort = cell(1,nc);
maxAcohortFirstN = cell(1,nc);
A0cohort = cell(1,nc);
ni = cell(1,nc);
niFirstN = cell(1,nc);
niFirst = cell(1,nc);
for k = 1:numel(lb)
    cidx = lb(k)<=lifetime_s & lifetime_s<ub(k);
    maxAcohort{k} = maxA(cidx);
    maxAcohortFirstN{k} = maxAFirstN(cidx);
    A0cohort{k} = A(cidx,1);
    
    ni0 = hist(maxAcohort{k}, xi);
    ni{k} = ni0/sum(ni0);
    ni0 = hist(maxAcohortFirstN{k}, xi);
    niFirstN{k} = ni0/sum(ni0);
    ni0 = hist(A0cohort{k}, xi);
    niFirst{k} = ni0/sum(ni0);
end

if isempty(ip.Results.YTick)
    di = 3;
    ymax = max(cellfun(@(i) max(i), niFirstN));
    mag = 10^floor(log10(ymax/di));
    dy = ceil(ymax/mag/di)*mag;
    ya = (0:di)*dy;
else
    ya = ip.Results.YTick;
end

cf0 = [1 1 1]*0.6;
ce0 = [1 1 1]*0.3;
fset = loadFigureSettings('print');

axPos = fset.axPos;
aw = axPos(3);
ah = 0.2*aw;
xo = 1.5;
yo = 1.5;
sh = ah/3;

% figure('Position', pos, 'Color', [1 1 1], 'PaperPositionMode', 'auto');
figure('Units', 'centimeters', 'Position', [5 5 8 ny*ah+(ny-1)*sh+2.5], 'Color', [1 1 1], 'PaperPositionMode', 'auto');
hbg = axes('Units', 'centimeters', 'Position', [0 0 2*xo+aw+2 2*yo+ny*ah+(ny-1)*sh]);
hold(hbg, 'on');
axis(hbg, [0 2*xo+aw+2 0 2*yo+ny*ah+(ny-1)*sh]);
for k = 1:nc
    
    iPos = [xo yo+(ny-k)*(ah+sh) aw ah];
    % plot grid below data
    hi = axes(fset.axOpts{:}, 'Position', iPos,...
        'XLim', [xa(1) xa(end)], 'XTick', xa, 'XTickLabel', [], 'YLim', [ya(1) ya(end)], 'YTickLabel', [],...
        'TickLength', [0 0], 'Color', 'none'); %#ok<LAXES>
    set(hi, 'XGrid', 'on', 'GridLineStyle', ':', 'LineWidth', 1);
    
    % data axes
    hi = axes(fset.axOpts{:}, 'Position', iPos,...
        'Color', 'none'); %#ok<LAXES>
    hold on;
    box off;
    
    if ip.Results.ShowFirstFrame
        bar(xi, niFirst{k}, 'BarWidth', 1, 'FaceColor', 0.9*[1 1 1], 'EdgeColor', 0.45*[1 1 1], 'LineWidth', 0.75);
    end
    bar(xi, niFirstN{k}, 'BarWidth', 1, 'FaceColor', cf0, 'EdgeColor', ce0, 'LineWidth', 0.75);
    bar(xi, ni{k}, 'BarWidth', 1, 'FaceColor', fset.cfB, 'EdgeColor', fset.ceB, 'LineWidth', 0.75);
    stairsXT(xi, niFirstN{k}, 'EdgeColor', ce0, 'LineWidth', 0.75);
    if ip.Results.ShowFirstFrame
        stairsXT(xi, niFirst{k}, 'EdgeColor', ce0, 'LineWidth', 0.75);
    end
    
    if ip.Results.ShowGaussians
        [mu_g,sigma_g,xg,g] = fitGaussianModeToHist(xi, niFirstN{k}/dxi);
        %[mu_g,sigma_g,xg,g] = fitGaussianModeToCDF(maxAcohortFirstN{k});
        plot(xg, g*dxi, 'Color', hsv2rgb([0 1 0.9]), 'LineWidth', 1);
        plot(norminv(0.99, mu_g, sigma_g)*[1 1], [0 0.15], '--', 'Color', hsv2rgb([0 1 0.9]), 'LineWidth', 1);
    end
    axis([xa(1) xa(end) ya(1) ya(end)]);
    
    % cohort label
    if lb(k)==ub(k)
        text(xa(end), ya(end), ['' num2str(lb(k)) ' s'], 'BackgroundColor', [1 1 1],...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', fset.sfont{:});
    else
        text(xa(end), ya(end), ['[' num2str(lb(k)) '...' num2str(ub(k)) '] s'],...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', fset.sfont{:},...
            'BackgroundColor', [1 1 1]);
    end
    
    if k==1
        text(xa(end), ya(end)*1.35, 'Lifetime cohort',...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', fset.sfont{:});
    end
    yal = [0 arrayfun(@(i) num2str(i, '%.2f'), ya(2:end), 'UniformOutput', false)];
    set(hi, 'YTick', ya, 'YTickLabel', yal, 'XTick', xa, 'XTickLabel', []);
    
    % plot vertical axis label (once)
    ncol = ceil(numel(lb)/ny);
    if mod(k,ny)==0
        set(hi, 'TickLength', [0.015 0], 'XTickLabel', xa);
        if floor(k/ny)==floor(ncol/2)+1
            hx = xlabel('Max. fluo. intensity (A.U.)', fset.lfont{:});
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
        hy = ylabel('Frequency', fset.lfont{:});
        ypos = get(hy, 'Position');
        if mod(ny,2)==0
            ypos(2) = 0;
        end
        set(hy, 'Position', ypos);
    end
    
    if k>1 && ip.Results.ShowSignificance % indicate which distributions are the same
        [hval, ~] = kstest2(maxAcohortFirstN{k-1}, maxAcohortFirstN{k});
        if hval==0
            % x: after box:
            x0 = xo + aw + 0.5;
            y0 = (ny-k)*(ah+sh) + yo + (ah+sh/2)/2;
            % line width in cm
            lw = 1/(72/2.54);
            plot(hbg, [x0 x0], y0+[0 ah+sh/2], 'Color', ce0, 'LineWidth', 1);
            plot(hbg, [x0-0.1 x0], [y0 y0]+ah+sh/2-lw/2, 'Color', ce0, 'LineWidth', 1);
            plot(hbg, [x0-0.1 x0], [y0 y0]+lw/2, 'Color', ce0, 'LineWidth', 1);
            text(x0+0.1, y0+(ah+sh/2)/2-0.1, '*', fset.lfont{:}, 'Color', ce0, 'Parent', hbg, 'VerticalAlignment', 'middle')
        end
    end
end

% bring background axes to front
% ch = get(gcf, 'Children');
% set(gcf, 'Children', [hbg; setdiff(ch, hbg)]);
axis(hbg, 'off');
