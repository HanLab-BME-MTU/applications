function plotMaxIntensityDistribution(data)

[~, figPath] = getCellDir(data(1));
figPath = [figPath 'Figures' filesep];
[~,~] = mkdir(figPath);
% bins: 5-10


% cohorts: 4, 5, 6, 7, 8, 9, 10, 11-20, 21-40, 41-60

lb = [4:10 11 21 41 61 81 101 141];
ub = [4:10 20 40 60 80 100 140 200]; 

nc = numel(lb);

dxi = 5;
xi = 0:dxi:400;


nd = numel(data);
res = struct([]);
for i = 1:nd
    load([data(i).source 'Tracking' filesep 'trackAnalysis.mat']);
    
    
    for k = 1:nc
        lft = [tracks.lifetime_s];
        itracks = tracks(lb(k)<=lft & lft<=ub(k));
        % max intensities
        maxA = arrayfun(@(t) max(t.A(1,:)), itracks);
        ni = hist(maxA, xi);
        ni = ni/sum(ni)/dxi;
        
        res(i).maxA{k} = maxA;
        res(i).ni{k} = ni;
        
    end
    
end

%%

% pool data sets
pres = struct([]);
for k = 1:numel(lb)
    tmp = vertcat(res.maxA);
    maxA = horzcat(tmp{:,k});
    pres(k).maxA = maxA;
    ni = hist(pres(k).maxA, xi);
    pres(k).ni = ni/sum(ni)/dxi;
    mu = mean(pres(k).maxA);
    mu3 = mean((maxA-mu).^3);
    mu4 = mean((maxA-mu).^4);
    sigma = std(pres(k).maxA);
    pres(k).skew = mu3/sigma^3;
    pres(k).kurt = mu4/sigma^4-3;
    
    
    if lb(k)==ub(k)
        pres(k).cohortLabel = num2str(lb(k));
    else
        pres(k).cohortLabel = [num2str(lb(k)) ' - ' num2str(ub(k))];
    end
    
    
end

[k0, nVec, xVec, fVec, aVec] = fitGammaDistN({pres.maxA});

%%
cb = [0 0.8 0];
cf3 = [1 1 1]*0.6;
ce3 = [1 1 1]*0.3;
tfont = {'FontName', 'Helvetica', 'FontSize', 14};
sfont = {'FontName', 'Helvetica', 'FontSize', 14};
lfont = {'FontName', 'Helvetica', 'FontSize', 18};

pos = get(0, 'DefaultFigurePosition');
pos(3:4) = [560 800];
figure('Position', pos, 'Color', [1 1 1], 'PaperPositionMode', 'auto');
hbg = axes('Units', 'pixels', 'Position', [0 0 pos(3:4)]);
% axis(hbg, [0 0 pos(3:4)])
hold(hbg, 'on');
for k = 1:numel(lb)
    aw = 5*40;
    hi = axes('Units', 'pixels', 'Position', [80+floor((k-1)/7)*(125+115) (7-mod(k-1,7)-1)*100+70 aw 80]);
    
    [kGamma(k), nGamma(k), x, f, a, kappa(k)] = fitGammaDist(pres(k).maxA);
    bar(xi, pres(k).ni*aVec(k), 'BarWidth', 1, 'FaceColor', cf3, 'EdgeColor', ce3, 'LineWidth', 1);
    hold on;
    plot(xVec, fVec{k}, 'r', 'LineWidth', 1.5);
    box off;
    ya = 0:0.01:0.05;
    xa = 0:40:320;
    set(hi, 'YTick', ya , 'LineWidth', 2, tfont{:}, 'Layer', 'top', 'TickDir', 'out',...
        'TickLength', [0.02 0], 'XTick', xa, 'XTickLabel', []);
    
    p95 = prctile(res(1).maxA{k}, 95);
    plot([p95 p95], [0 0.03], 'r--', 'LineWidth', 1.5);
    axis([xa(1) xa(end) ya(1) ya(end)]);
    
    if lb(k)==ub(k)
        text(xa(end), ya(end), ['t_L = ' num2str(lb(k))],...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', sfont{:});
    else
        text(xa(end), ya(end), ['t_L \in [' num2str(lb(k)) '...' num2str(ub(k)) ']'],...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', sfont{:});
    end
    if mod(k,7)==0
        xlabel('Fluo. intensity (A.U.)', lfont{:});
        set(hi, 'TickLength', [0.02 0], 'XTickLabel', xa);
    end
    if k>7
        set(hi, 'YTickLabel', []);
    end
    
    if k==4
       ylabel('Frequency', lfont{:}); 
    end
    
    if k>1
        [hval pval] = kstest2(pres(k-1).maxA, pres(k).maxA);
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


if nd==1
    filename = ['_' getCellDir(data)];
else
    filename = '_pooled';
end
print('-depsc2', '-loose', [figPath 'maxIntensityDist' filename '.eps']);
%%

meanLft = mean([lb; ub],1);

pos = get(0, 'DefaultFigurePosition');
pos(3) = 400;
pos(4) = 300;

figure('Units', 'pixels', 'Color', [1 1 1], 'PaperPositionMode', 'auto');
% axes('Units', 'Pixels', 'Position', [80 80 300 200]);
plot(meanLft, nVec, 'k.-', 'LineWidth', 2, 'MarkerSize', 16);
set(gca, 'Box', 'off', sfont{:}, 'LineWidth', 2, 'TickDir', 'out');
% set(gca, 'XTick', meanLft, 'XTickLabel', {pres.cohortLabel});
xlabel('Lifetime (s)', tfont{:});
ylabel('n', tfont{:});

print('-depsc2', '-loose', [figPath 'maxIntensityDist' filename '_gammaN.eps']);


pos(3) = 800;

figure('Units', 'pixels', 'Position', pos, 'Color', [1 1 1], 'PaperPositionMode', 'auto');
axes('Units', 'Pixels', 'Position', [80 80 300 200]);
plot(meanLft, 1./kGamma, 'k.-', 'LineWidth', 2, 'MarkerSize', 16);
set(gca, 'Box', 'off', sfont{:}, 'LineWidth', 2, 'TickDir', 'out');
% set(gca, 'XTick', meanLft, 'XTickLabel', {pres.cohortLabel});
xlabel('Lifetime (s)', tfont{:});
ylabel('k (s^{-1})', tfont{:});

axes('Units', 'Pixels', 'Position', [440 80 300 200]);
plot(meanLft, nGamma, 'k.-', 'LineWidth', 2, 'MarkerSize', 16);
set(gca, 'Box', 'off', sfont{:}, 'LineWidth', 2, 'TickDir', 'out');
xlabel('Lifetime (s)', tfont{:});
ylabel('n', tfont{:});

print('-depsc2', '-loose', [figPath 'maxIntensityDist' filename '_gParam.eps']);


