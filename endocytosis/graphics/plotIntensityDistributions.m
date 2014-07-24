function plotIntensityDistributions(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('XTick', []);
ip.addParamValue('CohortLB', [1 11 16 21 41 61]);
ip.addParamValue('CohortUB', [10 15 20 40 60 120]);
ip.addParamValue('ShowPct', false, @islogical);
ip.addParamValue('FigureName', '');
ip.addParamValue('LifetimeData', 'lifetimeData.mat');
ip.parse(varargin{:});
lb = ip.Results.CohortLB;
ub = ip.Results.CohortUB;

ivec = 0:3:120; % intensity range
tvec = 0:10; % time vector

lftData = getLifetimeData(data, 'LifetimeData', ip.Results.LifetimeData,...
    'Cutoff_f', 5, 'ExcludeVisitors', false, 'Scale', true);

A = vertcat(lftData.A);
lft = vertcat(lftData.lifetime_s);

fset = loadFigureSettings('print');
figure(fset.fOpts{:}, 'Name', ip.Results.FigureName);
if ip.Results.ShowPct
    set(gcf, 'Position', [5 5 8 7]);
    offset = 1.5;
else
    offset = 0;
end
iset = [fset.axOpts, 'XTick', 0:5:20, 'YTick', 0:40:120, 'XLim', [tvec(1)-0.5 tvec(end)+0.5], 'YLim', [ivec(1) ivec(end)], 'TickLength', fset.TickLength*6/1.8];
colormap(jet(256));
wx = 1.8;
wy = 1.6;
d0 = 0.3;
c = 1;

nc = numel(lb);
cv = jet(nc);
M = cell(1,nc);
for i = 2:-1:1
    for j = 1:3
        axes(iset{:}, 'Position', [1.5+(j-1)*(wx+d0) offset+1.5+(i-1)*(wy+d0) wx wy]); hold on;
        M{c} = A(lb(c)<=lft&lft<=ub(c),1:numel(tvec));
        T = repmat(tvec, [size(M{c},1),1]);
        mv = M{c}(:);
        tv = T(:);
        rmIdx = mv>ivec(end);
        mv(rmIdx) = [];
        tv(rmIdx) = [];
        hm = hist3([mv tv], {ivec, tvec});
        hm = hm./repmat(sum(hm,1), [numel(ivec) 1]);
        imagesc(tvec, ivec, hm);
        if ip.Results.ShowPct
            hold on;
            %stairsXT(tvec, prctile(M{c},95,1), 'bounds', 'open', 'EdgeColor', cv(c,:));
            stairsXT(tvec, prctile(M{c},50,1), 'bounds', 'open', 'EdgeColor', cv(c,:));
            %stairsXT(tvec, prctile(M{c},5,1), 'bounds', 'open', 'EdgeColor', cv(c,:));
        end
        text(tvec(end), 0.975*ivec(end), ['[' num2str(lb(c)) '...' num2str(ub(c)) '] s'],...
            'Color', 'w', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top',...
            fset.tfont{:}, 'FontWeight', 'bold')
        
        if c==1
            text(tvec(end), 1.2*ivec(end), 'Lifetime cohort', 'Color', 'k',...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', fset.tfont{:})
        end
        if i~=1
            set(gca, 'XTickLabel', []);
        end
        if j~=1
            set(gca, 'YTickLabel', []);
        end
        if i==1 && j==1
            yp = ylabel('Fluo. intensity (A.U.)', fset.lfont{:});
            ypos = get(yp, 'Position');
            ypos(2) = 140;
            set(yp, 'Position', ypos);
        end
        if i==1 && j==2
            xlabel('Time (s)', fset.lfont{:});
        end
        
        axes(iset{:}, 'Position', [1.5+(j-1)*(wx+d0) offset+1.5+(i-1)*(wy+d0) wx wy],...
            'XTick', [], 'YTick', [], 'Color', 'none');
        box on;
        
        c = c+1;
    end
end

if ip.Results.ShowPct
    axes(iset{:}, 'Position', [1.5+(1-1)*(wx+d0) 1.25 wx wy/3]); hold on;
    hold on;
    for c = 1:nc
        if ip.Results.ShowPct
            hold on;
            %stairsXT(tvec, prctile(M{c},95,1), 'bounds', 'open', 'EdgeColor', cv(c,:));
            stairsXT(tvec, prctile(M{c},50,1), 'bounds', 'open', 'EdgeColor', cv(c,:));
            %stairsXT(tvec, prctile(M{c},5,1), 'bounds', 'open', 'EdgeColor', cv(c,:));
        end
    end
    set(gca, 'YLim', [40 80], 'YTick', [40 80]);
end
