function h = plotLifetimes(lftRes, colorA, mode)

if nargin<3
    mode = '';
end

fset = loadFigureSettings(mode);

if isstruct(lftRes)
    
    figure(fset.fOpts{:}, 'Name', 'Lifetime dist. (intensity threshold)');
    axes(fset.axOpts{:});
    hold on;
    hp(4) = plot(lftRes.t, mean(vertcat(lftRes.lftHist_Ia), 1), 'Color', 0.6*[1 1 1], 'LineWidth', 2);
    hp(3) = plot(lftRes.t, mean(lftRes.pctVisit)*lftRes.meanLftHist_V, '-', 'Color', fset.cfB, 'LineWidth', 2);
    hp(2) = plot(lftRes.t, mean(lftRes.pctBelow)*lftRes.meanLftHist_B, '-', 'Color', hsv2rgb([1/3 0.3 0.9]), 'LineWidth', 2);
    hp(1) = plot(lftRes.t, mean(lftRes.pctAbove)*lftRes.meanLftHist_A, '-', 'Color', hsv2rgb([1/3 1 0.9]), 'LineWidth', 2);
    
    ya = 0:0.01:0.05;
    axis([0 min(120, lftRes.t(end)) 0 ya(end)]);
    set(gca, 'XTick', 0:20:200, 'YTick', ya, 'YTickLabel', ['0' arrayfun(@(x) num2str(x, '%.2f'), ya(2:end), 'UniformOutput', false)]);
    xlabel('Lifetime (s)', fset.lfont{:});
    ylabel('Frequency', fset.lfont{:});
    
    hl = legend(hp, ['Above threshold: ' num2str(mean(lftRes.pctAbove)*100, '%.1f') ' ± ' num2str(std(lftRes.pctAbove)*100, '%.1f') ' %'],...
        ['Below threshold: ' num2str(mean(lftRes.pctBelow)*100, '%.1f') ' ± ' num2str(std(lftRes.pctBelow)*100, '%.1f') ' %'],...
        ['Visitors: ' num2str(mean(lftRes.pctVisit)*100, '%.1f') ' ± ' num2str(std(lftRes.pctVisit)*100, '%.1f') ' %'],...
        'Raw distribution');
    set(hl, 'Box', 'off', fset.tfont{:});
    
    
%     h = figure;
%     hold on;
%     opts = {'.-', 'LineWidth', 2, 'MarkerSize', 16};
%     hp(2) = plot(lftRes.t, lftRes.meanLftHist_B, opts{:}, 'Color', 0.5*[1 1 1]);
%     hp(1) = plot(lftRes.t, lftRes.meanLftHist_A, opts{:}, 'Color', 'k');
% 
%     legendText = {['Above threshold (' num2str(mean(lftRes.pctAbove)*100,'%.1f') ' ± ' num2str(std(lftRes.pctAbove)*100,'%.1f') ' %)'],...
%         ['Below threshold (' num2str(mean(1-lftRes.pctAbove)*100,'%.1f') ' ± ' num2str(std(lftRes.pctAbove)*100,'%.1f') ' %)']};
%     
%     if isfield(lftRes, 'lftHist_Apos')
%         hp(6) = plot(lftRes.t, mean(lftRes.lftHist_Bneg,1), opts{:}, 'Color', hsv2rgb([0   0.4 0.9]));
%         hp(5) = plot(lftRes.t, mean(lftRes.lftHist_Bpos,1), opts{:}, 'Color', hsv2rgb([1/3 0.4 0.9]));
%         hp(4) = plot(lftRes.t, mean(lftRes.lftHist_Aneg,1), opts{:}, 'Color', hsv2rgb([0   1 0.9]));
%         hp(3) = plot(lftRes.t, mean(lftRes.lftHist_Apos,1), opts{:}, 'Color', hsv2rgb([1/3 1 0.9]));
%         legendText = [legendText ['Above, sign. (' num2str(mean(lftRes.pctAboveSignificant)*100,'%.1f') ' ± ' num2str(std(lftRes.pctAboveSignificant)*100,'%.1f') ' %)'],...
%             ['Above, not sign. (' num2str(mean(lftRes.pctAbove-lftRes.pctAboveSignificant)*100,'%.1f') ' ± ' num2str(std(lftRes.pctAboveSignificant)*100,'%.1f') ' %)'],...
%             ['Below, sign. (' num2str(mean(lftRes.pctBelowSignificant)*100,'%.1f') ' ± ' num2str(std(lftRes.pctBelowSignificant)*100,'%.1f') ' %)'],...
%             ['Below, not sign. (' num2str(mean(1-lftRes.pctAbove-lftRes.pctBelowSignificant)*100,'%.1f') ' ± ' num2str(std(lftRes.pctBelowSignificant)*100,'%.1f') ' %)']];
%     end
%     
%     axis([0 min(120, lftRes.t(end)) 0 0.05]);
%     set(gca, 'LineWidth', 2, fset.sfont{:}, fset.axOpts{:});
%     xlabel('Lifetime (s)', fset.lfont{:});
%     ylabel('Frequency', fset.lfont{:});
%     
%     
%     
%     
%     hl = legend(hp, legendText{:}, 'Location', 'NorthEast');
%     set(hl, 'Box', 'off', fset.ifont{:});
elseif iscell(lftRes)
    if nargin<2
        colorA = hsv2rgb([0.6 1 1;
            1/3 1 1;
            0 1 1;
            1/7 1 1;
            5/9 1 1]);
    end
    colorB = rgb2hsv(colorA);
    colorB(:,2) = 0.5;
    colorB = hsv2rgb(colorB);
    
    nd = numel(lftRes);
    h = figure;
    hold on;
    for i = 1:nd
        hp(2*(i-1)+2) = plot(lftRes{i}.t, lftRes{i}.meanLftHist_B, '.-', 'Color', colorB(i,:), 'LineWidth', 2, 'MarkerSize', 16);
        hp(2*(i-1)+1) = plot(lftRes{i}.t, lftRes{i}.meanLftHist_A, '.-', 'Color', colorA(i,:), 'LineWidth', 2, 'MarkerSize', 16);
        expName = getDirFromPath(getExpDir(lftRes{i}.data));
        legendText{2*(i-1)+1} = [expName ', above threshold (' num2str(mean(lftRes{i}.pctAbove)*100,'%.1f') ' ± ' num2str(std(lftRes{i}.pctAbove)*100,'%.1f') ' %)'];
        legendText{2*(i-1)+2} = [expName ', below threshold (' num2str(mean(1-lftRes{i}.pctAbove)*100,'%.1f') ' ± ' num2str(std(lftRes{i}.pctAbove)*100,'%.1f') ' %)'];
    end
    axis([0 min(120, lftRes{1}.t(end)) 0 0.05]);
    set(gca, 'LineWidth', 2, fset.sfont{:}, fset.axOpts{:});
    xlabel('Lifetime (s)', fset.lfont{:});
    ylabel('Frequency', fset.lfont{:});
    hl = legend(hp, legendText{:}, 'Location', 'NorthEast');
    set(hl, 'Box', 'off', fset.ifont{:}, 'Interpreter', 'none');
    
else
    error('Incompatible input');
end
