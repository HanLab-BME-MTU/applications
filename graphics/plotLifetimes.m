function h = plotLifetimes(lftRes, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('DisplayMode', '');
ip.addParamValue('ShowExpFits', false, @islogical);
ip.addParamValue('ShowStatistics', false, @islogical);
ip.addParamValue('ShowCargoDependent', false, @islogical);
ip.addParamValue('CargoName', 'cargo');
ip.addParamValue('YTick', 0:0.01:0.04);
ip.parse(varargin{:});
ya = ip.Results.YTick;
lw = 0.75;
h = [];

fmt = ['%.' num2str(ceil(abs(log10(ya(2))))) 'f'];
yal = ['0' arrayfun(@(x) num2str(x, fmt), ya(2:end), 'UniformOutput', false)];

fset = loadFigureSettings(ip.Results.DisplayMode);
ce = [0 0 0.6;
    1/3 0.3 0.9;
    1/3 1 0.9];
ce = hsv2rgb(ce);

if isstruct(lftRes)
    
    h(1) = figure(fset.fOpts{:}, 'Name', 'Lifetime dist. (intensity threshold)');
    axes(fset.axOpts{:});
    hold on;
    if isfield(lftRes, 'pctVisit')
        hp(4) = plot(lftRes.t, mean(lftRes.pctVisit)*lftRes.meanLftHist_V, '-', 'Color', fset.cfB, 'LineWidth', lw);
    end
    hp(3) = plot(lftRes.t, mean(vertcat(lftRes.lftHist_Ia), 1), 'Color', ce(1,:), 'LineWidth', lw);
    hp(2) = plot(lftRes.t, mean(lftRes.pctBelow)*lftRes.meanLftHist_B, '-', 'Color', ce(2,:), 'LineWidth', lw);
    
    if ip.Results.ShowExpFits
        ff = mean(lftRes.pctBelow)*lftRes.meanLftHist_B;
        %ff = mean(vertcat(lftRes.lftHist_Ia), 1);
        [mu,~,Aexp,expF] = fitExpToHist(lftRes.t(5:end), ff(5:end));
        tx = 0:0.1:lftRes.t(end);
        plot(tx, Aexp/mu*exp(-1/mu*tx), 'r--', 'LineWidth', 1)
        
        %ff = mean(lftRes.pctBelow)*lftRes.meanLftHist_B + mean(lftRes.pctVisit)*lftRes.meanLftHist_V;
        %plot(lftRes.t, ff, '--', 'Color', 'k', 'LineWidth', 2);
        %[mu,~,Aexp,expF] = fitExpToHist(lftRes.t, ff);
        %plot(tx, Aexp/mu*exp(-1/mu*tx), 'm--', 'LineWidth', 1)
    end
    hp(1) = plot(lftRes.t, mean(lftRes.pctAbove)*lftRes.meanLftHist_A, '-', 'Color', ce(3,:), 'LineWidth', lw+0.5);

    
    axis([0 min(120, lftRes.t(end)) 0 ya(end)]);
    set(gca, 'XTick', 0:20:200, 'YTick', ya, 'YTickLabel', yal);
    xlabel('Lifetime (s)', fset.lfont{:});
    ylabel('Frequency', fset.lfont{:});
    
    ltext = {[' Above threshold: ' num2str(mean(lftRes.pctAbove)*100, '%.1f') ' ± ' num2str(std(lftRes.pctAbove)*100, '%.1f') ' %'],...
        [' Below threshold: ' num2str(mean(lftRes.pctBelow)*100, '%.1f') ' ± ' num2str(std(lftRes.pctBelow)*100, '%.1f') ' %'],...
        ' Raw distribution'};
    lheight = 1;
    if isfield(lftRes, 'pctVisit')
        ltext = [ltext(1:2) [' Visitors: ' num2str(mean(lftRes.pctVisit)*100, '%.1f') ' ± ' num2str(std(lftRes.pctVisit)*100, '%.1f') ' %'] ltext(3)];
        hp = hp([1 2 4 3]);
        lheight = 1.25;
    end
    hl = legend(hp, ltext{:});
    set(hl, 'Box', 'off', fset.tfont{:});
    if strcmpi(ip.Results.DisplayMode, 'print')
        set(hl, 'Position', [4 5-lheight 1.75 lheight]); 
    end
    %%
    if ip.Results.ShowStatistics
        fs = loadFigureSettings('print');
        h(2) = figure(fs.fOpts{:}, 'Position', [5 5 5 6.5]);
        axes(fs.axOpts{:}, 'Position', [1.5 2 3 4]);
        boxplot2(lftRes.stats, 'AdjustFigure', false, 'XLabels', {'Raw', 'Below thresh', 'Above thresh'},...
            'FaceColor', ce, 'BarWidth', 0.6, 'LineWidth', 1);
        ylabel('Lifetime (s)', fset.sfont{:});
    end

    %%
    if ip.Results.ShowCargoDependent && isfield(lftRes, 'lftHist_Apos')
        figure(fset.fOpts{:}, 'Name', 'Lifetime dist.');
        axes(fset.axOpts{:});
        hold on;
        
        pAS = mean(lftRes.pctAboveSignificant);
        pANS = mean(lftRes.pctAboveNotSignificant);
        pBS = mean(lftRes.pctBelowSignificant);
        pBNS = mean(lftRes.pctBelowNotSignificant);
        
        % normalize w/o visitors (temporary fix!)
        tmp = pAS+pANS+pBS+pBNS;
        pAS = pAS/tmp;
        pANS = pANS/tmp;
        pBS = pBS/tmp;
        pBNS = pBNS/tmp;
                
        hp = zeros(1,7);
        % total distr
        mu = mean(vertcat(lftRes.lftHist_Ia), 1);
        %bAll = 1.96*getSE(lftRes, 'lftHist_Ia');
        %fill([lftRes.t lftRes.t(end:-1:1)], [mu+bAll mu(end:-1:1)-bAll(end:-1:1)], 'r', 'EdgeColor', 'none');
        hp(1) = plot(lftRes.t, mu, 'Color', 0.6*[1 1 1], 'LineWidth', lw);
        
        % Cargo-negative distributions
        hueCN = 0.6;
        hp(5) = plot(lftRes.t, (pANS+pBNS)*mean(vertcat(lftRes.lftHist_neg), 1), 'Color', hsv2rgb([hueCN 0.7 1]), 'LineWidth', lw);
        hp(7) = plot(lftRes.t, pBNS*mean(lftRes.lftHist_Bneg,1), '-', 'Color', hsv2rgb([hueCN-0.1 0.7 0.9]), 'LineWidth', lw);
        mu = mean(lftRes.lftHist_Aneg,1);
        %bAll = 1.96*getSE(lftRes, 'lftHist_Aneg');
        %fill([lftRes.t lftRes.t(end:-1:1)], pANS*[mu+bAll mu(end:-1:1)-bAll(end:-1:1)], 'r', 'EdgeColor', 'none');
        
        
        % Cargo-positive distributions
        hueCP = 0.33;
        hp(2) = plot(lftRes.t, (pAS+pBS)*mean(vertcat(lftRes.lftHist_pos), 1), 'Color', hsv2rgb([hueCP 0.7 1]), 'LineWidth', lw);
        hp(4) = plot(lftRes.t, pBS*mean(lftRes.lftHist_Bpos,1), '-',  'Color', hsv2rgb([hueCP-0.1 0.7 0.9]), 'LineWidth', lw);
        
        
        hp(6) = plot(lftRes.t, pANS*mu, '-', 'Color', hsv2rgb([hueCN 1 0.9]), 'LineWidth', lw+0.5);
        hp(3) = plot(lftRes.t, pAS*mean(lftRes.lftHist_Apos,1), '-', 'Color', hsv2rgb([hueCP 1 0.9]), 'LineWidth', lw+0.5);
        
        % All, above/below threshold
        %hp(2) = plot(lftRes.t, mean(lftRes.pctBelow)*lftRes.meanLftHist_B, '-', 'Color', hsv2rgb([2/3 0.3 0.9]), 'LineWidth', 2);
        %hp(1) = plot(lftRes.t, mean(lftRes.pctAbove)*lftRes.meanLftHist_A, '-', 'Color', hsv2rgb([2/3 1 0.9]), 'LineWidth', 2);
        cargo = ip.Results.CargoName;
        fmt = '%.1f';
        legendText = {' All',...
            [' + ' cargo ' (' num2str((pAS+pBS)*100, fmt) '%)'],...
            [' + ' cargo ', max. int. > T (' num2str(pAS*100, fmt) '%)'],...
            [' + ' cargo ', max. int. < T (' num2str(pBS*100, fmt) '%)'],...
            [' - ' cargo ' (' num2str((pANS+pBNS)*100, fmt) '%)'],...
            [' - ' cargo ', max. int. > T (' num2str(pANS*100, fmt) '%)'],...
            [' - ' cargo, ', max. int. < T (' num2str(pBNS*100, fmt) '%)']};
        
        axis([0 min(120, lftRes.t(end)) 0 ya(end)]);
        set(gca, 'XTick', 0:20:200, 'YTick', ya, 'YTickLabel', yal);
        
        xlabel('Lifetime (s)', fset.lfont{:});
        ylabel('Frequency', fset.lfont{:});
        
        hl = legend(hp, legendText{:}, 'Location', 'NorthEast');
        set(hl, 'Box', 'off', fset.tfont{:}, 'Position', [4 3 2.5 2]);
        
        
        figure(fset.fOpts{:}, 'Name', 'Lifetime dist.');
        axes(fset.axOpts{:});
        hold on;
    
        hp = zeros(1,4);
        % total distr
        mu = mean(vertcat(lftRes.lftHist_Ia), 1);
        hp(1) = plot(lftRes.t, mu, 'Color', 0.6*[1 1 1], 'LineWidth', lw);
        hp(2) = plot(lftRes.t, pANS*mean(lftRes.lftHist_Aneg,1)+pAS*mean(lftRes.lftHist_Apos,1), '-', 'Color', hsv2rgb([hueCN 1 0]), 'LineWidth', lw);
        hp(4) = plot(lftRes.t, pANS*mean(lftRes.lftHist_Aneg,1), '-', 'Color', hsv2rgb([hueCN 1 0.9]), 'LineWidth', lw+0.5);
        hp(3) = plot(lftRes.t, pAS*mean(lftRes.lftHist_Apos,1), '-', 'Color', hsv2rgb([hueCP 1 0.9]), 'LineWidth', lw+0.5);
        
        % All, above/below threshold
        cargo = ip.Results.CargoName;
        fmt = '%.1f';
        legendText = {' All',...
            [' Max. int. > T (' num2str((pAS+pANS)*100, fmt) '%)'],...
            [' + ' cargo ', max. int. > T (' num2str(pAS*100, fmt) '%)'],...
            [' - ' cargo ', max. int. > T (' num2str(pANS*100, fmt) '%)']};
        
        axis([0 min(120, lftRes.t(end)) 0 ya(end)]);
        set(gca, 'XTick', 0:20:200, 'YTick', ya, 'YTickLabel', yal);
        
        xlabel('Lifetime (s)', fset.lfont{:});
        ylabel('Frequency', fset.lfont{:});
        
        hl = legend(hp, legendText{:}, 'Location', 'NorthEast');
        set(hl, 'Box', 'off', fset.tfont{:}, 'Position', [4 4 2.5 1.2]);

        
        
        
    end
  
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
    axis([0 min(120, lftRes{1}.t(end)) 0 ya(end)]);
    xlabel('Lifetime (s)', fset.lfont{:});
    ylabel('Frequency', fset.lfont{:});
    hl = legend(hp, legendText{:}, 'Location', 'NorthEast');
    set(hl, 'Box', 'off', fset.ifont{:}, 'Interpreter', 'none');
    
else
    error('Incompatible input');
end


function SE = getSE(lftData, fieldName)
M = lftData.(fieldName);
nd = size(M,1);
meanM = zeros(size(M));
for i = 1:nd
    meanM(i,:) = mean(M(setdiff(1:nd,i),:),1);
end
SE = std(meanM,[],1) / sqrt(nd);
