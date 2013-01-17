%[a c medIdx] = rescaleEDFs(samples, varargin) computes the x-scaling factor between the EDFs of the input sample sets
%
% Outputs:
%          a : scaling factor
%          c : estimated fraction of missing data
%

% Francois Aguet, 03/06/2012

function [a, c, refIdx] = rescaleEDFs(samples, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Display', false, @islogical);
ip.addParamValue('Reference', 'med', @(x) isscalar(x) || any(strcmpi(x, {'max', 'med'})));
ip.addParamValue('FigureName', 'EDF scaling');
ip.addParamValue('XTick', []);
ip.parse(varargin{:});

nd = numel(samples);
samples = cellfun(@(i) i(:), samples, 'UniformOutput', false);

if nd==1
    a = 1;
    c = 0;
    refIdx = 1;
else
    
    opts = optimset('Jacobian', 'off', ...
        'MaxFunEvals', 1e4, ...
        'MaxIter', 1e4, ...
        'Display', 'off', ...
        'TolX', 1e-6, ...
        'Tolfun', 1e-6);
    
    % Generate EDF for each set of samples
    fEDF = cell(1,nd);
    xEDF = cell(1,nd);
    
    for i = 1:nd
        [fEDF{i}, xEDF{i}] = ecdf(samples{i});
        %fEDF{i} = fEDF{i}(2:end);
        %xEDF{i} = xEDF{i}(2:end);
    end
    
    % Now, interpolate on f(x)
    fi = 0:0.001:1;
    f = cell(1,nd);
    for i = 1:nd
        f{i} = interp1(fEDF{i}, xEDF{i}, fi);
    end
    
    % scale to reference distribution, with offset for missing data
    switch ip.Results.Reference
        case 'max' % highest-valued (highest median) distribution
            mu = cellfun(@(i) median(i), samples);
            refIdx = find(mu==max(mu),1,'first');
        case 'med' % median distribution
            M = vertcat(f{:});
            medianEDF = median(M,1);
            J = nansum((M-repmat(medianEDF, [nd 1])).^2, 2);
            refIdx = find(J==min(J),1,'first');
        otherwise
            refIdx = ip.Results.Reference;
    end
    idx = setdiff(1:nd, refIdx);
    
    %x0 = linspace(0,max(vertcat(samples{:})),1000);
    %x0 = linspace(min(samples{refIdx}),max(samples{refIdx}),1000); % robust
    x0 = linspace(prctile(vertcat(samples{:}),1), prctile(vertcat(samples{:}),99), 1000);
    
    %     % Generate EDFs
    %     fEDF = cell(1,nd);
    %     xEDF = cell(1,nd);
    %     for i = 1:nd
    %         [fEDF{i}, xEDF{i}] = ecdf(samples{i});
    %         %fEDF{i} = fEDF{i}(2:end);
    %         %xEDF{i} = xEDF{i}(2:end);
    %     end
    
    
    a = ones(1,nd);
    c = zeros(1,nd);
    refEDF = interpEDF(xEDF{refIdx}, fEDF{refIdx}, x0);
    %refEDF = interpEDF(x, medEDF, x0);
    for i = 1:nd-1
        p = lsqnonlin(@cost, [1 0], [0 -1], [Inf 1], opts, xEDF{idx(i)}, fEDF{idx(i)}, refEDF, x0);
        a(idx(i)) = p(1);
        c(idx(i)) = p(2);
    end
    
    [~,idxa] = sort(a);
    [~,idxa] = sort(idxa);
    idxa(refIdx) = []; % reference shown in black

    if ip.Results.Display
        colorV = hsv(nd);
        fset = loadFigureSettings('print');
        if isempty(ip.Results.XTick)
            T99 = prctile(samples{refIdx}, 99.9);
            xa = 0:50:T99;
        else
            xa = ip.Results.XTick;
            T99 = xa(end);
        end
        lw = 1;
        axPos = fset.axPos;
        dx = 0.3*fset.axPos(4);
        
        figure(fset.fOpts{:}, 'Position', [5 5 8 1.5+axPos(4)*2+dx+1], 'Color', 'w', 'Name', ip.Results.FigureName);
        axPos(2) = axPos(2)+dx+axPos(4);
        axes(fset.axOpts{:}, 'Position', axPos);
        hold on;
        for i = nd-1:-1:1
            plot(xEDF{idx(i)}, fEDF{idx(i)}, '-', 'Color', colorV(idxa(i),:), 'LineWidth', lw);
        end
        hp = plot(xEDF{refIdx}, fEDF{refIdx}, 'k', 'LineWidth', lw+0.5);
        axis([0 T99 0 1.01]);
        set(gca, 'YTick', 0:0.2:1, 'XTick', xa, 'XTickLabel', []);
        formatTickLabels(gca);
        ylabel('Cumulative frequency', fset.lfont{:});
        text(0, 1.1, 'Raw distributions', 'HorizontalAlignment', 'left', fset.lfont{:});
        hl = legend(hp, ' Median distr.', 'Location', 'SouthEast');
        set(hl, 'Box', 'off', fset.sfont{:}, 'Position', [5 6 1.5 1]);
        
        axes(fset.axOpts{:});
        hold on;
        plot(xEDF{refIdx}, fEDF{refIdx}, 'k', 'LineWidth', lw+0.5);
        for i = nd-1:-1:1
            ci = c(idx(i));
            plot(xEDF{idx(i)}*a(idx(i)), ci+(1-ci)*fEDF{idx(i)}, 'Color', colorV(idxa(i),:), 'LineWidth', lw);
        end
        axis([0 T99 0 1.01]);
        set(gca, 'YTick', 0:0.2:1, 'XTick', xa);
        formatTickLabels(gca);
        xlabel('Max. fluo. intensity (A.U.)', fset.lfont{:});
        ylabel('Cumulative frequency', fset.lfont{:});
        text(0, 1.1, 'Scaled distributions', 'HorizontalAlignment', 'left', fset.lfont{:});
        
        
        % Plot inset with scales
        axPos = fset.axPos;
        axes(fset.axOpts{:}, 'Position', [axPos(1)+0.85*axPos(3) 1.5*axPos(2) axPos(3)/8 axPos(4)*0.5], 'Layer', 'bottom');
        hold on;
        plot(zeros(numel(a)), a, 'o', 'Color', 0.4*[1 1 1], 'LineWidth', 1, 'MarkerSize', 5);
        he = errorbar(0, mean(a), std(a), 'Color', 0*[1 1 1], 'LineWidth', 1.5);
        plot(0.1*[-1 1], mean(a)*[1 1], 'Color', 0*[1 1 1], 'LineWidth', 1.5);
        setErrorbarStyle(he, 0.15);
        ylim = [floor(min(a)/0.2) ceil(max(a)/0.2)]*0.2;
        axis([-0.5 0.5 ylim]);
        ya = linspace(ylim(1), ylim(2), 5);
        set(gca, 'TickLength', fset.TickLength*3, 'XTick', [], 'YTick', ya, 'XColor', 'w');
        formatTickLabels(gca);
        ylabel('Relative scale', fset.sfont{:});
    end
end


function v = cost(p, xEDF, fEDF, refEDF, x0)
a = p(1);
c = p(2);

f_i = interpEDF(xEDF, fEDF, x0/a);
v = c+(1-c)*f_i - refEDF;
v(f_i==0 | f_i==1 | refEDF==0 | refEDF==1) = 0;


function f = interpEDF(xEDF, fEDF, x)
f = interp1(xEDF(2:end), fEDF(2:end), x(2:end), 'linear');
f(1:find(~isnan(f),1,'first')-1) = 0;
f(find(isnan(f),1,'first'):end) = 1;
f = [0 f];
