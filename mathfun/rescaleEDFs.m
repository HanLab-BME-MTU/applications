%[a c medIdx] = rescaleEDFs(samples, varargin) computes the x-scaling factor between the EDFs of the input sample sets
%
% Outputs:
%          a : scaling factor
%          c : estimated fraction of missing data
%

% Francois Aguet, 03/06/2012

function [a c refIdx] = rescaleEDFs(samples, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Display', false, @islogical);
ip.addParamValue('Reference', 'med', @(x) any(strcmpi(x, {'max', 'med'})));
ip.addParamValue('FigureName', 'EDF scaling');
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
    f_edf = cell(1,nd);
    x_edf = cell(1,nd);
    f = cell(1,nd);
    
    x = 0:0.001:1;
    for i = 1:nd
        [f_edf{i}, x_edf{i}] = ecdf(samples{i});
        f_edf{i} = f_edf{i}(2:end);
        x_edf{i} = x_edf{i}(2:end);
        f{i} = interp1(f_edf{i}, x_edf{i}, x);
    end
    
    % scale to reference distribution, with offset for missing data
    switch ip.Results.Reference
        case 'max' % highest-valued (highest mean) distribution
            %mu = cellfun(@(i) mean(i), samples);
            mu = cellfun(@(i) median(i), samples);
            refIdx = find(mu==max(mu),1,'first');
        case 'med' % median distribution
            M = vertcat(f{:});
            medianEDF = median(M,1);
            J = nansum((M-repmat(medianEDF, [nd 1])).^2, 2);
            refIdx = find(J==min(J),1,'first');
    end
    idx = setdiff(1:nd, refIdx);
    
    x0 = linspace(0,max(vertcat(samples{:})),1000);
    
    % Generate EDFs
    fEDF = cell(1,nd);
    xEDF = cell(1,nd);
    for i = 1:nd
        [fEDF{i}, xEDF{i}] = ecdf(samples{i});
    end
    
    
    a = ones(1,nd);
    c = zeros(1,nd);
    refEDF = interpEDF(xEDF{refIdx}, fEDF{refIdx}, x0);
    for i = 1:nd-1
        p = lsqnonlin(@cost, [1 0], [0 -1], [Inf 1], opts, xEDF{idx(i)}, fEDF{idx(i)}, refEDF, x0);
        a(idx(i)) = p(1);
        c(idx(i)) = p(2);
    end


    if ip.Results.Display
        %colorV = rand(nd,3);
        colorV = zeros(nd,3);
        
        fset = loadFigureSettings('print');
        fset.axOpts = [fset.axOpts, {'TickLength', [0.015 0]}];
        T99 = prctile(samples{refIdx}, 99.9);
        lw = 1;

        %figure('Position', pos, 'PaperPositionMode', 'auto', 'Color', 'w', 'Name', ip.Results.FigureName);
        figure;
        axPos = fset.axPos;
        dx = 0.3*fset.axPos(4);
        axPos(2) = axPos(2)+dx+axPos(4);
        
        
        axes(fset.axOpts{:}, 'Position', axPos);
        hold on;
        for i = 1:nd-1
            plot(x_edf{idx(i)}, f_edf{idx(i)}, '-', 'Color', colorV(i,:), 'LineWidth', lw);
        end
        hp = plot(x_edf{refIdx}, f_edf{refIdx}, 'r', 'LineWidth', lw);

        axis([0 T99 0 1.01]);
        set(gca, 'YTick', 0:0.2:1, 'YTickLabel', ['0' arrayfun(@(x) num2str(x, '%.1f'), 0.2:0.2:1, 'UniformOutput', false)], 'XTickLabel', []);
        ylabel('Cumulative frequency', fset.lfont{:});
        text(0, 1.1, 'Raw EDF', 'HorizontalAlignment', 'left', fset.lfont{:});
        hl = legend(hp, 'Median distr.', 'Location', 'SouthEast');
        set(hl, 'Box', 'off', fset.sfont{:});
        
        
        axes(fset.axOpts{:});
        hold on;
        for i = 1:nd-1
            ci = c(idx(i));
            plot(x_edf{idx(i)}*a(idx(i)), ci+(1-ci)*f_edf{idx(i)}, 'Color', colorV(i,:), 'LineWidth', lw);
        end
        plot(x_edf{refIdx}, f_edf{refIdx}, 'r', 'LineWidth', lw);
        axis([0 T99 0 1.01]);
        set(gca, 'YTick', 0:0.2:1, 'YTickLabel', ['0' arrayfun(@(x) num2str(x, '%.1f'), 0.2:0.2:1, 'UniformOutput', false)]);
        xlabel('Max. fluo. intensity (A.U.)', fset.lfont{:});
        ylabel('Cumulative frequency', fset.lfont{:});
        text(0, 1.1, 'Scaled EDF', 'HorizontalAlignment', 'left', fset.lfont{:});
         
        
        % Plot inset with scales
        axPos = fset.axPos;
        axes(fset.axOpts{:}, 'Position', [axPos(1)+0.85*axPos(3) 1.2*axPos(2) axPos(3)/8 axPos(4)*0.75]);
        hold on;
        plot(zeros(numel(a)), a, 'o', 'Color', 0.4*[1 1 1], 'LineWidth', 1, 'MarkerSize', 5);
        he = errorbar(0, mean(a), std(a), 'Color', 0*[1 1 1], 'LineWidth', 1.5);
        plot(0.1*[-1 1], mean(a)*[1 1], 'Color', 0*[1 1 1], 'LineWidth', 1.5);
        setErrorbarStyle(he, 0.15);
        ymax = ceil(max(a)/0.2)*0.2;
        axis([-0.5 0.5 0 ymax]);
        ya = 0:0.2:ymax;
        set(gca, 'TickLength', fset.TickLength/0.75, 'XTick', [], 'YTick', ya, 'YTickLabel', ['0' arrayfun(@(x) num2str(x, '%.1f'), ya(2:end), 'UniformOutput', false)]);
        ylabel('Relative scale', fset.sfont{:});
    end
end


function v = cost(p, xEDF, fEDF, refEDF, x0)
a = p(1);
c = p(2);

f_i = interpEDF(xEDF, fEDF, x0/a);
v = c+(1-c)*f_i - refEDF;
v(f_i==0 | f_i==1 | refEDF==0 | refEDF==1) = 0;


%function f = interpEDF(samples, x)
function f = interpEDF(xEDF, fEDF, x)
f = interp1(xEDF(2:end), fEDF(2:end), x(2:end), 'linear');
f(1:find(~isnan(f),1,'first')-1) = 0;
f(find(isnan(f),1,'first'):end) = 1;
f = [0 f];
