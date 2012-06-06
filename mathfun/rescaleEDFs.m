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
ip.addParamValue('Reference', 'max', @(x) any(strcmpi(x, {'max', 'med'})));
ip.parse(varargin{:});

nd = numel(samples);

if nd==1
    a = 1;
    c = 0;
    refIdx = 1;
else
    
    opts = optimset('Jacobian', 'off', ...
        'MaxFunEvals', 1e4, ...
        'MaxIter', 1e4, ...
        'Display', 'off', ...
        'TolX', 1e-8, ...
        'Tolfun', 1e-8);
    
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
    
    a = ones(1,nd);
    c = zeros(1,nd);
    for i = 1:nd-1
        p = lsqnonlin(@cost, [1 0], [0 0], [Inf 1], opts, samples{idx(i)}, samples{refIdx});
        a(idx(i)) = p(1);
        c(idx(i)) = p(2);
    end
    
    if ip.Results.Display
        fset = loadFigureSettings();
        T99 = prctile(samples{refIdx}, 99.9);

        pos = get(0, 'DefaultFigurePosition');
        pos(3) = 800;
        pos(4) = 400;
        figure('Position', pos, 'PaperPositionMode', 'auto', 'Color', 'w');
        
        axes('Units', 'pixels', 'Position', [80 80 300 280]);
        hold on;
        plot(x_edf{refIdx}, f_edf{refIdx}, 'r', 'LineWidth', 1.5);
        for i = 1:nd-1
            plot(x_edf{idx(i)}, f_edf{idx(i)}, 'k', 'LineWidth', 1);
        end
        axis([0 T99 0 1.01]);
        set(gca, fset.axOpts{:}, 'LineWidth', 2, fset.tfont{:});
        xlabel('Max. fluo. intensity (A.U.)', fset.sfont{:});
        ylabel('P(X \leq x)', fset.sfont{:});
        title('Raw EDF', fset.sfont{:});
        
        axes('Units', 'pixels', 'Position', [440 80 300 280]);
        hold on;
        plot(x_edf{refIdx}, f_edf{refIdx}, 'r', 'LineWidth', 1.5);
        for i = 1:nd-1
            ci = c(idx(i));
            %plot(p(i)*x_edf{idx(i)}, f_edf{idx(i)}, 'k', 'LineWidth', 1);
            plot(x_edf{idx(i)}*a(idx(i)), ci+(1-ci)*f_edf{idx(i)}, 'k', 'LineWidth', 1);
        end
        axis([0 T99 0 1.01]);
        set(gca, fset.axOpts{:}, 'LineWidth', 2, fset.tfont{:}, 'YTick', [], 'YColor', 'w');
        xlabel('Max. fluo. intensity (A.U.)', fset.sfont{:});
        title('Scaled EDF', fset.sfont{:});
        
        
        % Histograms
        figure('Position', pos, 'PaperPositionMode', 'auto', 'Color', 'w');
        
        axes('Units', 'pixels', 'Position', [80 80 300 280]);
        hold on;
        
        %dx = 10;
        %xi = 0:dx:x_edf{refIdx}(end);
        %ni = hist(samples{refIdx}, xi);
        %ni = ni/sum(ni)/dx;
        [ni,xi] = ksdensity(samples{refIdx}, 'npoints', 1000);
        %sum(ni)
        plot(xi, ni, 'r-', 'LineWidth', 1.5);
        for i = 1:nd-1
            %ni = hist(samples{idx(i)}, xi);
            %ni = ni/sum(ni)/dx;
            [ni,xi] = ksdensity(samples{idx(i)}, 'npoints', 1000);
            plot(xi, ni, 'k-', 'LineWidth', 1);
        end
        set(gca, fset.axOpts{:}, 'LineWidth', 2, fset.tfont{:}, 'XLim', [0 T99]);
        xlabel('Max. fluo. intensity (A.U.)', fset.sfont{:});
        ylabel('P(X \leq x)', fset.sfont{:});
        title('Raw kernel density', fset.sfont{:});
        YLim = get(gca, 'YLim');
        
        axes('Units', 'pixels', 'Position', [440 80 300 280]);
        hold on;
        %dx = 10;
        %xi = 0:dx:x_edf{refIdx}(end);
        %ni = hist(samples{refIdx}, xi);
        %ni = ni/sum(ni)/dx;
        [ni,xi] = ksdensity(samples{refIdx}, 'npoints', 1000);
        plot(xi, ni, 'r-', 'LineWidth', 1.5);
        for i = 1:nd-1
            %ni = hist(samples{idx(i)}*a(idx(i)), xi);
            %ni = ni/sum(ni)/dx * (1-c(idx(i))); 
            [ni,xi] = ksdensity(samples{idx(i)}*a(idx(i)), 'npoints', 1000);
            ni = ni*(1-c(idx(i)));
            plot(xi, ni, 'k-', 'LineWidth', 1);
        end
        axis([0 T99 YLim]);
        set(gca, fset.axOpts{:}, 'LineWidth', 2, fset.tfont{:}, 'YTick', [], 'YColor', 'w');
        xlabel('Max. fluo. intensity (A.U.)', fset.sfont{:});
        title('Scaled kernel density', fset.sfont{:});
        
    end
end



function v = cost(p, samples, refSamples)
a = p(1);
c = p(2);

x0 = linspace(0,max(refSamples),1000);

refEDF = interpEDF(refSamples, x0);

f_i = interpEDF(samples, x0/a);
v = c+(1-c)*f_i - refEDF;
v(f_i==0 | f_i==1) = 0;



% function v = cost(p, samples, refSamples)
% nd = numel(samples);
% a = p(1:2:end);
% c = p(2:2:end);
% 
% x0 = linspace(0,max(refSamples),1000);
% 
% refEDF = interpEDF(refSamples, x0);
% 
% % a = p(1);
% % c = p(2);
% % % f_i = interpEDF(samples{1}, x0/a);
% % v = c+(1-c)*f_i - refEDF;
% % v(f_i==0 | f_i==1) = 0;
% 
% v = cell(1,nd);
% for i = 1:nd
%     iSamples = samples{i};
%     f_i = interpEDF(iSamples, x0/a(i));
%     
%     v{i} = c(i)+(1-c(i))*f_i - refEDF;
%     v{i}(f_i==0 | f_i==1) = 0;
% end    
% v = vertcat(v{:});
% v(isnan(v)) = [];


function f = interpEDF(samples, x)
[fEDF, xEDF] = ecdf(samples);
f = interp1(xEDF(2:end), fEDF(2:end), x(2:end), 'linear');
f(1:find(~isnan(f),1,'first')-1) = 0;
f(find(isnan(f),1,'first'):end) = 1;
f = [0 f];




