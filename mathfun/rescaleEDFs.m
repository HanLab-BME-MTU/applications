%function [a medIdx] = rescaleEDFs(samples, varargin) computes the x-scaling factor between the EDFs of the input sample sets

% Francois Aguet, 03/06/2012

function [a medIdx] = rescaleEDFs(samples, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Display', false, @islogical);
ip.parse(varargin{:});

nd = numel(samples);

if nd==1
    a = 1;
    medIdx = 1;
else
    
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
    M = vertcat(f{:});
    medianEDF = median(M,1);
    J = nansum((M-repmat(medianEDF, [nd 1])).^2, 2);
    medIdx = find(J==min(J),1,'first');
    
    
    opts = optimset('Jacobian', 'off', ...
        'MaxFunEvals', 1e4, ...
        'MaxIter', 1e4, ...
        'Display', 'off', ...
        'TolX', 1e-8, ...
        'Tolfun', 1e-8);
    
    T99 = prctile(samples{medIdx}, 99.9);
    
    lb = zeros(1,nd-1);
    ub = inf(1,nd-1);
    idx = setdiff(1:nd, medIdx);
    [p,resnorm,~,~,~,~,J] = lsqnonlin(@cost, ones(1,nd-1), lb, ub, opts, x_edf(idx), f_edf(idx), x_edf{medIdx}, f_edf{medIdx}, T99);
    a = ones(1,nd);
    a(idx) = p;
    
    if ip.Results.Display
        fset = loadFigureSettings();
        
        pos = get(0, 'DefaultFigurePosition');
        pos(3) = 800;
        pos(4) = 400;
        figure('Position', pos, 'PaperPositionMode', 'auto', 'Color', 'w');
        
        axes('Units', 'pixels', 'Position', [80 80 300 280]);
        hold on;
        plot(x_edf{medIdx}, f_edf{medIdx}, 'r', 'LineWidth', 3);
        for i = 1:nd-1
            plot(x_edf{idx(i)}, f_edf{idx(i)}, 'k', 'LineWidth', 1);
        end
        axis([0 T99 0 1.01]);
        set(gca, 'LineWidth', 2, 'TickDir', 'out', fset.tfont{:});
        xlabel('Max. fluo. intensity (A.U.)', fset.sfont{:});
        ylabel('P(X \leq x)', fset.sfont{:});
        title('Raw EDF', fset.sfont{:});
        
        axes('Units', 'pixels', 'Position', [440 80 300 280]);
        hold on;
        plot(x_edf{medIdx}, f_edf{medIdx}, 'r', 'LineWidth', 3);
        for i = 1:nd-1
            plot(p(i)*x_edf{idx(i)}, f_edf{idx(i)}, 'k', 'LineWidth', 1);
        end
        axis([0 T99 0 1.01]);
        set(gca, 'LineWidth', 2, 'TickDir', 'out', fset.tfont{:}, 'YTick', [], 'YColor', 'w');
        xlabel('Max. fluo. intensity (A.U.)', fset.sfont{:});
        title('Scaled EDF', fset.sfont{:});
        
        
        % Histograms
        figure('Position', pos, 'PaperPositionMode', 'auto', 'Color', 'w');
        
        axes('Units', 'pixels', 'Position', [80 80 300 280]);
        hold on;
        
        dx = 10;
        xi = 0:dx:x_edf{medIdx}(end);
        ni = hist(samples{medIdx}, xi);
        ni = ni/sum(ni)/dx;
        %[ni,xi] = ksdensity(samples{medIdx}, 'npoints', 1000);
        plot(xi, ni, 'r-', 'LineWidth', 3);
        for i = 1:nd-1
            ni = hist(samples{idx(i)}, xi);
            ni = ni/sum(ni)/dx;
            %[ni,xi] = ksdensity(samples{idx(i)}, 'npoints', 1000);
            plot(xi, ni, 'k-', 'LineWidth', 1);
        end
        axis([0 T99 0 0.02]);
        set(gca, 'LineWidth', 2, 'TickDir', 'out', fset.tfont{:});
        xlabel('Max. fluo. intensity (A.U.)', fset.sfont{:});
        ylabel('P(X \leq x)', fset.sfont{:});
        title('Raw histogram', fset.sfont{:});
        
        axes('Units', 'pixels', 'Position', [440 80 300 280]);
        hold on;
        ni = hist(samples{medIdx}, xi);
        ni = ni/sum(ni)/dx;
        %[ni,xi] = ksdensity(samples{medIdx}, 'npoints', 1000);
        plot(xi, ni, 'r-', 'LineWidth', 3);
        for i = 1:nd-1
            ni = hist(samples{idx(i)}*p(i), xi);
            ni = ni/sum(ni)/dx;
            %[ni,xi] = ksdensity(samples{idx(i)}*p(i), 'npoints', 1000);
            plot(xi, ni, 'k-', 'LineWidth', 1);
        end
        axis([0 T99 0 0.02]);
        set(gca, 'LineWidth', 2, 'TickDir', 'out', fset.tfont{:}, 'YTick', [], 'YColor', 'w');
        xlabel('Max. fluo. intensity (A.U.)', fset.sfont{:});
        title('Scaled histogram', fset.sfont{:});
        
    end
end


function v = cost(p, x_edf, f_edf, xRef, fRef, T)
nd = numel(p);

v = cell(1,2);
for i = 1:nd
    v{i} = xRef - p(i)*interp1(f_edf{i}, x_edf{i}, fRef);
    %v{i} = p(i)*x_edf{i} - interp1(fRef, xRef, f_edf{i});
    v{i}(xRef>T) = []; % ignore differences above Tth percentile
end
v = vertcat(v{:});
v(isnan(v)) = [];
