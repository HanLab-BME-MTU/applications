% Francois Aguet, 06/08/2012

function outlierIdx = detectEDFOutliers(samples, varargin)

ns = numel(samples);

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('samples', @iscell);
ip.addOptional('offset', zeros(1,ns));
ip.addOptional('refIdx', []);
ip.addParamValue('Display', 'on', @(x) any(strcmpi(x, {'on', 'off'})));
ip.addParamValue('FigureName', '');
ip.parse(samples, varargin{:});
offset = ip.Results.offset;
medIdx = ip.Results.refIdx;
xEDF = cell(1,ns);
fEDF = cell(1,ns);
for i = 1:ns
    [fEDF{i}, xEDF{i}] = ecdf(samples{i});
end

% perform detection over values with full overlap only
minV = max(cellfun(@(i) min(i), xEDF));
maxV = min(cellfun(@(i) max(i), xEDF));
x = vertcat(xEDF{:});
x = unique(x(minV<=x & x<=maxV))';

for i = 1:ns
    fEDF{i} = interp1(xEDF{i}(2:end), offset(i) + (1-offset(i))*fEDF{i}(2:end), x);
end
M = vertcat(fEDF{:});

% median ECDF: majority of median indexes
% medIdx = sum(repmat(median(M,1), [ns 1])==M,2);
% medIdx = find(medIdx==max(medIdx),1,'first');
medianEDF = median(M,1);
if isempty(medIdx)
    J = nansum((M-repmat(medianEDF, [ns 1])).^2, 2);
    medIdx = find(J==min(J),1,'first');
end
medEDF = M(medIdx,:);

% MAD: median(abs(X-median(X)))
% sigma = 1/norminv(0.75) * mad(M, 1, 1);
xMAD = median(abs(M-repmat(medEDF,[ns 1])),1);
sigma = 1/norminv(0.75) * xMAD;

xi = linspace(x(1), x(end), 1000);
si = filterGauss1D(interp1(x, sigma, xi, 'cubic'), 15);
si = interp1(xi, si, x);

% 3*sigma bounds
ub = medEDF + 3*si;
lb = medEDF - 3*si;

% idx = find(chub==numel(x));
% chub = interp1(x(chub(idx:end)), ub(chub(idx:end)), x, 'cubic');
% du = chub - M(medIdx,:);
% ub = M(medIdx,:)+du;
% lb = M(medIdx,:)-du;

% ub(ub>1) = 1;
% lb(lb<0) = 0;

% for each EDF, test whether 95% of its points fall within median ± 3*sigma
outlierIdx = zeros(1,ns);
for i = 1:ns
    tmp = M(i,:)<lb | ub<M(i,:);
    outlierIdx(i) = sum(tmp)/numel(tmp) > 0.05;    
end
outlierIdx = find(outlierIdx);

if ip.Results.Display
    fset = loadFigureSettings();
    hp = NaN(1,4);
    figure('Name', ip.Results.FigureName);
    hold on;
    hp(4) = fill([x x(end:-1:1)], [lb ub(end:-1:1)], 0.6*[1 1 1], 'EdgeColor', 'none');
    
    if ~isempty(outlierIdx)
        hx = plot(x, M(outlierIdx,:), 'Color', [1 0 0], 'LineWidth', 2);
        hp(2) = hx(1);
    end
    hx = plot(x, M(setdiff(1:ns, [outlierIdx medIdx]),:), 'Color', [0 1 0], 'LineWidth', 2);
    hp(1) = hx(1);
    hp(3) = plot(x, medEDF, 'k', 'LineWidth', 2);

    set(gca, fset.sfont{:}, fset.axOpts{:},...
        'YTick', 0:0.1:1, 'YLim', [0 1.01]);
    xlabel('Time (s)');
    ylabel('F(t)');
    if isnan(hp(2))
        hp(2) = [];
        hl = legend(hp, ' Inlier distribution', ' Median distribution', ' Critical area', 'Location', 'SouthEast');
    else
        hl = legend(hp, ' Inlier distribution', ' Outlier distribution', ' Median distribution', ' Critical area', 'Location', 'SouthEast');
    end
    set(hl, fset.tfont{:}, 'Box', 'off');
    pos = get(hl, 'Position');
    pos(2) = 0.2;
    set(hl, 'Position', pos);
    drawnow;
end

