% Francois Aguet, 02/17/2012

function mask = maskFromFirstMode(img, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('img');
ip.addParamValue('Connect', true, @islogical);
ip.addParamValue('Display', false, @islogical);
ip.addParamValue('ModeRatio', 0.6, @isscalar);
ip.parse(img, varargin{:});

[ny,nx] = size(img);

g = filterGauss2D(img, 5);

v = g(:);
v(isnan(v)) = [];

pct = prctile(v, [1 99]);
v(v<pct(1) | pct(2)<v) = [];

[f,xi] = ksdensity(v, 'npoints', 100);

% local max/min
lmax = locmax1d(f, 3);
lmin = locmin1d(f, 3);

% max value
% hmax = find(f==max(f), 1, 'first');
dxi = xi(2)-xi(1);

% identify min after first mode
if ~isempty(lmin)
    idx = find(lmin>lmax(1), 1, 'first');
    if ~isempty(idx) && sum(f(1:lmin(idx(1))))*dxi < ip.Results.ModeRatio
        min0 = lmin(idx);
        T = xi(min0);
        mask = g>T;
        
        % retain largest connected component
        if ip.Results.Connect
            CC = bwconncomp(mask, 8);
            compsize = cellfun(@(i) numel(i), CC.PixelIdxList);
            mask = zeros([ny,nx]);
            mask(CC.PixelIdxList{compsize==max(compsize)}) = 1;
        end
    else
        mask = ones(ny,nx);
    end
else
    mask = ones(ny,nx);
end

if ip.Results.Display
    dx = xi(2)-xi(1);
    ni = hist(v, xi);
    ni = ni/sum(ni)/dx;
    
    figure;
    hold on;
    plot(xi, ni, 'k.-', 'LineWidth', 3, 'MarkerSize', 20);
    plot(xi, f, 'r-', 'LineWidth', 1);
    set(gca, 'LineWidth', 2, 'FontSize', 18);
end