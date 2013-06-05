%[mask, clusterID, nn] = cellMaskFromCCPDensity(data, varargin)

% Francois Aguet, May 2013

function [mask, clusterID, nn, R] = cellMaskFromCCPDensity(data, varargin)%frame, connect, showHist, modeRatio)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParamValue('Close', true, @islogical);
ip.addParamValue('MinSize', 5);
ip.addParamValue('Radius', []);
ip.addParamValue('Display', false, @islogical);
ip.parse(data, varargin{:});
R = ip.Results.Radius;

% if single input: data
tracks = loadTracks(data, 'Category', 'all', 'Cutoff_f', 5);
x = arrayfun(@(i) nanmean(i.x(1,:)), tracks);
y = arrayfun(@(i) nanmean(i.y(1,:)), tracks);

X = [x' y'];

if isempty(R)
    kdtreeobj = KDTree(X);
    np = size(X,1);
    dist = zeros(np,1);
    for i = 1:np
        [~,idist] = kdtreeobj.knn(X(i,:),2);
        dist(i) = idist(2);
    end
    [fEDF,xEDF] = ecdf(dist);

    opts = optimset('Jacobian', 'off', ...
        'MaxFunEvals', 1e4, ...
        'MaxIter', 1e4, ...
        'Display', 'off', ...
        'TolX', 1e-8, ...
        'Tolfun', 1e-8);
    
    p = lsqnonlin(@cost, [1 1], [0 0], [Inf Inf], opts, xEDF, fEDF);
    %m = chiCDF(xEDF, p(1), p(2));
    %figure; plot(xEDF,fEDF,'k');
    %hold on; plot(xEDF,m,'r');
    R = gammaincinv(0.95, p(1), 'lower')*p(2);
end
[clusterID, nn] = dbscan(X, ip.Results.MinSize, R);
ny = data.imagesize(1);
nx = data.imagesize(2);

mask = zeros(ny,nx);
mask(sub2ind([ny nx], round(y(clusterID~=0)), round(x(clusterID~=0)))) = 1;
mask = imclose(mask, strel('disk',round(R)));
mask = imdilate(mask, strel('disk',round(R/2)));

if ip.Results.Close
    mask = imfill(mask, 'holes');
end

if ip.Results.Display
    figure; imagesc(mask); axis image; colormap(gray(256)); colorbar;
    hold on;
    nc = numel(setdiff(unique(clusterID),0));
    plot(x(clusterID==0), y(clusterID==0), 'o', 'Color', 0.6*[1 1 1]);
    for c = 1:nc
        plot(x(clusterID==c), y(clusterID==c), 'o', 'Color', 'm');
    end
end



function v = cost(p, x, f)
v = chiCDF(x, p(1), p(2)) - f;

function f = chiCDF(x, k, s)
f = gammainc(0.5*x.^2/s^2, 0.5*k, 'lower');

