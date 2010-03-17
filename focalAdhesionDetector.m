function [params, Im] = focalAdhesionDetector(I, mask, sigmaPSF)

% make sure I is type double
I = double(I);

% Get a first crude segmentation
BW = logical(blobSegmentThreshold(I,1,0));

% Restrict analysis to cell footprint
BW = logical(BW .* mask);

% Remove noise by filtering image with a Gaussian whose sigma = sigmaPSF
% NOTE: this is not optimal since this step has already been done in the
% blobSegmentThreshold function.
If = Gauss2D(I,sigmaPSF,1);

% Estimate background by filtering image with a Gaussian whose sigma = 10
% NOTE: this is not optimal since this step has already been done in the
% blobSegmentThreshold function.
BG = Gauss2D(I,10,1);

% Update the image
Iu = If - BG;

% Get initial rod parameter from the the connected component properties
CCstats = regionprops(BW, Iu, 'Centroid','Orientation','MajorAxisLength','PixelValues');

% params is a Nx5 matrix where each column stands for Xc, Yc, A, l, theta.
params = zeros(numel(CCstats), 5);
params(:,1:2) = vertcat(CCstats(:).Centroid);
params(:,3) = cellfun(@(px) mean(px), {CCstats(:).PixelValues});
params(:,4) = vertcat(CCstats(:).MajorAxisLength);
params(:,5) = -vertcat(CCstats(:).Orientation) * pi/180;

% Fit
initParams = params;
options = optimset('Jacobian', 'on', 'MaxFunEvals', 1e4, 'MaxIter', 1e4, ...
    'Display', 'off', 'TolX', 1e-6, 'Tolfun', 1e-6);
fun = @(x) fitDLSegment2D(x, Iu, sigmaPSF);
params = lsqnonlin(fun, initParams, [], [], options);

imshow(Iu, []); hold on;
quiver(initParams(:,1), initParams(:,2), (initParams(:,4) / 2) .* cos(initParams(:,5)), ...
    (initParams(:,4) / 2) .* sin(initParams(:,5)), 0, 'g');
quiver(initParams(:,1), initParams(:,2), (initParams(:,4) / 2) .* cos(initParams(:,5) + pi), ...
    (initParams(:,4) / 2) .* sin(initParams(:,5) + pi), 0, 'g');
plot(initParams(:,1), initParams(:,2), 'r.');

quiver(params(:,1), params(:,2), (params(:,4) / 2) .* cos(params(:,5)), ...
    (params(:,4) / 2) .* sin(params(:,5)), 0, 'g');
quiver(params(:,1), params(:,2), (params(:,4) / 2) .* cos(params(:,5) + pi), ...
    (params(:,4) / 2) .* sin(params(:,5) + pi), 0, 'g');
plot(params(:,1), params(:,2), 'r.'); hold off;

function [F J] = fitDLSegment2D(x, Iu, sigmaPSF)

F = reshape(Iu - imageSegmentModel(size(Iu), sigmaPSF, x), numel(Iu), 1);

if nargout > 1
    [n p] = size(x);
    m = numel(Iu);
    
    indPixels = cell(n,1);
    indParams = cell(n,1);
    val = cell(n,1);    
    
    for i = 1:n
        % Define the support of segment i
        % NOTE 1: these are the same ranges calculated in imageSegmentModel.
        % Make this computation only once.
        % NOTE 2: IMPORTANT: check in Mathematica that the support of the
        % partial derivatives of the segment is the same as for the segment
        % itself.
        xRange = [];
        yRange = [];
        
        % Compute all partial derivatives of segment parameters (except
        % sigmaPSF, since it is fixed).
        dFdXc = dlSegment2D_dFdXc(xRange, yRange, x(i,1), x(i,2), x(i,3), sigmaPSF, x(i,4), x(i,5));
        dFdYc = dlSegment2D_dFdYc(xRange, yRange, x(i,1), x(i,2), x(i,3), sigmaPSF, x(i,4), x(i,5));
        dFdA = dlSegment2D_dFdA(xRange, yRange, x(i,1), x(i,2), x(i,3), sigmaPSF, x(i,4), x(i,5));
        dFdl = dlSegment2D_dFdl(xRange, yRange, x(i,1), x(i,2), x(i,3), sigmaPSF, x(i,4), x(i,5));
        dFdt = dlSegment2D_dFdt(xRange, yRange, x(i,1), x(i,2), x(i,3), sigmaPSF, x(i,4), x(i,5));
        
        [X Y] = meshgrid(xRange, yRange);
        ind = sub2ind(size(Iu), Y(:), X(:));
        
        indPixels{i} = repmat(ind, p, 1);
        indParams{i} = vertcat(arrayfun(@(k) ones(numel(ind), 1) * i + ...
            k * n, 0:p-1, 'UniformOutput', 'false'));        
        val{i} = vertcat(dFdXc(:), dFdYc(:), dFdA(:), dFdl(:), dFdt(:));
    end
    
    indPixels = vertcat(indPixels{:});
    indParams = vertcat(indParams{:});
    val = vertcat(val{:});
    J = sparse(indPixels, indParams, val, m, n * p, length(val));
end

