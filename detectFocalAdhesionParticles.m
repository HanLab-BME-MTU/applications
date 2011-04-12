function featuresInfo = detectFocalAdhesionParticles(ima, mask, sigmaPSF, kSigma, alpha, minDist)

ima = double(ima);
[nrows ncols] = size(ima);

% Filter image with laplacian
bandPassIso = filterLoG(ima,sigmaPSF);
bandPassIso(bandPassIso < 0) = 0;
bandPassIso(~mask) = 0;

% Filter image with steerable filter
[~,T] = steerableFiltering(ima,2,sigmaPSF);
 
% Compute the local maxima of the bandpass filtered images
locMaxIso = locmax2d(bandPassIso, [5 5]);
 
bw = blobSegmentThreshold(bandPassIso,0,0,mask);

locMaxIso(~bw) = 0;

indMax = find(locMaxIso);
[y x] = ind2sub(size(ima), indMax);

P = zeros(size(y, 1), 7);
P(:,1) = x;
P(:,2) = y;
P(:,3) = ima(indMax);
P(:,4:5) = sigmaPSF; % sigmaX, sigmaY
P(:,6) = T(indMax);

% Subresolution detection
radius = kSigma * sigmaPSF;
hside = ceil(radius);
xmin = x - hside;
xmax = x + hside;
ymin = y - hside;
ymax = y + hside;

isValid = find(xmin >= 1 & xmax <= ncols & ymin >= 1 & ymax <= nrows);

xmin = xmin(isValid);
xmax = xmax(isValid);
ymin = ymin(isValid);
ymax = ymax(isValid);
P = P(isValid,:);

stdP = zeros(size(P));
stdR = zeros(size(P,1),1);

[X,Y] = meshgrid(-hside:hside);
disk = X.^2 + Y.^2 - hside^2 <= 0;
numDegFreedom = nnz(disk) - 6;

success = false(numel(xmin),1);

for iFeature = 1:numel(xmin)
        
    crop = ima(ymin(iFeature):ymax(iFeature), xmin(iFeature):xmax(iFeature));
    P(iFeature,7) = min(crop(:)); % background
    P(iFeature,3) = P(iFeature,3) - P(iFeature,7); % amplitude above background
    crop(~disk) = NaN;
    
    [params, stdParams, ~, R] = fitAnisoGaussian2D(crop, ...
        [0, 0, P(iFeature,3), 3 * P(iFeature,4), P(iFeature,5), P(iFeature,6), P(iFeature,7)], 'xyArtC');
    
    % TEST 1: on parameter values:
    % - position must remain in the disk
    % - amplitude > 0
    % - sigmaX and sigmaY > 1
    % - sigmaX < 10 * sigmaPSF
    isValid = ...
        max(abs(params(1:2))) < radius && ...
        params(3) > 0 && ...
        min(params(4:5)) > 1 && ...
        params(4) < 10 * sigmaPSF;
    
    % TEST 2: model goodness-of-fit (K-S test)
    validRes = R(disk);
    stdR(iFeature) = std(validRes);
    isValid = isValid & ~kstest(validRes ./ stdR(iFeature), [], alpha);
    
    % TEST 3: test sigmaX
    testStat = params(4) / stdParams(4);
    pValue = 1-tcdf(testStat, numDegFreedom);
    isValid = isValid & pValue < alpha;
    
    success(iFeature) = isValid;
    
    P(iFeature,1) = P(iFeature,1) + params(1);
    P(iFeature,2) = P(iFeature,2) + params(2);
    P(iFeature,3) = params(3);
    P(iFeature,4) = params(4);
    P(iFeature,5) = params(5);
    P(iFeature,6) = params(6);
    P(iFeature,7) = params(7);
    
    stdP(iFeature,1) = stdParams(1);
    stdP(iFeature,2) = stdParams(2);
    stdP(iFeature,3) = stdParams(3);
    stdP(iFeature,4) = stdParams(4);
    stdP(iFeature,6) = stdParams(5);
    stdP(iFeature,7) = stdParams(6);
end

P = P(success,:);
stdP = stdP(success,:);

% Remove any detection which has been localised at the same position
ind = KDTreeBallQuery(P(:,1:2), P(:,1:2), repmat(minDist, size(P,1), 1));
ind = cellfun(@(c) c(2:end), ind, 'UniformOutput', false);
isInCluster = cellfun(@(c) ~isempty(c), ind);
indInCluster = find(isInCluster);

isValid = true(size(P,1),1);

for iiFeature = 1:numel(indInCluster)
    iFeature = indInCluster(iiFeature);
    
    if stdR(iFeature) < max(stdR(ind{iFeature}))
        isValid(iFeature) = false;
    end
end

P = P(isValid,:);
stdP = stdP(isValid,:);

featuresInfo.xCoord = [P(:,1), stdP(:,1)];
featuresInfo.yCoord = [P(:,2), stdP(:,2)];
featuresInfo.amp = [P(:,3), stdP(:,3)];
featuresInfo.stdAlong = [P(:,4), stdP(:,4)];
featuresInfo.stdAside = [P(:,5), stdP(:,5)];
featuresInfo.theta = [P(:,6), stdP(:,6)];
featuresInfo.bkg = [P(:,7), stdP(:,7)];
