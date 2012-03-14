function nearestNeighborClutterRemoval(obj,k,mode)
% function nearestNeighborClutterRemoval(obj,k,varargin)
% SYNOPSIS:
% Paper: Nearest-Neighbor clutter Removal for Estimating Features in
% Spatial Point Processes
%
% REQUIRED INPUTS:
% - k
% The number of nearest neighbors
%
% - mode
% The density estimation mode: = 'Gaussian' (Gaussian mixture)
%                              = 'Poisson' (Poisson mixture)
%                              = 'Display' (Display the Gaussian and Poisson mixture)
%                              = scalar (Hard threshold)
%
% OPTIONAL INPUTS:
%
% NEEDED PROPERTIES:
% obj.data.points
% obj.data.nPoints
% obj.data.intensity
%
% MODIFIED PROPERTIES:
% obj.data.points
% obj.data.intensity
%
% OUTPUTS:
%
% Pascal Bérard, October 2011

% Determine the dimensionality of the points
assert(size(obj.data.points,2)==3,'This method has only been tested with 3D points!');
d = nnz(any(obj.data.points ~= 0,1));

% Estimate the size of the area in which every point should have at least k neighbors
extent = max(obj.data.points,[],1) - min(obj.data.points,[],1);
vol = prod(extent(1:d)); % Total data hypervolume
dens = obj.data.nPoints/vol; % Point density
avVol = k/dens; % Hypervolume per point

dist = cell(obj.data.nPoints,1);
haveEnoughNeighbors = false(obj.data.nPoints,1);

maxAttempts = 100;
for i=1:maxAttempts
    vol = avVol * 10^i; % Use a security margin
    rad = nthroot(vol,d)/2; % Query radius
    
    % TODO: Use mirror boundaries!
    [~,dist(~haveEnoughNeighbors)] = KDTreeBallQuery(obj.data.points,obj.data.points(~haveEnoughNeighbors,:),rad);
    haveEnoughNeighbors = cellfun(@(a) numel(a)>=k+1,dist); % The point itself is also found => k+1
    
    if all(haveEnoughNeighbors)
        break;
    elseif i == maxAttempts
        assert(false,'Process: %d points have not k nearest neighbors!\n',nnz(~haveEnoughNeighbors));
    end
end

% Compute the kth nearest neighbor distances
kNNdist = cellfun(@(a) a(k+1),dist(haveEnoughNeighbors)); % kth nearest neighbor distance

if strcmp(mode,'Gaussian')
    isClutter = useGaussian(kNNdist);
elseif strcmp(mode,'Poisson')
    isClutter = usePoisson(kNNdist);
elseif isscalar(mode)
    isClutter = useScalar(kNNdist,mode);
elseif strcmp(mode,'Display')
    [~,p,lambda_feature,lambda_clutter] = usePoisson(kNNdist);
    [~,gm] = useGaussian(kNNdist);
    figure(333);
    [n,xout] = hist(kNNdist,round(3*sqrt(obj.data.nPoints)));
    bar(xout,n/(sum(n)*mean(diff(xout))));
    hold on
    X = linspace(1,max(kNNdist),500)';
    plot(X,arrayfun(@(a) p*fDK(a,lambda_feature),X),'r','LineWidth',2);
    plot(X,arrayfun(@(a) (1-p)*fDK(a,lambda_clutter),X),'r','LineWidth',2);
    plot(X,gm.pdf(X),'g','LineWidth',2)
    hold off
else
    assert(0,'Process: Invalid mode specified!');
end

% Remove clutter
if ~strcmp(mode,'Display')
    nPointsOld = obj.data.nPoints;
    obj.data.points = obj.data.points(haveEnoughNeighbors,:);
    obj.data.points = obj.data.points(~isClutter,:);
    obj.data.intensity = obj.data.intensity(haveEnoughNeighbors);
    obj.data.intensity = obj.data.intensity(~isClutter);
    
    fprintf('Process: Removed %d clutter points out of %d points!\n',nPointsOld-obj.data.nPoints,nPointsOld);
end

    function lambda = poissonRate(di,delta)
        % lambda = k*sum(delta)/(pi*sum(di.^2.*delta));
        lambda = k*sum(delta)/(pi*sum(di.^2.*delta))/40;
    end

    function p = mixtureCoeff(delta)
        p = sum(delta)/numel(delta);
    end

    function out = fDK(di,lambda) % The probability density function
        term1 = d*exp(-(lambda*pi^(d/2)*di^d)/(gamma((2+d)/2)));
        term2 = ((lambda*pi^(d/2)*di^d)/gamma((2+d)/2))^k;
        out = term1*term2/(di*gamma(k));
    end

    function [isClutter,gm] = useGaussian(kNNdist)
        s = 2*ones(numel(kNNdist),1);
        s(mean(kNNdist) > kNNdist) = 1;
        gm = gmdistribution.fit(kNNdist,2,'Start',s);
        componentIdx = gm.cluster(kNNdist);
        isClutter = logical(componentIdx-1);
    end

    function isClutter = useScalar(kNNdist,threshold)
        isClutter = false(numel(kNNdist),1);
        isClutter(kNNdist > threshold) = true;
    end

    function [isClutter,p,lambda_feature,lambda_clutter] = usePoisson(kNNdist)
        % Initialize the labels
        % delta_i == 1 => feature
        % delta_i == 0 => clutter
        delta = ones(numel(kNNdist),1);
        delta(mean(kNNdist) > kNNdist) = 0;
        
        % Initialize the poisson rates
        lambda_feature = poissonRate(kNNdist,delta);
        lambda_clutter = poissonRate(kNNdist,1-delta);
        
        % Initialize the mixture coefficients
        p = mixtureCoeff(delta);
        
        % TODO: Define EM stopping criterion
        % disp('Process: WARNING: No EM stopping criterion has been defined!');
        
        % EM algorithm
        maxIter = 20;
        for iter=1:maxIter
            % The E-step
            fDK_feature = arrayfun(@(a) fDK(a,lambda_feature),kNNdist);
            fDK_clutter = arrayfun(@(a) fDK(a,lambda_clutter),kNNdist);
            delta = (p*fDK_feature)./(p*fDK_feature+(1-p)*fDK_clutter);
            mask = isnan(delta);
            delta(mask) = 0.5;
            
            % The M-step
            lambda_feature_updated = poissonRate(kNNdist,delta);
            lambda_clutter_updated = poissonRate(kNNdist,1-delta);
            
            % Test for convergence
            if (lambda_feature-lambda_feature_updated)/lambda_feature < 1e-5 && (lambda_clutter-lambda_clutter_updated)/lambda_clutter < 1e-5
                break;
            end
            
            lambda_feature = lambda_feature_updated;
            lambda_clutter = lambda_clutter_updated;
            p = mixtureCoeff(delta);
            if iter == maxIter
                disp('Process: WARNING: Maximum number of EM iteration has been reached without satisfying the stopping criterion!');
            end
        end
        
        isClutter = (delta > 0.5);
    end
end



