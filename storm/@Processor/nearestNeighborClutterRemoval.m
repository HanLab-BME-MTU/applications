function nearestNeighborClutterRemoval(obj,k)
% function nearestNeighborClutterRemoval(obj,k)
% SYNOPSIS:
% Same notation is used as in the paper, ref to paper
%
% REQUIRED INPUTS:
% - k
% The number of nearest neighbors
%
% OPTIONAL INPUTS:
%
% NEEDED PROPERTIES:
% obj.data.points
% obj.data.nPoints
% obj.data.clusters
% obj.data.nullCluster
% obj.data.neighbors
% obj.data.parents
% obj.data.error
% obj.data.modelBezCP
% obj.data.modelType
% obj.data.modelLength
% obj.data.modelRes
% obj.data.modelProj
% obj.data.modelIsOutOfDate
% obj.data.modelVar
% obj.data.roiPosition
% obj.data.roiSize
% obj.data.intensity
% obj.data.frame
% obj.data.weights
%
% MODIFIED PROPERTIES:
% obj.data.points
% obj.data.clusters
% obj.data.nullCluster
% obj.data.modelIsOutOfDate
% obj.data.orientation
% obj.data.magnitude
% obj.data.roiPosition
% obj.data.roiSize
% obj.data.intensity
% obj.data.frame
% obj.data.weights
%
% OUTPUTS:
%
% Pascal Bérard, October 2011

% Determine the dimensionality of the points
d = nnz(any(obj.data.points ~= 0,1));

% Estimate the size of the area in which every point should have at least k neighbors
extent = max(obj.data.points,[],1) - min(obj.data.points,[],1);
vol = prod(extent(1:d)); % Total data hypervolume
dens = obj.data.nPoints/vol; % Point density
avVol = k/dens; % Hypervolume per point

maxAttempts = 1;
for i=1:maxAttempts
    vol = avVol * 100^i; % Use a security margin
    rad = nthroot(vol,d)/2 % Query radius
    
    % TODO: Mirror the boundary points
    
    [~,dist] = KDTreeBallQuery(obj.data.points,obj.data.points,rad);
    haveEnoughNeighbors = cellfun(@(a) numel(a)>=k+1,dist); % The point itself is also found => k+1
    
    % TODO: Only increase the radius of the points that have not enough
    % neighbors
    
    if all(haveEnoughNeighbors)
        break;
    elseif i == maxAttempts
        % % %         assert(false,'Process: %d points have not k nearest neighbors!\n',nnz(~haveEnoughNeighbors));
    end
end


% Compute the kth nearest neighbor distances
kNNdist = cellfun(@(a) a(k+1),dist(haveEnoughNeighbors)); % kth nearest neighbor distance

% Initialize the labels
% delta_i == 1 => feature
% delta_i == 0 => clutter
delta = ones(numel(kNNdist),1);

% % % % Find the first local minimum after the first local maximum
% % [count,xout] = hist(kNNdist,round(sqrt(obj.data.nPoints)));
aaaa=obj.data.nPoints;
% % stepUp = (diff(count)>=0);
% % idx = find(~[true,stepUp(1:end-1)] & stepUp,3);
% % firstLocMin = xout(idx(3))
% % % % 


mean(kNNdist)
delta(mean(kNNdist) > kNNdist) = 0;

% Initialize the poisson rates
lambda_feature = poissonRate(kNNdist,delta);
lambda_clutter = poissonRate(kNNdist,1-delta);

% Initialize the mixture coefficients
p = mixtureCoeff(delta);

% TODO: Define stopping criterion

% EM algorithm
% % % % % for i=1:10
% % % % %     i
% % % % %     % The E-step
% % % % %     fDK_feature = arrayfun(@(a) fDK(a,lambda_feature),kNNdist);
% % % % %     fDK_clutter = arrayfun(@(a) fDK(a,lambda_clutter),kNNdist);
% % % % %     delta = (p*fDK_feature)./(p*fDK_feature+(1-p)*fDK_clutter);
% % % % %
% % % % %     %%%%%aaa=(p*fDK_feature)./(p*fDK_feature)
% % % % %     mask = isnan(delta);
% % % % %     delta(mask) = 0.5;
% % % % %     % The M-step
% % % % %     lambda_feature = poissonRate(kNNdist,delta);
% % % % %     lambda_clutter = poissonRate(kNNdist,1-delta);
% % % % %     p = mixtureCoeff(delta);
% % % % % end

lambda_feature
lambda_clutter
isClutter = (delta > 0.5);
obj.data.points = obj.data.points(haveEnoughNeighbors,:);
obj.data.points = obj.data.points(~isClutter,:);
obj.data.intensity = obj.data.intensity(haveEnoughNeighbors);
obj.data.intensity = obj.data.intensity(~isClutter);
xx=linspace(min(kNNdist),max(kNNdist),200);
fDK_feature = arrayfun(@(a) fDK(a,lambda_feature),xx);
fDK_clutter = arrayfun(@(a) fDK(a,lambda_clutter),xx);
plot(xx,fDK_feature,'r','LineWidth',2)
plot(xx,fDK_clutter,'r','LineWidth',2)
% figure(1)
[n,xout] = hist(kNNdist,round(sqrt(aaaa)));
bar(xout,n/obj.data.nPoints);
hold on
X = 1:0.1:300;
plot(X,arrayfun(@(a) fDK(a,lambda_feature),X),'r','LineWidth',2);
plot(X,arrayfun(@(a) fDK(a,lambda_clutter),X),'r','LineWidth',2);
hold off

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

end



