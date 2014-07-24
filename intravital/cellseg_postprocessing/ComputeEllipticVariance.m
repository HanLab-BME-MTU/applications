function [ ellipticVariance ] = ComputeEllipticVariance( ptRegion, spacing )

    if exist( 'spacing', 'var' )       
       ptRegion = ptRegion * diag(spacing);
    end

    numPoints = size(ptRegion, 1);
    
    % compute mean
    ptMean = mean( ptRegion );

    % center the points at the mean
    ptRegionCentered = ptRegion - repmat(ptMean, [numPoints, 1]);
    
    % compute covariance
    CovMat = cov( ptRegionCentered );    
    
    % compute mahalanobis distance
    mahalDist = sqrt( dot(ptRegionCentered / CovMat, ptRegionCentered, 2) );
    
    meanMahalDist = mean( mahalDist );
    
    % compute normalized elliptic variance 
    ellipticVariance = mean( (mahalDist / meanMahalDist - 1).^2 );    
    ellipticVariance = 1 / (1 + ellipticVariance);
    
end