function [radii] = ComputeEllipsoidRadii(imObjectMask, spacing)

    if ~exist( 'spacing', 'var' )
        spacing = ones(1, ndims(imObjectMask));
    end
    
    % object pixel coords
    ptObject = ind2submat( size(imObjectMask), find(imObjectMask > 0) ); 
    
    % convert to physical space
    ptObject = bsxfun(@times, ptObject, spacing);
    
    % shift origin to mean
    ptObject = bsxfun(@minus, ptObject, mean(ptObject)); 
    
    % pca
    [V, D] = eig(cov(ptObject));
    radii = 2 * sqrt(diag(D));

end