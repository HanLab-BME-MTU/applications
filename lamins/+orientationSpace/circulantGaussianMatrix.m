function [ matrix ] = circulantGaussianMatrix( N, s , tol)
%circulantGaussianMatrix Obtain the circulant gaussian matrix
    x = 0:N-1;
    if(nargin < 2 || isempty(s))
        % default
        s = x';
    elseif(isscalar(s))
        % sub-sample into N space
        s = (0:s-1)'*N/s;
    elseif(~iscolumn(s))
        % otherwise make sure it's a column
        s = s(:);
    end
    xx = bsxfun(@minus,x,s)';
    xx = wraparoundN(xx,-N/2,N/2);
    matrix = exp(-xx.^2/2);

    % do not suppress if tol not given or tol == 0
    if(nargin > 2 && tol)
        % suppress zeros less than tolerance
        matrix(abs(matrix) < tol) = 0;
    end
end

