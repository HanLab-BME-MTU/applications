function [ vmfmult ] = vonMisesFischer2dMult( theta, data)
%vonMisesFischer2dMult Compound distribution consisting of the sum of
%multiple von Mises distributions with different mu (mean,median,mode) and
%a constant offset
vmf = vonMisesFischer2d(theta);
% vmfres = @(x) data - x(2:2:end-2)*vmf(x(1:2:end-2)',x(end-1)) - x(end);
if(nargin < 2)
    data = 0;
end

    vmfmult = @vonMisesFischer2dMult_impl;

    function [F,J] = vonMisesFischer2dMult_impl(x)
        mus = x(1:2:end-2)';
        scaleFactor = x(2:2:end-2);
        kappa = x(end-1);
        offset = x(end);
        if(nargout > 1)
            [vmfMat, vmfJacobian] = vmf(mus,kappa);
            F = data - scaleFactor*vmfMat - offset;
            J = zeros(length(theta),length(x));
            J(:,1:2:end-2) = bsxfun(@times,-scaleFactor,vmfJacobian(:,:,1)');
            J(:,2:2:end-2) = -vmfMat';
            J(:,end-1) = -scaleFactor*vmfJacobian(:,:,2);
            J(:,end) = -1;
        else
            vmfMat = vmf(mus,kappa);
            F = data - scaleFactor*vmfMat - offset;
        end
    end
end

