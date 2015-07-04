function [ vmfmult ] = vonMisesFischer2dMult( theta, data, consolidateBessel)
%vonMisesFischer2dMult Compound distribution consisting of the sum of
%multiple von Mises distributions with different mu (mean,median,mode) and
%a constant offset

if(nargin < 3)
    consolidateBessel = true;
end

% consolidate the bessel function use here
vmf = vonMisesFischer2d(theta,false,~consolidateBessel);
% vmf = vonMisesFischer2d(theta,false,true);
% vmfres = @(x) data - x(2:2:end-2)*vmf(x(1:2:end-2)',x(end-1)) - x(end);
if(nargin < 2)
    data = 0;
end

    vmfmult = @vonMisesFischer2dMult_impl;

    function [F,J] = vonMisesFischer2dMult_impl(x)
        % mode of von Mises distribution
        mus = x(1:2:end-2)';
        % constant scale factor
        scaleFactor = x(2:2:end-2);
        % width
        kappa = x(end-1);
        % theta invariant offset
        offset = x(end);
        % normalization factor such that int(vmf) == 1
        if(consolidateBessel)
            bessel0_kappa = besseli(0,kappa);
        end
        if(nargout > 1)
            [vmfMat, vmfJacobian] = vmf(mus,kappa);
            if(consolidateBessel)
                vmfMat = vmfMat./bessel0_kappa;
                vmfJacobian = vmfJacobian ./bessel0_kappa;
            end
            F = data - scaleFactor*vmfMat - offset;
            J = zeros(length(theta),length(x));
            J(:,1:2:end-2) = bsxfun(@times,-scaleFactor,vmfJacobian(:,:,1)');
            J(:,2:2:end-2) = -vmfMat';
            if(consolidateBessel)
                bessel1_kappa = besseli(1,kappa);
                J(:,end-1) = -scaleFactor*(vmfJacobian(:,:,2) - bessel1_kappa/bessel0_kappa * vmfMat);
            else
                J(:,end-1) = -scaleFactor*vmfJacobian(:,:,2);
            end
            J(:,end) = -1;
        else
            vmfMat = vmf(mus,kappa);
            if(consolidateBessel)
                vmfMat = vmfMat./bessel0_kappa;
            end
            F = data - scaleFactor*vmfMat - offset;
        end
    end
end

