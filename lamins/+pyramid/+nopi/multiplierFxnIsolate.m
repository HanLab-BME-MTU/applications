function [ M_isolated ] = multiplierFxnIsolate( omega, N )
%MULTIPLIERFXNISOLATE Produces an isolated bump
%
% Omega, default = 2
% N is total number of multipliers
%
    if(nargin < 1)
        omega = 2;
    end
    if(nargin < 2)
        N = 9;
    end
    
    M = pyramid.nopi.multiplierFxn(omega, N);
    h = pyramid.nopi.meyerH;
    M_isolated = @(p,n) M(p,n).*h(bsxfun(@times,p,2.^(omega*shiftdim(n(:),-ndims(p))/N)));

end

