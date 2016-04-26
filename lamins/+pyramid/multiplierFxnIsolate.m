function [ M_isolated ] = multiplierFxnIsolate( omega, N )
%MULTIPLIERFXNISOLATE Produces an isolated bump
    if(nargin < 1)
        omega = 2;
    end
    if(nargin < 2)
        N = 9;
    end
    
    M = pyramid.multiplierFxn(omega, N);
    h = pyramid.meyerH;
    M_isolated = @(p,n) M(p,n).*h(p.*2^(omega*n/N));

end

