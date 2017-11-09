function [ M_h, m_h , alpha ] = multiplierFxn( omega, N )
%MULTIPLIERFXN Summary of this function goes here
%   Detailed explanation goes here

    if(nargin < 1)
        omega = 2;
    end
    if(nargin < 2)
        N = 9;
    end

    alpha = sqrt(4685)/14055*[125 101*sqrt(2) 53*sqrt(2) 16*sqrt(2) 2*sqrt(2)];
    M_h = @M;
    m_h = @m;
    function v = m(p)
        v = alpha(1)/sqrt(N)+sqrt(2/N)*alpha(2:end)*cos(2*pi*(1:4)'*p(1:end)/omega);
        v = reshape(v,size(p));
    end
    function v = M(p,n)
        if(nargin < 2)
            n = 0;
        end
        n = shiftdim(n(:),-ndims(p));
        v = m(log2(bsxfun(@times,2.^(omega*n/N),p)));
        v(repmat(p == 0,[1 1 numel(n)])) = 0;
    end

end

