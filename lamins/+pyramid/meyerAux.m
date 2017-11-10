function [ nu, G ] = meyerAux( n, epsilon )
%MEYERAUX returns a meyer auxillary function of order n

%     p = 1:2:n;
    
    G = pyramid.meyerGnum(n);

    nu = @nu_impl;

%     function G = G_impl(x)
%     end

    function nu = nu_impl(w)
        w_plus = (w+1)/epsilon;
        w_minus = (w-1)/epsilon;
        
        nu = cos(G(w_plus) - pi/2 + G(w_minus));
    end



end

