function [ G ] = meyerGnum( n )
%MEYERGNUM Obtain the finite smoothness Meyer window G function

if(n == Inf)
    G = @pyramid.meyerGInf;
    return;
end

powers = 1:2:2*n+1;
coefficients = pyramid.binomCoefficients(n) ./ powers .* (-1).^(powers/2-0.5);
xintercept = sum(coefficients);
scale = 1;
scale = pi/2/(G_impl(1) - G_impl(-1));

G = @G_impl;
   
    function y = G_impl(x)
        y = zeros(size(x));
        y(x < -1) = 0;
        y(x > 1) = pi/2;
        xin = x >= -1 & x <= 1;
        x = x(xin);
        y(xin) = scale*(xintercept + coefficients * bsxfun(@power,x(:)',powers'));
    end

end

