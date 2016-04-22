function [ b ] = binomCoefficients( n )
%BINOMCOEFFICIENTS Calculates all the binomial coefficients for n trials

% validateattributes(n,{'numeric'},{'positive','scalar','integer'});

if(n)
    len = n + 1;
    b = zeros(1,len);
    b(end-1:end) = 1;
    for iter=1:n-1
        b = conv(b,[1  1],'same');
    end
else
    b = 1;
end

end

