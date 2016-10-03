function [ G ] = meyerG( n )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    syms xx
    GG = matlabFunction(xx*hypergeom([1/2, -n], 3/2, xx^2));
    G = @(x) (GG(x) - GG(-1))/(GG(1) - GG(-1))*pi/2.*(x >= -1).*(x <= 1) + (x > 1)*pi/2;

end

