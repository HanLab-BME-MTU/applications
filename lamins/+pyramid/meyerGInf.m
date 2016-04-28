function [ y ] = meyerGInf( w )
%MEYERGINF Summary of this function goes here
%   Detailed explanation goes here

lambda_plus = zeros(size(w));
filter = 1+w > 0;
lambda_plus(filter) =  exp(-1./(1+w(filter)).^2);

lambda_minus = zeros(size(w));
filter = 1-w > 0;
lambda_minus(filter) = exp(-1./(1-w(filter)).^2);

y = pi/2*lambda_plus./(lambda_plus+lambda_minus);


end

