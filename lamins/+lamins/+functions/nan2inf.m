function [ in ] = nan2inf( in )
%NAN2INF Summary of this function goes here
%   Detailed explanation goes here

in(isnan(in)) = Inf;

end

