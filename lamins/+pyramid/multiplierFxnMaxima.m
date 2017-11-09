function [ maxima ] = multiplierFxnMaxima( n, N, octave )
%MULTIPLERFXNMAXIMA Summary of this function goes here
%   Detailed explanation goes here

if(nargin < 3)
%     octave = [0 -1 -2 -3];
    octave = 0;
end

maxima = 4.^(-n/N-1+octave)*2*pi;

end

