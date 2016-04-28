function [ maxima ] = multiplierFxnMaxima( n, N, octave )
%MULTIPLERFXNMAXIMA Summary of this function goes here
%   Detailed explanation goes here

if(nargin < 3)
    octave = [0 -1 -2 -3];
end

maxima = 4^(-n/N)*2*4.^octave;

end

