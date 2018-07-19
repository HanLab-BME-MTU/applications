function [ h, H, G ] = meyerH( epsilon, n_G )
%MEYERH Meyer wavelet in frequency domain
%
% epsilon - approximation parameter that governs transition speed
% n_G - order of the meyerG smoothing polynomial

if(nargin < 1)
    epsilon = 0.1;
end
if(nargin < 2)
    n_G = 4;
end

G = pyramid.meyerGnum(n_G);
H = @(g) G((g+1)/epsilon)-pi/2+G((g-1)/epsilon);
h = @(g) cos(H(log2(2^(1+epsilon)*2*g)))/sqrt(2);

end

