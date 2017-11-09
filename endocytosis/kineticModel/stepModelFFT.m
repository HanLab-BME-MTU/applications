%[pdf] = stepModelODE(k, t) calculates the lifetime distribution for a 
% multi-step model of the form:
%
%    k1     k2      kn-1
% S1 --> S2 --> ... --> Sn
%
% INPUTS
%     k : vector of rate constants
%     t : time vector
%
% OUTPUTS
%   pdf : pdf of lifetimes
%
% Note: FFT-based implementation

% Francois Aguet, 2013

function [pdf] = stepModelFFT(k, t)
N = numel(t);
dx = t(2)-t(1);
w = ((0:N-1)-floor(N/2))/dx/N*2*pi;
F = ones(size(t));
for i = 1:numel(k)
    F = F.*k(i)./(k(i)+1j*w);
end
pdf = abs(ifft(ifftshift(F)));