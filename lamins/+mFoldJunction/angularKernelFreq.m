function [ K ] = angularKernelFreq( M, sz_or_coords, S )
%mFoldJunction.angularKernelFreq Implement M-fold junction filter kernel from Puposki
%2016
%
% Zsuzsanna Püspöki ;; Virginie Uhlmann ; Cédric Vonesch ; Michael Unser.
% Design of Steerable Wavelets to Detect Multifold Junctions.
% IEEE Transactions on Image Processing  (Volume:25 ,  Issue: 2 )
% Page(s):
% 643 - 657
% ISSN :
% 1057-7149
% INSPEC Accession Number:
% 15683740
% DOI:
% 10.1109/TIP.2015.2507981
% Date of Publication :
% 11 December 2015
% Date of Current Version :
% 29 December 2015
% Issue Date :
% Feb. 2016
% Sponsored by :
% IEEE Signal Processing Society
% Publisher:
% IEEE


if(isstruct(sz_or_coords))
    coords = sz_or_coords;
%     sz = size(coords.theta);
else
    coords = orientationSpace.getFrequencySpaceCoordinates(sz_or_coords);
end
if(nargin < 3)
    S = (0:9)*M;
end
phi = -pi:pi/100:pi-pi/100;
w_hat = fft(wraparoundN(phi,[-pi/M pi/M]).^2);
N = bsxfun(@minus,S,S');
N(N < 0) = N(N < 0)+length(phi);
W = w_hat(N+1);
[V,D] = eig(W);
[~,idx] = min(abs(diag(D)));
K = reshape(exp(1j*bsxfun(@times,coords.theta(:),S))*V(:,idx),size(coords.theta));


end

