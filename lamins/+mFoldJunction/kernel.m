function [ K ] = kernel( M, sz, S )
%mFoldJunction.kernel Implement M-fold junction filter kernel from Puposki
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
% 
if(nargin < 3)
    S = (0:9)*M;
end

coords = orientationSpace.getFrequencySpaceCoordinates(sz);
% phi = -pi:pi/100:pi-pi/100;
% w_hat = fft(wraparoundN(phi,[-pi/M pi/M]).^2);
% N = bsxfun(@minus,S,S');
% N(N < 0) = N(N < 0)+length(phi);
% W = w_hat(N+1);
% [V,D] = eig(W);
% [~,idx] = min(abs(diag(D)));
% radialComp = orientationSpace.radialKernel(0.05,[],sz);
% K = ifft2(reshape(exp(1j*bsxfun(@times,coords.theta(:),S))*V(:,idx),size(coords.theta)).*Fr);

coords.f = coords.f*4;
area = coords.f <= 0.5 & coords.f > 0.125;
radialComp = zeros(size(coords.f));
radialComp(area) = cos(log2(coords.f(area)*4)*pi/2);
angularComp = mFoldJunction.angularKernelFreq(M,coords,S);
K = ifft2(radialComp.*angularComp);


end

