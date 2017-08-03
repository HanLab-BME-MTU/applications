function [ K_out ] = newtonK( x, K, R, r, c)
%NEWTONK Use Newton's method to find K where root exists
%
% INPUT
% x - theta coordinate
% K - initial K guess, scalar
% R - OrientationSpaceResponse object
% r - row of point of interest
% c - column of point of interest
%
% OUTPUT
% K_out - K where the root exists

D = 2*pi^2;

% x = 0:0.01:2*pi;
K = repmat(K,1,length(x));
K_out = K;
for i=1:10
    rhorho = R.getResponseAtOrderFTatPoint(r,c,K_out);
    rhorhod = cell(1,3);
    for j=1:3
        rhorhod{j} = ifft(fft(rhorho).*(ifftshift(-8:8)*1i).'.^j);
    end;
    K_out = K_out + interpft1([0 2*pi],rhorhod{1},x)./(interpft1([0 2*pi],rhorhod{3},x).*D.*4./(2*K+1).^3);
end
hold on; plot(K_out,x,'g.')

end

