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
    v = cell(1,3);
    for j=[1 3];
        rhorhod_hat = bsxfun(@times,fft(rhorho),(ifftshift(-8:8)*1i).'.^j);
        v{j} = interpft1([0 2*pi],rhorhod_hat,x,'horner_freq');
    end;
    K_out = K_out + v{1}./(v{3}.*D.*4./(2*K+1).^3);
%     hold on; plot(K_out,x,'.')
end

if(nargout < 1)
    hold on; plot(K_out,x,'g.')
end

end

