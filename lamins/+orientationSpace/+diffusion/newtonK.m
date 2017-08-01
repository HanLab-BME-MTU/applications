function [ output_args ] = newtonK( x, K, R)
%NEWTONK Summary of this function goes here
%   Detailed explanation goes here

% x = 0:0.01:2*pi;
KK = repmat(K,1,length(x));
KKK = KK;
rhorho = R.getResponseAtOrderFTatPoint(r,c,KKK);
for i=1:3; rhorhod{i} = ifft(fft(rhorho).*(ifftshift(-8:8)*1i).'.^i); end;
KKK = KKK + interpft1([0 2*pi],rhorhod{1},x)./(interpft1([0 2*pi],rhorhod{3},x).*D.*4./(2*KK+1).^3);
hold on; plot(KKK,x,'.')

end

