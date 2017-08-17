function [ K_out ] = halleyK( x, K, R, r, c)
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
% for i=1:10
vk{1} = 1;
last = Inf;
while(abs(vk{1}) > 1e-12)
    rhorho = R.getResponseAtOrderFTatPoint(r,c,K_out);
    vt = cell(1,5);
    vk = cell(1,5);
    for j=[1 3 5];
        rhorhod_hat = bsxfun(@times,fft(rhorho),(ifftshift(-8:8)*1i).'.^j);
        n = (j-1)/2;
        vt{j} = interpft1([0 2*pi],rhorhod_hat,x,'horner_freq').*D.^n;
    end;
    dtdK = -4./(2*K_out+1).^3;
    d2tdK2 = 24./(2*K_out+1).^4;
    vk{1} = vt{1};
    vk{3} =  vt{3}.*dtdK;
    vk{5} =  vt{5}.*(dtdK).^2+vt{3}.*d2tdK2;
    if(abs(vk{1}) < abs(last))
        K_out = K_out - 2*vk{1}.*vk{3}./(2*vk{3}.^2-vk{1}.*vk{5});
        last = vk{1};
    else
        break;
    end
%     hold on; plot(K_out,x,'.')
end

if(nargout < 1)
    hold on; plot(K_out,x,'g.')
end

end

