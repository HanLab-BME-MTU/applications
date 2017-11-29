function [ K_out , vd2] = halleyK( x, K, R, r, c)
%NEWTONK Use Newton's method to find K where root exists
%
% INPUT
% x - theta coordinate
% K - initial K guess, scalar
% R - OrientationSpaceResponse object / N x X matrix
% r - row of point of interest / K_org
% c - column of point of interest
%
% OUTPUT
% K_out - K where the root exists

D = 2*pi^2;

if(isempty(x))
    K_out = zeros(size(x));
    vd2 = zeros(size(x));
    return;
end

if(isscalar(K))
    K_out = repmat(K,size(x));
else
    K_out = K;
end
x_sz = size(x);

if(nargin > 4 && isa(R,'OrientationSpaceResponse'))
    response = squeeze(real(R.a(r,c,:)));
    response_sz = size(response);
    response_hat = fft(response);
    K_org = R.filter.K;
else
    K_org = r;
    if(nargin > 4)
        freq = c;
    else
        freq = false;
    end
    response_sz = size(R);
    if(freq)
        response_hat = R;
    else
        response = real(R(:,:));
        response_hat = fft(response);
    end
end

if(isscalar(x))
    x_sz = [1 prod(response_sz(2:end))];
    x = repmat(x,[1 response_sz(2:end)]);
end

x = repmat(x,[1 1 3]);


freqM = [0:floor(response_sz(1)/2) -floor(response_sz(1)/2):-1];
freqM = shiftdim(freqM,1)*1i;
freqM = bsxfun(@power,freqM,shiftdim([1 3 5 2],-1));
% freqM = freqM.^shiftdim([1 3 5 2],-1);
freqM(:,:,2) = freqM(:,:,2)*D;
freqM(:,:,3) = freqM(:,:,3)*D^2;
response_hat = bsxfun(@times,response_hat,freqM);
response_hat_d2 = response_hat(:,:,4);
response_hat = response_hat(:,:,1:3);

last = Inf(x_sz);
TOL = 1e-12;
notDone = ':';

K_out_notDone = K_out;
last_notDone = last;
new_notDone = true;


while(any(new_notDone))
    response_hat_at_K_out = getResponseAtOrderFT(response_hat(:,notDone,:),K_org,K_out_notDone);
    vt = interpft1([0 2*pi],response_hat_at_K_out,x(1,notDone,:),'horner_freq');

    dtdK = -4./(2*K_out_notDone+1).^3;
    d2tdK2 = 24./(2*K_out_notDone+1).^4;
    
    vk = vt;
  % vk(:,:,1) = vt(:,:,1);
    vk(:,:,2) =  vt(:,:,2).*dtdK;
    vk(:,:,3) =  vt(:,:,3).*(dtdK).^2+vt(:,:,2).*d2tdK2;
    is_better = abs(vk(:,:,1)) < last_notDone;
    new_notDone = is_better & abs(vk(:,:,1)) > TOL;
    vk = vk(:,is_better,:);
    K_out_notDone(is_better) = K_out_notDone(is_better) - 2*vk(:,:,1).*vk(:,:,2)./(2*vk(:,:,2).^2-vk(:,:,1).*vk(:,:,3));
    last_notDone(is_better) = abs(vk(:,:,1));
    
    K_out(1,notDone) = K_out_notDone;
    last(1,notDone) = last_notDone;
%     new_notDone = is_better & abs(vk(:,:,1)) > TOL;
    if(ischar(notDone))
        notDone = new_notDone;
    else
        notDone(notDone) = new_notDone;
    end
    
    K_out_notDone = K_out_notDone(new_notDone);
    last_notDone  = last_notDone(new_notDone);
    
end

K_out(last > TOL) = NaN;

if(nargout < 1)
    hold on; plot(K_out,x(:,:,1),'g.')
end

if(nargout > 1)
    vd2 = interpft1([0 2*pi],getResponseAtOrderFT(response_hat_d2,K_org,K_out),x(:,:,1),'horner_freq');
    vd2 = reshape(vd2,x_sz);
end

K_out = reshape(K_out,x_sz);

end

function responseAtOrder = getResponseAtOrderFT(response_hat,Korg,Kg)
    if(isempty(Kg))
        responseAtOrder = NaN(size(response_hat));
        return;
    end
%     x = [0:ceil(R.filter.K)*R.filter.sampleFactor -ceil(R.filter.K)*R.filter.sampleFactor:-1];
%     x = [0:8 -8:-1];
    x = [0:floor(size(response_hat,1)/2) -floor(size(response_hat,1)/2):-1];
    n_org = 2*Korg+1;

    n_new = 2*Kg+1;
    s_inv = sqrt(n_org.^2.*n_new.^2./(n_org.^2-n_new.^2));
    s_hat = s_inv/(2*pi);
    f_hat = exp(-0.5 * bsxfun(@rdivide,x(:),s_hat).^2);
    responseAtOrder = bsxfun(@times,response_hat,f_hat);
end

