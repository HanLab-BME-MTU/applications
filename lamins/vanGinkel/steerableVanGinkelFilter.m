function [ response, nms, theta, angularResponse, angularLocalMaxima  ] = steerableVanGinkelFilter( I, f_c, b_f, K, n )
%steerableVanGinkelFilter Summary of this function goes here
%   Detailed explanation goes here

I = double(I);
% I = gpuArray(I);
If = fft2(I);

if(nargin < 2 || isempty(f_c))
    f_c = 0.1;
end

if(nargin < 3 || isempty(b_f))
    b_f = f_c*0.8;
end

if(nargin < 4 || isempty(K))
    K = 8;
end

if(nargin < 5 || isempty(n))
    n = 2*K+1;
end

f = cell(1,n);
angularResponse = cell(1,n);
% nms = cell(1,n);

angles = 0:pi/n:pi-pi/n;
outAngles = angles;
wrap = outAngles > pi/2;
outAngles(wrap) = outAngles(wrap) - pi;

for t=1:n
    f{t} = fftshift(steerableVanGinkelKernel(f_c, b_f, K ,angles(t),length(If)/2));
    angularResponse{t} = apply_freq_filter(If,f{t});
%     angularResponse{t} = gather(imfilter(I,f{t},'symmetric'));
%     nms{t} = nonMaximumSuppression(response{t},theta);
end

angularResponse = cat(3,angularResponse{:});
[response,theta] = max(real(angularResponse),[],3);
[response_i,theta_i] = max(imag(angularResponse).^2,[],3);
response = response + 1j*response_i;
theta = theta + 1j*theta_i;
if(nargout > 2)
    theta = outAngles(real(theta)) + 1j*outAngles(imag(theta));
%     theta = theta*pi/n;
    nms =   nonMaximumSuppression(real(response),real(theta)) ...
          + 1j .* nonMaximumSuppression(imag(response),imag(theta));
%     theta = theta - pi/2;
end
if(nargout > 4)
    angularLocalMaxima = real(angularResponse) > real(angularResponse(:,:,[end 1:end-1])) ...
                       & real(angularResponse) > real(angularResponse(:,:,[2:end 1])) ...
                       & real(angularResponse(:,:,[2:end 1])) > real(angularResponse(:,:,[3:end 1 2])) ...
                       & real(angularResponse(:,:,[end 1:end-1])) > real(angularResponse(:,:,[end-1 end 1:end-2])) ...
                       + 1j.*(imag(angularResponse) > imag(angularResponse(:,:,[end 1:end-1])) ...
                       & imag(angularResponse) > imag(angularResponse(:,:,[2:end 1])) ...
                       & imag(angularResponse(:,:,[2:end 1])) > imag(angularResponse(:,:,[3:end 1 2])) ...
                       & imag(angularResponse(:,:,[end 1:end-1])) > imag(angularResponse(:,:,[end-1 end 1:end-2])));
end
% nms = cat(3,nms{:});


end

function R = apply_freq_filter(If,f)
    R = ifft2(If.*real(f)) + 1j*ifft2(If.*imag(f).* -1j);
end