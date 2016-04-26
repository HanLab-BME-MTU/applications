function [ response, theta, nms, angularResponse, angularLocalMaxima  ] = steerableOrientationSpaceFilter( I, f_c, b_f, K, nn, n )
%steerableOrientationSpaceFilter
%
% Based on the thesis by Michael van Ginkel. Chapter 3
% "Image Analysis using Orientation Sapce based on Steerable Filters".
% Delft University of Technology. October 7th, 2002.
%
% f_c: maximum frequency for the radial filter
% b_f: frequency bandwidth for the radial filter
% K: number of rotation angles through 360 degrees
% nn: number of rotation filter bank samples, default: n 
% n: number of rotation points, usually 2K+1

% Mark Kittisopikul, August 22nd, 2015
% Jaqaman Lab
% UT Southwestern

% Rather than import, explicitly state the full name of the function
%     import orientationSpace.*;

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

if(nargin < 6 || isempty(n))
     n = 2*K+1;
end

if(nargin < 5 || isempty(nn))
    nn = n;
end

% f = cell(1,n);
% angularResponse = cell(1,n);
% nms = cell(1,n);

angles = 0:pi/n:pi-pi/n;



F = orientationSpace.kernel(f_c, b_f, K ,angles,size(If));
% F = fftshift(fftshift(F,1),2);

angularResponse = apply_freq_filter(If,F);

% for t=1:n
% %     f{t} = fftshift(steerableVanGinkelKernel(f_c, b_f, K ,angles(t),length(If)/2));
%     angularResponse{t} = apply_freq_filter(If,fftshift(F(:,:,t)));
% %     angularResponse{t} = gather(imfilter(I,f{t},'symmetric'));
% %     nms{t} = nonMaximumSuppression(response{t},theta);
% end

% angularResponse = cat(3,angularResponse{:});
if(isfinite(nn))
    
    % outAngles = 0:pi/nn:pi-pi/nn;
    % wrap = outAngles > pi/2;
    % outAngles(wrap) = outAngles(wrap) - pi;
    outAnglesRidge = wraparoundN(0:  nn-1,-nn/2,nn/2) * pi/nn;
    outAnglesEdge  = wraparoundN(0:2*nn-1,-nn,  nn)   * pi/nn;
    
    if(n ~= nn)
        angularResponse = orientationSpace.upsample(angularResponse,pi/nn);
    end
    [response,theta] = max(real(angularResponse),[],3);

    % [response_i,theta_i] = max(abs(imag(angularResponse)),[],3);
    % [response_i2] = max(imag(angularResponse),[],3);
    % negRes = response_i ~= response_i2;
    % response_i(negRes) = -response_i(negRes);

    [response_i,theta_i] = max(cat(3,imag(angularResponse),-imag(angularResponse)),[],3);

    response = response + 1j*response_i;
    theta = theta + 1j*theta_i;
    if(nargout > 1)
        theta = outAnglesRidge(real(theta)) + 1j*outAnglesEdge(imag(theta));
    %     theta = theta*pi/n;
        nms =   nonMaximumSuppression(real(response),real(theta)) ...
               + 1j .* nonMaximumSuppression(imag(response),imag(theta));
    %           + 1j .* nonMaximumSuppression(abs(imag(response)),imag(theta)) .* sign(imag(response));
    %     theta = theta - pi/2;
    end
else
    [theta, response, nms] = orientationSpace.maxima(angularResponse);
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
%     R = ifft2(If.*real(f)) + 1j*ifft2(If.*imag(f).* -1j);
    R = real(ifft2(bsxfun(@times,If,real(f)))) + 1j*real(ifft2(bsxfun(@times,If.*-1j,imag(f))));
end