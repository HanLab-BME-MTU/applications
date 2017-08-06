function [ xg, yg, dtn_dnm ] = quickSearch( xg, yg, R, r, c)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

import orientationSpace.diffusion.*;
D = 2*pi^2;

rho_hat = fft(shiftdim(real(R.a(r,c,:)),2));

%             n_new = 2*R.filter.sampleFactor*K_new+1;

            % The convolution of two Gaussians results in a Gaussian
            % The multiplication of two Gaussians results in a Gaussian
            % The signal has been convolved with a Gaussian with sigma = pi/obj.n
            % Compute the signal convoled with a Gaussian with sigma = pi/n_new
            % Note if n == n_new, we divide by 0. Then s_inv = Inf
%             s_inv = sqrt(obj.n^2.*n_new.^2./(obj.n.^2-n_new.^2));
%             s_hat = s_inv/(2*pi);
%             x = -ceil(R.filter.K)*R.filter.sampleFactor:ceil(R.filter.K)*R.filter.sampleFactor;
            x = [0:ceil(R.filter.K)*R.filter.sampleFactor -ceil(R.filter.K)*R.filter.sampleFactor:-1];
            
%             tic
            xi = shiftdim(x*1i,1);
            nDerivs = 2*1+1;
            % second derivative
%             rho_hat = bsxfun(@times,rho_hat,xi.*xi);
%             for n=3:nDerivs
%                 rho_hat(:,:,n-1) = bsxfun(@times,rho_hat(:,:,n-2),xi);
%             end
%             toc
%             rho_hat = fft(shiftdim(real(R.a(r,c,:)),2));
%             tic
%                 xi = bsxfun(@power,xi,shiftdim(2:nDerivs,-1));
%                 rho_hat = bsxfun(@times,rho_hat,xi);
%             toc
%             rho_hat = fft(shiftdim(real(R.a(r,c,:)),2));
%             xi = shiftdim(x*1i,1);
%             tic
                xip = xi.*xi;
                for n=3:nDerivs
                    xip(:,:,n-1) = xip(:,:,n-2).*xi;
                end
                rho_hat = bsxfun(@times,rho_hat,xip);
%             toc
            
            
            % Each column represents a Gaussian with sigma set to s_hat(column)
%             f_hat = exp(-0.5 * bsxfun(@rdivide,x(:),s_hat).^2); % * obj.n/n_new;
%             f_hat = ifftshift(f_hat,1);
            
            % Angular response will be in a column
%             a_hat = fft(squeeze(real(obj.a(r,c,:))));
            % Each column represents an angular response with order K_new(column)
%             a_hat = bsxfun(@times,a_hat,f_hat);
%             response = ifft(a_hat);

% xg = out(3,43); yg = K(43);
% tic;
for i=1:5;
%     tic
    n_new = 2*R.filter.sampleFactor*yg+1;
    s_inv = sqrt(R.n^2.*n_new.^2./(R.n.^2-n_new.^2));
    s_hat = s_inv/(2*pi);
    f_hat = exp(-0.5 * bsxfun(@rdivide,x(:),s_hat).^2);
    a_hat = bsxfun(@times,rho_hat,f_hat);
    
%     a_hat = R.getResponseAtOrderFTatPoint(r,c,yg);
    
    dtn_dnm = orientationMaximaTimeDerivatives(a_hat,yg,1,xg,2*pi,true);
    xg = xg + dtn_dnm*D;
    yg = newtonK(xg,yg,R,r,c);
%     toc
end;
% toc
% dtn_dnm = orientationMaximaTimeDerivatives(R.getResponseAtOrderFTatPoint(r,c,yg),yg,1,xg)
% dtn_dnm = orientationMaximaTimeDerivatives(rho_hat,yg,1,xg,2*pi,true)


end

