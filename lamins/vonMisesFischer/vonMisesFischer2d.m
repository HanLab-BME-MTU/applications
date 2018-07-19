function [ vmf ] = vonMisesFischer2d( theta , sampleNorm, useBessel)
%vonMisesFischer2d returns a function of mu and kappa

if(isscalar(theta))
    theta = -pi/2:pi/theta:pi/2-pi/theta;
end
if(nargin < 2)
    sampleNorm = false;
end
if(nargin < 3)
    useBessel = true;
end

if(sampleNorm)
    % if true, then sum(vmf(mu,kappa)) == 1
    norm = 1./length(theta);
else
    % otherwise sum(vmf(mu,kappa)) == length(theta)./(2*pi)
    norm = 1./(2*pi);
end

doubleTheta = theta*2;

% vmf = @(mu,kappa) 1./(2*pi*besseli(0,kappa)).*exp(kappa*cos(bsxfun(@minus,doubleTheta,mu*2)));
vmf = @vonMisesFischer2d_impl;

cos_theta = cos(doubleTheta);
sin_theta = sin(doubleTheta);

    function [F,J] = vonMisesFischer2d_impl(mu,kappa)
        if(useBessel)
            bessel0_kappa = besseli(0,kappa);
        else
            bessel0_kappa = 1;
        end
        
%         angle = bsxfun(@minus,doubleTheta,mu*2);
%         cos_angle = cos(angle);
        two_mu = 2*mu;
        cos_mu = cos(two_mu);
        sin_mu = sin(two_mu);
        cos_angle = zeros(length(mu),length(theta));
        for j=1:length(mu)
            cos_angle(j,:) = cos_theta*cos_mu(j) + sin_theta*sin_mu(j);
        end
        
        %F = norm*1./(bessel0_kappa).*exp(kappa*cos_angle);
        F = exp(kappa*cos_angle);
        % F = F*36/sum(F);
        F = norm.*F;
        F = F./bessel0_kappa;
        
        if(nargout > 1)
%             sin_angle = sin(angle);
            sin_angle = zeros(length(mu),length(theta));
            for j=1:length(mu)
                sin_angle(j,:) = sin_theta*cos_mu(j) - cos_theta*sin_mu(j);
            end
            if(useBessel)
                bessel1_kappa = besseli(1,kappa);
            else
                % derivative of 1
                bessel1_kappa = 0;
            end
            J = zeros([size(F) 2]);
%             J = cat(3,F,F);
            % partial derivative with respect to mu
            J(:,:,1) = 2*kappa*sin_angle.*F;
            % partial derivative with respect to kappa
            % derrivative of besseli(0,kappa) is besseli(1,kappa)
            J(:,:,2) = (-bessel1_kappa./bessel0_kappa+cos_angle).*F;
        end
    end

end

