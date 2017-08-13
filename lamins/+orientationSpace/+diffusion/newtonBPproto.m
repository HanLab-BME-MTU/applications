function [ xg, Kg ] = newtonBPproto( R, n, coords , xg, Kg)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

import orientationSpace.diffusion.*;

out = interpft_extrema(R.getResponseAtOrderFTatPoint(coords.r(n),coords.c(n),8:-0.01:1)); out = orientationSpace.diffusion.alignExtrema(out);
figure; plot(8:-0.01:1,out.');
out(:,1:3);
hold on
grid on

dnK_dmn(:,:,1) = 1;

% xg = maxima2(n)*2; Kg = 8;

last = Inf;

D = 2*pi.^2;

maxKdelta = 4;
jump_threshold = 0.5;

while(abs(dnK_dmn(:,:,1)) > 1e-6)

xold = xg; Kold = Kg;
% [~,dtn_dnm,dnK_dmn] = orientationMaximaTimeDerivatives(R.getResponseAtOrderFTatPoint(coords.r(n),coords.c(n),Kg),Kg,3,xg,2*pi,false);
%  || abs(dnK_dmn(:,:,1)./dnK_dmn(:,:,2)) < 0.2 && abs(dnK_dmn(:,:,1)) < 10
% if(abs(dnK_dmn(:,:,1)./dnK_dmn(:,:,2)) < 0.1)
xgd = interpft1_derivatives(R.getResponseAtOrderFTatPoint(coords.r(n),coords.c(n),Kg),xg,2:4);
d_p2rho_pm2_K = -D.*xgd(:,:,3-1).^2./xgd(:,:,2-1) + D.*xgd(:,:,4-1);
d_p2rho_pm2_K = d_p2rho_pm2_K.*-4./(2*Kg+1).^3;
K_jump_estimate = xgd(:,:,2-1) ./ d_p2rho_pm2_K ;
if( K_jump_estimate < jump_threshold && K_jump_estimate > 0)
% if(true)
    [~,dtn_dnm,dnK_dmn] = orientationMaximaTimeDerivatives(R.getResponseAtOrderFTatPoint(coords.r(n),coords.c(n),Kg),Kg,3,xg,2*pi,false);
    xg_newton = dnK_dmn(:,:,1)./dnK_dmn(:,:,2);
    xg_halley = 2*dnK_dmn(:,:,1).*dnK_dmn(:,:,2)./(2*dnK_dmn(:,:,2).^2-dnK_dmn(:,:,1).*dnK_dmn(:,:,3));
    if(abs(xg_newton) < abs(xg_halley))
        xg = xg - xg_newton;
    else
        xg = xg - xg_halley;
    end
%      2*v.*vd./(2*vd.^2-v.*vdd);
    Kgpd = (xg-xold)*dnK_dmn(:,:,1)+(xg-xold).^2*dnK_dmn(:,:,2)/2+(xg-xold).^3*dnK_dmn(:,:,3)/6;
    Kg = NaN;
    if(Kgpd < 0 && abs(Kgpd) < 1)
        test = linspace(xold,xg,100);
        plot(Kg+(test-xold)*dnK_dmn(:,:,1)+(test-xold).^2*dnK_dmn(:,:,2)/2+(test-xold).^3*dnK_dmn(:,:,3)/6,test,'--');
        Kg = halleyK(xg,Kold+Kgpd,R,coords.r(n),coords.c(n));
    end
    if(isnan(Kg))
        Kg = halleyK(xg,Kold,R,coords.r(n),coords.c(n));
    end
    if(abs(Kg - Kold) > jump_threshold)
        Kg = NaN;
    end
    if(Kg > Kold)
        Kg = NaN;
    end
    if(isnan(Kg))
        jump_threshold = jump_threshold/2;
        Kg = Kold;
        xg = xold;
    end
    plot(Kg,xg,'.');
    if(xgd(:,:,2-1) < last)
        break;
    else
        last = xgd(:,:,2-1);
    end
%     if(abs(last) < abs(dnK_dmn(:,:,1)))
%         dnK_dmn(:,:,1)
% 
%         break;
%     else
%         last = dnK_dmn(:,:,1);
%     end
else
%     dnm_dKn(:,:,1) = 1./dnK_dmn(:,:,1);
%     dnm_dKn(:,:,2) = -dnK_dmn(:,:,2).*dnm_dKn(:,:,1).^3;
%     dnm_dKn(:,:,3) = (-3*dnK_dmn(:,:,2).*dnm_dKn(:,:,2)-dnK_dmn(:,:,3).*dnK_dmn(:,:,1).^(-2)).*dnK_dmn(:,:,1).^(-2);
%     Kg = Kg - 0.5;
%     xg = halleyft( R.getResponseAtOrderFTatPoint(coords.r(n),coords.c(n),Kg), xg - 0.5./dnK_dmn(:,:,1),false,1);
    if(K_jump_estimate < 0)
        K_jump_estimate = maxKdelta;
    end
    Kg = Kg - min(K_jump_estimate,maxKdelta);
    Kg = max(Kg,0);
    xg = halleyft( R.getResponseAtOrderFTatPoint(coords.r(n),coords.c(n),Kg), xg,false,1);
    xg(interpft1_derivatives(R.getResponseAtOrderFTatPoint(coords.r(n),coords.c(n),Kg),xg,2,2*pi) > 0) = NaN;
    xg_change = abs(xg - xold);
    if(isnan(xg) || min(xg_change,2*pi-xg_change) > pi/6)
        maxKdelta = maxKdelta/2;
        Kg = Kold;
        xg = xold;
    end
%     if(Kg < 1)
%         Kg = Kold;
%         xg = xold;
%         break;
%     end
    plot(Kg,xg,'.');
    last = -Inf;
end

if(Kg < 1)
    dnK_dmn(:,:,1);

    break;
end




dnK_dmn(:,:,1);

end

plot(Kg,xg,'s');



end

