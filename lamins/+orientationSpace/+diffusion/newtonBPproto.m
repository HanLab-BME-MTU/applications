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

while(abs(dnK_dmn(:,:,1)) > 1e-6)

xold = xg; Kold = Kg;
[~,dtn_dnm,dnK_dmn] = orientationMaximaTimeDerivatives(R.getResponseAtOrderFTatPoint(coords.r(n),coords.c(n),Kg),Kg,3,xg,2*pi,false);
if(abs(dnK_dmn(:,:,1)./dnK_dmn(:,:,2)) < 0.02 || abs(dnK_dmn(:,:,1)./dnK_dmn(:,:,2)) < 0.2 && abs(dnK_dmn(:,:,1)) < 10)
    xg = xg-dnK_dmn(:,:,1)./dnK_dmn(:,:,2)
    test = linspace(xold,xg,100);
    plot(Kg+(test-xold)*dnK_dmn(:,:,1)+(test-xold).^2*dnK_dmn(:,:,2)/2+(test-xold).^3*dnK_dmn(:,:,3)/6,test,'--');
    Kg = halleyK(xg,Kg+(xg-xold)*dnK_dmn(:,:,1)+(xg-xold).^2*dnK_dmn(:,:,2)/2+(xg-xold).^3*dnK_dmn(:,:,3)/6,R,coords.r(n),coords.c(n));
    plot(Kg,xg,'.');
    if(abs(last) < abs(dnK_dmn(:,:,1)))
        dnK_dmn(:,:,1)

        break;
    else
        last = dnK_dmn(:,:,1);
    end
else
    dnm_dKn(:,:,1) = 1./dnK_dmn(:,:,1);
    dnm_dKn(:,:,2) = -dnK_dmn(:,:,2).*dnm_dKn(:,:,1).^3;
    dnm_dKn(:,:,3) = (-3*dnK_dmn(:,:,2).*dnm_dKn(:,:,2)-dnK_dmn(:,:,3).*dnK_dmn(:,:,1).^(-2)).*dnK_dmn(:,:,1).^(-2);
    Kg = Kg - 0.5;
    xg = halleyft( R.getResponseAtOrderFTatPoint(coords.r(n),coords.c(n),Kg), xg - 0.5./dnK_dmn(:,:,1),false,1);
    plot(Kg,xg,'.');
    last = Inf;
end

if(Kg < 1)
    dnK_dmn(:,:,1)

    break;
end





dnK_dmn(:,:,1)

end


end

