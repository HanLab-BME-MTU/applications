function [ xg, Kg ] = refineBifurcation( x, K, response, Korg )
%refineBifurcation Summary of this function goes here
%   Detailed explanation goes here

    import orientationSpace.diffusion.*;

    jump_threshold = repmat(0.5,size(x));
    
    notDone = K > 1;
    Kold = K;
    xold = x;
    
    response_hat = fft(response);
    

out = interpft_extrema(getResponseAtOrderFT(response_hat,8,8:-0.01:1),1,true,1e-12,false); out = orientationSpace.diffusion.alignExtrema(out);
figure; plot(8:-0.01:1,out.');
% title(sprintf('Local maxima trace for r=%d, c=%d, m=%d',coords.r(n),coords.c(n),coords.m(n)));
xlabel('K');
ylabel('2\theta (Orientation, radians)');
out(:,1:3);
hold on
grid on
plot(K,x,'ko');

last = Inf(size(x));

while(any(notDone))
    fprintf('%02.1f%% Complete\n',100-sum(notDone)./numel(notDone)*100);
    response_at_K = getResponseAtOrderFT(response_hat,Korg,Kold);
    [~,dtn_dnm,dnK_dmn] = orientationMaximaTimeDerivatives(response_at_K,Kold,3,xold,2*pi,true);
    xg_newton = dnK_dmn(:,:,1)./dnK_dmn(:,:,2);
    xg_halley = 2*dnK_dmn(:,:,1).*dnK_dmn(:,:,2)./(2*dnK_dmn(:,:,2).^2-dnK_dmn(:,:,1).*dnK_dmn(:,:,3));
    s = abs(xg_newton) < abs(xg_halley);
    xg_delta = xg_halley;
    xg_delta(s) = xg_newton(s);
    xg = xold - xg_delta;
    plot([Kold Kold],[xold xg],'r--');
%      2*v.*vd./(2*vd.^2-v.*vdd);

    Kgpd = (xg-xold).*dnK_dmn(:,:,1)+(xg-xold).^2.*dnK_dmn(:,:,2)/2+(xg-xold).^3.*dnK_dmn(:,:,3)/6;
    Kg = NaN(size(xg));
    s = Kgpd < 0 & abs(Kgpd) < 1;
        test = linspace(xold,xg,100);
        plot(Kold+(test-xold)*dnK_dmn(:,:,1)+(test-xold).^2*dnK_dmn(:,:,2)/2+(test-xold).^3*dnK_dmn(:,:,3)/6,test,'g:');

    Kg( s) = halleyK(xg( s),Kold( s)+Kgpd(s),response(:, s),Korg);
    Kg(~s) = halleyK(xg(~s),Kold(~s)        ,response(:,~s),Korg);
%     if(Kgpd < 0 && abs(Kgpd) < 1)
%         test = linspace(xold,xg,100);
%         plot(Kold+(test-xold)*dnK_dmn(:,:,1)+(test-xold).^2*dnK_dmn(:,:,2)/2+(test-xold).^3*dnK_dmn(:,:,3)/6,test,'g:');
%         Kg = halleyK(xg,Kold+Kgpd,R,coords.r(n),coords.c(n));
%         plot([Kold Kold+Kgpd],[xg xg],'y:');
%         plot([Kold+Kgpd Kg],[xg xg],'m:');
%     end
%     if(isnan(Kg))
%         Kg = halleyK(xg,Kold,R,coords.r(n),coords.c(n));
%         plot([Kold Kg],[xg xg],'c:');
%     end

    K_change = Kg - Kold;
    Kg(abs(K_change) > jump_threshold | K_change > 1e-10) = NaN;

%     if(abs(Kg - Kold) > jump_threshold)
%         Kg = NaN;
%     end
%     if(Kg - Kold > 1e-10)
%         Kg = NaN;
%     end

    s = isnan(Kg);
    jump_threshold(s) = jump_threshold(s)/2;
    Kg(s) = Kold(s);
    xg(s) = xg(s);

%     if(isnan(Kg))
%         jump_threshold = jump_threshold/2;
%         fprintf('jump threshold now %d\n',jump_threshold);
%         Kg = Kold;
%         xg = xold;
%     end
%     plot(Kg,xg,'.');

    s = ~s;
    notDone(s) = abs(dnK_dmn(:,s,1)) > 1e-10 & abs(dnK_dmn(:,s,1)) < last(s) & Kg(s) > 1;
    Kold(s) = Kg(s);
    xold(s) = xg(s);
    last(s) = abs(dnK_dmn(:,s,1));
    dnK_dmn(:,:,:)

%     if(abs(xgd(:,:,2-1)) > abs(last))
%         break;
%     else
%         last = xgd(:,:,2-1);
%     end
end


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
    s_inv = sqrt(n_org^2.*n_new.^2./(n_org.^2-n_new.^2));
    s_hat = s_inv/(2*pi);
    f_hat = exp(-0.5 * bsxfun(@rdivide,x(:),s_hat).^2);
    responseAtOrder = bsxfun(@times,response_hat,f_hat);
end

