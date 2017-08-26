function [ maxima, coords, response, xgd ] = descendHalleyK( maxima, coords, response, n )
%descendHalleyK Summary of this function goes here
%   Detailed explanation goes here

if(nargin > 3)
    maxima = maxima(:,n);
    coords.K = coords.K(n);
    coords.r = coords.r(n);
    coords.c = coords.c(n);
    coords.m = coords.m(n);
    response = response(:,n);
end

D = 2*pi.^2;
K_jump_threshold = 0.1;
xg_change_threshold = pi/6;
maxKdelta = 4;
Korg = coords.K(1);
response_hat = fft(response);

maxKdelta = repmat(maxKdelta,size(maxima));

xgd = interpft1_derivatives(response_hat,maxima,2:4,2*pi,true);
d_p2rho_pm2_K = -D.*xgd(:,:,3-1).^2./xgd(:,:,2-1) + D.*xgd(:,:,4-1);
d_p2rho_pm2_K = d_p2rho_pm2_K.*-4./(2*coords.K+1).^3;
K_jump_estimate = xgd(:,:,2-1) ./ d_p2rho_pm2_K ;

    invalid_K_jump_estimate = K_jump_estimate < 0;
    K_jump_estimate(invalid_K_jump_estimate) = maxKdelta(invalid_K_jump_estimate);
%     Kg = coords.K - min(K_jump_estimate,maxKdelta);
%     Kg = max(Kg,0);
%     responseAtKg = getResponseAtOrderFT(response,Korg,Kg);

% first = true;
K_jump_estimate_above_threshold = true(size(maxima));

while(any(K_jump_estimate_above_threshold))
    s = K_jump_estimate_above_threshold;
    [K_jump_estimate_above_threshold(s),maxima(s),coords.K(s),K_jump_estimate(s),maxKdelta(s)] = doJumpIteration(response_hat(:,s), maxima(1,s), coords.K(s), K_jump_estimate(1,s), maxKdelta(1,s), Korg, K_jump_threshold, xg_change_threshold, D);

%     invalid_K_jump_estimate = K_jump_estimate < 0;
%     K_jump_estimate(invalid_K_jump_estimate) = maxKdelta(invalid_K_jump_estimate);
%     xg_change = -Inf;
%     xg = NaN;
%     while(isnan(xg) || min(xg_change,2*pi-xg_change) > pi/6)
% %         Kg = coords.K - min(K_jump_estimate,maxKdelta);
% %         Kg = max(Kg,0);
% %         responseAtKg = getResponseAtOrderFT(response_hat,Korg,Kg);
% %         [xg, xgd] = halleyft( responseAtKg, maxima,true,1);
% %         xg(xgd(:,:,2) > 0) = NaN;
% %         xg_change = abs(xg - maxima);
% %         failed = isnan(xg) | min(xg_change,2*pi-xg_change) > xg_change_threshold;
% %         succeeded = ~failed;
% %         maxKdelta(failed) = maxKdelta(failed)/2;
% %         coords.K(succeeded) = Kg(succeeded);
% %         maxima(succeeded) = xg(succeeded);
% % 
% %         responseAtKg(:,succeeded) = getResponseAtOrderFT(response_hat(:,succeeded),Korg,coords.K(succeeded));
% % %         responseAtKg = getResponseAtOrderFT(response_hat,Korg,coords.K);
% %                 
% %         xgd(:,:,4) = interpft1_derivatives(responseAtKg,maxima,4,2*pi,true);
% %         d_p2rho_pm2_K = -D.*xgd(:,:,3).^2./xgd(:,:,2-1) + D.*xgd(:,:,4);
% %         d_p2rho_pm2_K = d_p2rho_pm2_K.*-4./(2*coords.K+1).^3;
% %         K_jump_estimate = xgd(:,:,2) ./ d_p2rho_pm2_K ;
% %         invalid_K_jump_estimate = K_jump_estimate < 0;
% %         K_jump_estimate(invalid_K_jump_estimate) = maxKdelta(invalid_K_jump_estimate);
% %         K_jump_estimate_above_threshold = K_jump_estimate(:) > K_jump_threshold;
%         
%         Kg = coords.K - min(K_jump_estimate,maxKdelta);
%         Kg = max(Kg,0);
%         responseAtKg = getResponseAtOrderFT(response,Korg,Kg);
%         if(isnan(xg) || min(xg_change,2*pi-xg_change) > xg_change_threshold)
%             maxKdelta = maxKdelta/2;
%             fprintf('maxKdelta now %d\n',maxKdelta)
%             Kg = Kold;
% %             xg = xold;
%         else
%             plot([Kold Kg],[xold xold],'b-.');
%             plot([Kg Kg],[xold xg],'r-.');
%         end
%     end


%     [ refined, refined_derivs ] = halleyft( v, guess, freq, deriv, TOL, maxIter, avoidNaN );
    
end


end

function [K_jump_estimate_above_threshold,maxima,K,K_jump_estimate,maxKdelta] = doJumpIteration(response_hat, maxima, K, K_jump_estimate, maxKdelta, Korg, K_jump_threshold, xg_change_threshold, D)
    K_jump = max(min(K_jump_estimate,maxKdelta),min(K_jump_threshold,abs(K)));
    Kg = K - K_jump;
    Kg = max(Kg,0);
%     K_jump = K - Kg;
    responseAtKg = getResponseAtOrderFT(response_hat,Korg,Kg);
    [xg, xgd] = halleyft( responseAtKg, maxima,true,1);
    xg(xgd(:,:,2) > 0) = NaN;
    xg_change = abs(xg - maxima);
    failed = isnan(xg); % | min(xg_change,2*pi-xg_change) > xg_change_threshold;
    succeeded = ~failed;
    maxKdelta(failed) = maxKdelta(failed)/2;
    K(succeeded) = Kg(succeeded);
    maxima(succeeded) = xg(succeeded);

    responseAtKg(:,succeeded) = getResponseAtOrderFT(response_hat(:,succeeded),Korg,K(succeeded));
%         responseAtKg = getResponseAtOrderFT(response_hat,Korg,coords.K);

    xgd(:,:,4) = interpft1_derivatives(responseAtKg,maxima,4,2*pi,true);
    d_p2rho_pm2_K = -D.*xgd(:,:,3).^2./xgd(:,:,2) + D.*xgd(:,:,4);
    d_p2rho_pm2_K = d_p2rho_pm2_K.*-4./(2*K+1).^3;
    K_jump_estimate(succeeded) = xgd(:,succeeded,2) ./ d_p2rho_pm2_K(succeeded) ;
    invalid_K_jump_estimate = K_jump_estimate < 0;
    K_jump_estimate(invalid_K_jump_estimate) = maxKdelta(invalid_K_jump_estimate);
%     K_jump_estimate_above_threshold = K_jump_estimate(:) > K_jump_threshold;
    K_jump_estimate_above_threshold = (K_jump > K_jump_threshold | succeeded) & K > 0;
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

 