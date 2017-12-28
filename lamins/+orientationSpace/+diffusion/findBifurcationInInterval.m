function [ maxima_bp, K_bp ] = findBifurcationInInterval( response, response_K, maxima, K_high, K_low )
%findBifurcationInInterval Find bifurcation point in K interval
% 
% INPUT
% response - NxM matrix or OrientationResponse object
% response_K - 1xM matrix or empty to derive from response
% maxima - 1xM matrix of location of local orientation maxima at K_high
% K_high - 1xM matrix of higher K of interval
% K_low  - 1xM matrix of lower  K of interval
%
% OUTPUT
% maxima_bp - 1xM matrix of local orientation maxima at bifurcation point
% K_bp - 1xM matrix of K at bifurcation point

K_high_temp = max(K_high,K_low);
K_low = min(K_high,K_low);
K_high = K_high_temp;
clear K_high_temp;


%% Process Response
K = floor(size(v,1)/2);
freqM = ifftshift(-K:K).'*1i;

% Transform into Fourier space
response_hat = fft(response);
% Take 1st derivative
response_hat = bsxfun(@times,response_hat,freqM);

% Do something if pool is available

while(any(~notDone))
end




end

function [maxima_high, K_high, K_low] = halveBifurcationIntervalHat(response_hat, response_K, maxima, K_high, K_low)
    % Calculate the midpoint
    K_mid = (K_high+K_low)/2;
    [maxima_high, K_high, K_low] = cutBifurcationIntervalHat(response_hat, response_K, maxima, K_high, K_low, K_mid);
end
function [maxima_high, K_high, K_low] = halveBifurcationIntervalHatT(response_hat, response_K, maxima, K_high, K_low)
    % Calculate the midpoint
    t_high = 1./(2*K_high+1).^2;
    t_low = 1./(2*K_low+1).^2;
    t_mid = (t_high + t_low)/2;
    K_mid = 1./sqrt(t_mid)./2-1;
    [maxima_high, K_high, K_low] = cutBifurcationIntervalHat(response_hat, response_K, maxima, K_high, K_low, K_mid);
end
function [maxima_high, K_high, K_low] = cutBifurcationIntervalHat(response_hat, response_K, maxima, K_high, K_low, K_mid)
    response_mid = orientationSpace.getResponseAtOrderVecHat(response_hat,response_K,K_mid);
    
    % Find continuation at midpoint
    [maxima_high, maxima_high_deriv] = halleyft(response_mid,maxima);
    maxima_high(maxima_high_deriv(:,:,3) > 0) = NaN;
    
    % Decide whether to keep high or low K of old interval
    keep_high = isnan(maxima_high);
    keep_low = ~keep_high;
    
    % If keeping high, use old maxima
    maxima_high(keep_high) = maxima;
    K_high(keep_low) = K_mid(keep_low);
    K_low(keep_high) = K_mid(keep_high);
    
end