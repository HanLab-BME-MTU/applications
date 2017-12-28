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

% Do something if pool is available

while(any(~notDone))
end




end

