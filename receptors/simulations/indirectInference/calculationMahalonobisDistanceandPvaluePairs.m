function [s,p]=calculationMahalonobisDistanceandPvaluePairs(thetaTarget12,thetaProbe12,vTarget12,vProbe12)

% This fuction calculates the Mahalanobis distance value using the target
%and probe intermediate statistics given.
%  INPUT:   
%      paramMatrixTarget    : target intermediate statistics array
% 
%      paramMatrixProbe     : probe intermediate statistics array
%    
%
%   OUTPUT:
%       
%      s                    : calculated Mahalanobis distance    
%      
%      p                    : p value
%
%      dof                  : degress of freedon    
% Luciana de Oliveira, December 2016

%% calculate Mahalonobis distance


% s: Mahalanobis distance

s = transpose(thetaProbe12-thetaTarget12)*(inv(vProbe12+vTarget12))*(thetaProbe12-thetaTarget12);

%p Value   

% the degrees of freedom are calculated as the number of intermediate
% statistis, in our case will be: on,off number of clusters plus density number of clusters.
% we can be calculated as the size of theta. 
  dof= size(thetaTarget12,1);
  p = 1 - chi2cdf(s,dof);
end