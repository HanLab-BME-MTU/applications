function [s,p]=calculationMahalonobisDistanceandPvalue(paramMatrixTarget,paramMatrixProbe)

% This fuction calculates the Mahalanobis distance value using the target
%and probe intermediate statistics given.
%
%
%  INPUT:   
%     paramMatrixTarget    : target intermediate statistics array
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
%      
% Luciana de Oliveira, July 2016


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reform of the parameters and calculation of theta and v

[thetaTarget,thetaProbe,vTarget,vProbe]=calcThetaVaramCovMatrixForMahalonobisDistance(paramMatrixTarget,paramMatrixProbe);

%% calculation of Mahalanobis distance

% s: Mahalanobis distance

s = transpose(thetaProbe-thetaTarget)*(pinv(vProbe+vTarget))*(thetaProbe-thetaTarget);

%p Value   

% the degrees of freedom are calculated as the number of intermediate
% statistis, in our case will be: on,off number of clusters plus density number of clusters.
% we can be calculated as the size of theta. 
  dof= size(thetaTarget,1);
  p = 1 - chi2cdf(s,dof);

end