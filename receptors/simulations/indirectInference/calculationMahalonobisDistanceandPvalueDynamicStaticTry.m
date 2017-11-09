 function [s,p]=calculationMahalonobisDistanceandPvalueDynamicStaticTry(thetaTargetFinal,thetaProbeFinal,vTargetFinal,vProbeFinal)
% function [s]=calculationMahalonobisDistanceandPvalueDynamicStatic(thetaTargetFinal,thetaProbeFinal,vTargetFinal,vProbeFinal)



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


%% calculation of Mahalanobis distance

% s: Mahalanobis distance

s = transpose(thetaProbeFinal-thetaTargetFinal)*(pinv(vProbeFinal+vTargetFinal))*(thetaProbeFinal-thetaTargetFinal);

%p Value   

% the degrees of freedom are calculated as the number of intermediate
% statistis, in our case will be: on,off number of clusters plus density number of clusters.
% we can be calculated as the size of theta. 
    dof= size(thetaTargetFinal,1);
    p = 1 - chi2cdf(s,dof);

end