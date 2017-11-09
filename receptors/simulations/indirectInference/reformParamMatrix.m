function [paramMatrixTargetReform,paramMatrixProbeReform]=reformParamMatrix(paramMatrixTarget,paramMatrixProbe)

% This fuction reform the probe and target paramMatrix to have the same size 
% for the calculation of Mahalanobis distance and pValue and also calculate
% the values of theta and v for target and probe
%
%  INPUT:   
%      paramMatrixTarget    : target intermediate statistics array
% 
%      paramMatrixProbe     : probe intermediate statistics array
%    
%
%   OUTPUT:
%       
%     thetaTarget: mean value of the paramMatrix target intermediate statistics
%     thetaProbe: mean value of the paramMatrix probe intermediate statistics
%     vTarget: variance covariance target paramMatrix
%     vProbe: variance covariance probe paramMatrix
%
% Luciana de Oliveira, december 2016
% 


% In this version when target and probe param matrix have different sizes 
%it takes off the extra rows of on/off rate and fill zeros in the density



%%
%calculates the maximum cluster size
 maxClusterSizeTarget=size(paramMatrixTarget,1)/3;
 maxClusterSizeProbe=size(paramMatrixProbe,1)/3;
 %call on/off rates and density matrix

%target
paramMatrixTargetOnRate = paramMatrixTarget(1:maxClusterSizeTarget,:);
paramMatrixTargetOffRate = paramMatrixTarget(maxClusterSizeTarget+1:2*maxClusterSizeTarget,:);
paramMatrixTargetDensity = paramMatrixTarget(2*maxClusterSizeTarget+1:end,:);

%probe
paramMatrixProbeOnRate = paramMatrixProbe(1:maxClusterSizeProbe,:); 
paramMatrixProbeOffRate = paramMatrixProbe(maxClusterSizeProbe+1:2*maxClusterSizeProbe,:); 
paramMatrixProbeDensity = paramMatrixProbe(2*maxClusterSizeProbe+1:end,:);    


%If there is differences in target and probe sizes, ajust parameter matrix
%% association rate

 if maxClusterSizeTarget ~= maxClusterSizeProbe
    
     %calculates the difference between the sizes
     clusterDiffOnRate=abs(maxClusterSizeTarget-maxClusterSizeProbe);

     % modify the association matrix to match probe and target sizes 
    if (maxClusterSizeTarget > maxClusterSizeProbe);
    paramMatrixTargetOnRate = paramMatrixTargetOnRate(1:maxClusterSizeTarget-clusterDiffOnRate,:);
    elseif (maxClusterSizeTarget < maxClusterSizeProbe)
    paramMatrixProbeOnRate = paramMatrixProbeOnRate(1:maxClusterSizeProbe-clusterDiffOnRate,:);    
    end
 end
       
%% dissociation rate

 if maxClusterSizeTarget ~= maxClusterSizeProbe
     
     %calculates the difference between the sizes
     
     clusterDiffOffRate=abs(maxClusterSizeTarget-maxClusterSizeProbe);
     
    if (maxClusterSizeTarget > maxClusterSizeProbe);
    paramMatrixTargetOffRate = paramMatrixTargetOffRate(1:maxClusterSizeTarget-clusterDiffOffRate,:);
    elseif (maxClusterSizeTarget < maxClusterSizeProbe)
    paramMatrixProbeOffRate = paramMatrixProbeOffRate(1:maxClusterSizeProbe-clusterDiffOffRate,:);    
    end
end

%% density

 if maxClusterSizeTarget ~= maxClusterSizeProbe
     
      %calculates the difference between the sizes
     
     clusterDiffDensity=abs(maxClusterSizeTarget-maxClusterSizeProbe);
    
    if (maxClusterSizeTarget > maxClusterSizeProbe);
    paramMatrixProbeDensity = [paramMatrixProbeDensity;zeros( clusterDiffDensity,size(paramMatrixProbeDensity,2))];
    elseif (maxClusterSizeTarget < maxClusterSizeProbe)
    paramMatrixTargetDensity=[paramMatrixTargetDensity; zeros(clusterDiffDensity,size(paramMatrixTargetDensity,2))];
    end
 end

%% removing rowns with only NaNs
 
%cutOffNaN is the maximum number off NaNs permited in a row to this row enter
%in the calculation of Mahalanobis distance

cutOffNaN=6;

 %OnRate
 
%target:
sumNaN = sum(~isnan(paramMatrixTargetOnRate),2);  
iRowNaN=find (sumNaN>cutOffNaN);
paramMatrixTargetOnRate=paramMatrixTargetOnRate(iRowNaN,:);  
paramMatrixProbeOnRate=paramMatrixProbeOnRate(iRowNaN,:);

%probe
sumNaN = sum(~isnan(paramMatrixProbeOnRate),2);  
iRowNaN=find (sumNaN>cutOffNaN);
paramMatrixProbeOnRate=paramMatrixProbeOnRate(iRowNaN,:);
paramMatrixTargetOnRate=paramMatrixTargetOnRate(iRowNaN,:);  


%OffRate

%target:
sumNaN = sum(~isnan(paramMatrixTargetOffRate),2);  
iRowNaN=find (sumNaN>cutOffNaN);
paramMatrixTargetOffRate=paramMatrixTargetOffRate(iRowNaN,:);  
paramMatrixProbeOffRate=paramMatrixProbeOffRate(iRowNaN,:);

%probe
sumNaN = sum(~isnan(paramMatrixProbeOffRate),2);  
iRowNaN=find (sumNaN>cutOffNaN); 
paramMatrixProbeOffRate=paramMatrixProbeOffRate(iRowNaN,:);
paramMatrixTargetOffRate=paramMatrixTargetOffRate(iRowNaN,:);  


%Density

%target:
sumNaN = sum(~isnan(paramMatrixTargetDensity),2);  
iRowNaN=find (sumNaN>cutOffNaN);
paramMatrixTargetDensity=paramMatrixTargetDensity(iRowNaN,:);  
paramMatrixProbeDensity=paramMatrixProbeDensity(iRowNaN,:);

%probe
sumNaN = sum(~isnan(paramMatrixProbeDensity),2);  
iRowNaN=find (sumNaN>cutOffNaN);
paramMatrixProbeDensity=paramMatrixProbeDensity(iRowNaN,:);
paramMatrixTargetDensity=paramMatrixTargetDensity(iRowNaN,:);  


    
    %update the paramMatrix
    paramMatrixTargetReform=[paramMatrixTargetOnRate;paramMatrixTargetOffRate;paramMatrixTargetDensity];
    paramMatrixProbeReform=[paramMatrixProbeOnRate;paramMatrixProbeOffRate;paramMatrixProbeDensity];


end