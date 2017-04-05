function [ thetaTargetFinal,thetaProbeFinal,vTargetFinal,vProbeFinal ] = combineStaticDynamicFor(cellArrayVprobeS,cellArrayVtargetS,Vprobe,Vtarget)
%combineStaticDynamicParamMatrix will combine paramMatrix static and
%Dynamic in a unique matrix. Following the correct format for the variance
%covariance matrix.
%
%Luciana de Oliveira, February, 2017.

%Input

%            cellArrayThetaD: cell array containing theta for target
%            and probe for the daynamic data. The dimensions follow
%            (2* length(lR)) for the probe input.
%            cellArrayThetaS: cell array containing theta for target
%            and probe for the static data. The dimensions follow
%            (2* length(lR)) for the probe input.
%                   
% 
%            cellArrayVD: cell array containing matrix variance-covariance
%            for target and probe for the daynamic data. The dimensions follow
%            (2* length(lR)) for the probe input.
%             
%            cellArrayVS: cell array containing matrix variance-covariance
%            and probe for the static data. The dimensions follow (2* length(lR))
%            for the probe input.
%
% Outputs:
%             thetaTarget,thetaProbe,vTarget,vProbe


    % determine the size of the new vMatrix
    
      for indexCellArray=1:length(cellArrayVprobeS)
         sizeNewMatrix{indexCellArray}= size(cellArrayVprobeS{indexCellArray},1);
      end
  
