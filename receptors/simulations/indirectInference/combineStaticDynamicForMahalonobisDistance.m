function [ thetaTargetFinal,thetaProbeFinal,vTargetFinal,vProbeFinal ] = combineStaticDynamicForMahalonobisDistance(cellArrayThetaTarget,cellArrayThetaProbe,...
    cellArrayVtarget,cellArrayVprobe,sizeV)
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
%            vTargetFinal: varamCov target matrix with dynamic and static data
%            together in the final format.
%            vProbeFinal: varamCov probe matrix with dynamic and static data
%            together in the final format.
%
% Luciana de Oliveira, February 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% reserve memory for theta
% thetaTargetFinal=zeros(sum(sizeVTarget),1);

% combine theta

%target

thetaTargetFinal=cellArrayThetaTarget{1};

for indexCellArray=2:length(cellArrayThetaTarget)
    thetaTargetFinal=[thetaTargetFinal;cellArrayThetaTarget{indexCellArray}];
end

%probe

thetaProbeFinal=cellArrayThetaProbe{1};

for indexCellArray=2:length(cellArrayThetaProbe)
    thetaProbeFinal=[thetaProbeFinal;cellArrayThetaProbe{indexCellArray}];
end


% reserve space for the new vMatrix
% target
vTargetFinal=zeros(sum(sizeV));

% probe
vProbeFinal=zeros(sum(sizeV));

VcumSum= cumsum(sizeV);

%Fill values for vTarget

% Fill the first value
vTargetFinal(1:length(cellArrayVtarget{1}),1:length(cellArrayVtarget{1}))=cellArrayVtarget{1};
% fill the first value for probe
vProbeFinal(1:length(cellArrayVprobe{1}),1:length(cellArrayVprobe{1}))=cellArrayVprobe{1};

for indexCellArray=1:length(cellArrayVtarget)-1
    vTargetFinal(VcumSum(indexCellArray)+1:VcumSum(indexCellArray+1),VcumSum(indexCellArray)+1:VcumSum(indexCellArray+1))=cellArrayVtarget{indexCellArray+1};
    
    vProbeFinal(VcumSum(indexCellArray)+1:VcumSum(indexCellArray+1),VcumSum(indexCellArray)+1:VcumSum(indexCellArray+1))=cellArrayVprobe{indexCellArray+1};
end


