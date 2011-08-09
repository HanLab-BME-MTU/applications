function [constrForceFieldCorrected,forceFieldCorrected]=correctSysErr(constrForceField,forceField,displField)
% this function only works for square lattices. For other lattices, to
% obtain the correctionPerForceNode one has to multiply the
% errorForce_perpix with 1/3*(support of the base function). The latter is
% the volume of the uni-base function.

if nargin<1 || isempty(constrForceField)
    load('cellCellForces.mat');
end

if nargin<1 || isempty(forceField)
    load('forceField.mat');
end

if nargin<1 || isempty(displField)
    load('displField.mat');
end

% first estimate the systematic error:
numFrames=length(constrForceField);
errorStress       =zeros(numFrames,2);
errorStress_perpix=zeros(numFrames,2);
area        = zeros(numFrames,1);
gridSpacing = zeros(numFrames,1);

for iframe=1:numFrames
    if ~isempty(constrForceField{iframe})
        factor_Pa2nN=constrForceField{iframe}.par.gridSpacing^2*constrForceField{iframe}.par.pixSize_mu^2*10^(-3);
        % read out the error force and convert it to stress:
        errorStress(iframe,:)            = constrForceField{iframe}.errorSumForce.vec/factor_Pa2nN;
        area(iframe)                     = sum(constrForceField{iframe}.segmRes.maskDilated(:));
        errorStress_perpix(iframe,:)     = errorStress(iframe,:)/area(iframe);        
        gridSpacing(iframe)              = constrForceField{1}.par.gridSpacing;
    end
end

% the estimate for the systematic error is:
sysErrStressEst_perpix=nanmean(errorStress_perpix,1); % this is the correction needed per pix in nN!

for iframe=1:numFrames
    crrPFN=sysErrStressEst_perpix*gridSpacing(iframe)^2;
    forceFieldCorrected(iframe).pos=forceField(iframe).pos;
    forceFieldCorrected(iframe).vec=horzcat(forceField(iframe).vec(:,1)-crrPFN(:,1),forceField(iframe).vec(:,2)-crrPFN(:,2));
    forceFieldCorrected(iframe).par=forceField(iframe).par;
    forceFieldCorrected(iframe).par=forceField(iframe).par;
end


% This is the same step as used in TFM_part_3_calcForces
[forceFieldCorrected]=shiftForceField(forceFieldCorrected,displField);


% This is the same step as used in cutOutForceFieldManyCells
[constrForceFieldCorrected]=updateConstrForceField(constrForceField,forceFieldCorrected,'pixIntp');

% save the preliminary results:
save('forceFieldCorrected.mat', 'forceFieldCorrected');
save('cellCellForcesCorrected.mat','constrForceFieldCorrected','-v7.3');


% update the elastic energy
% [constrForceFieldCorrected]=TFM_part_5_calcElEnergies(constrForceFieldCorrected,forceFieldCorrected,displField,2^11,[],1);

% update the elastic energy
[forceFieldCorrected]=TFM_part_6_clusterAnalysis(constrForceFieldCorrected,forceFieldCorrected,[],[],1);