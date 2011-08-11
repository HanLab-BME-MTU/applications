function [constrForceFieldCorrected,forceFieldCorrected]=correctSysErrTmp(constrForceField,forceField,displField)
% this function only works for square lattices. For other lattices, to
% obtain the correctionPerForceNode one has to multiply the
% errorForce_perpix with 1/3*(support of the base function). The latter is
% the volume of the uni-base function.

if nargin<1 || isempty(constrForceField)
    load('constrForceFieldCorrected.mat');
    % copyfile('cellCellForces.mat','cellCellForces_tmp.mat')
end

if nargin<1 || isempty(forceField)
    load('forceFieldCorrected.mat');
    % copyfile('forceField.mat','forceField_tmp.mat')
end

if nargin<1 || isempty(displField)
    % load('displField.mat');
    % copyfile('displField.mat','displField_tmp.mat')
end

% update the cluster analysis:
[constrForceFieldCorrected]=TFM_part_6_clusterAnalysis(constrForceFieldCorrected,forceFieldCorrected,[],[],1);