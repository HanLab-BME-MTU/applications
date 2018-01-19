
%% Given the features and precalculated PCA transformation, calculate the PCA transformation
function [score] = LCH_PcaPrecalc(feats,precalcPcaCoeff)
% features = features';
% [nObs mFeats] = size(features);
% curFeatsNorm = (features - repmat(precalcPCA.meanFeats,[nObs 1])) ./ repmat(precalcPCA.stdFeats,[nObs 1]); % zscore(A)

score = feats * precalcPcaCoeff;
end