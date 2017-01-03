% Uses control-based PCA to output PCA of all experiments
function  screenPCA = coeffSpatiotemporalPCA2016(allFeatures,cntrlPCA)
cntrlPCA.cntrlInds;
screenPCA.speed.score = allFeatures.speedFeats.features' * cntrlPCA.speed.coeff;
screenPCA.directional.score = allFeatures.directionalityFeats.features' * cntrlPCA.directional.coeff;
screenPCA.coordination.score = allFeatures.coordinationFeats.features' * cntrlPCA.coordination.coeff;
error('BUG HERE!');
end