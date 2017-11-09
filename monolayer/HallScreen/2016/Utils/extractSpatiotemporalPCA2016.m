% Takes only control experiments and performs PCA on all measures
function [cntrlPCA] = extractSpatiotemporalPCA2016(allFeatures, strLabels)

assert(isControlStr2016(strLabels{1}));

cntrlInds = getControlInds(strLabels);

cntrlPCA.cntrlInds =  cntrlInds;
cntrlPCA.healingRate = allFeatures.healingRate(cntrlInds);
cntrlPCA.speed.feats = allFeatures.speedFeats.features(:,cntrlInds)';
cntrlPCA.directional.feats = allFeatures.directionalityFeats.features(:,cntrlInds)';
cntrlPCA.coordination.feats = allFeatures.coordinationFeats.features(:,cntrlInds)';

[cntrlPCA.speed.coeff,cntrlPCA.speed.score,cntrlPCA.speed.latent] = pca(cntrlPCA.speed.feats);
cntrlPCA.speed.accVariance = cumsum(cntrlPCA.speed.latent)./sum(cntrlPCA.speed.latent);
[cntrlPCA.directional.coeff,cntrlPCA.directional.score,cntrlPCA.directional.latent] = pca(cntrlPCA.directional.feats);
cntrlPCA.directional.accVariance = cumsum(cntrlPCA.directional.latent)./sum(cntrlPCA.directional.latent);
[cntrlPCA.coordination.coeff,cntrlPCA.coordination.score,cntrlPCA.coordination.latent] = pca(cntrlPCA.coordination.feats);
cntrlPCA.coordination.accVariance = cumsum(cntrlPCA.coordination.latent)./sum(cntrlPCA.coordination.latent);
end

%%
function cntrlInds = getControlInds(strLabels)
cntrlInds = [];
for i = 1 : length(strLabels)
    if isControlStr2016(strLabels{i})
        cntrlInds = [cntrlInds i];
    end
end
end