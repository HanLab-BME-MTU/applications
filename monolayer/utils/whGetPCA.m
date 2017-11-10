%% Given the features, calculate the PCA transformation
function [out] = whGetPCA(features)
out.feats = features;
features = features';
[out.nObs out.mFeats] = size(features);
out.meanFeats = mean(features);
out.stdFeats = std(features);
out.curFeatsNorm = (features - repmat(out.meanFeats,[out.nObs 1])) ./ repmat(out.stdFeats,[out.nObs 1]); % zscore(A)

[coeff,score,out.latent] = pca(out.curFeatsNorm);

out.coeff = fixCoeff(coeff);

[out.score] = whGetPCAPrecalc(features',out);

out.accVariance = cumsum(out.latent)./sum(out.latent);
end

%%
function [outCoeff] = fixCoeff(inCoeff)
outCoeff = nan(size(inCoeff));
for i = 1 : size(inCoeff,2)
    curCoeff = inCoeff(:,i);
    if curCoeff(1) < 0 
        curCoeff = (-1) * curCoeff;
    end
    outCoeff(:,i) = curCoeff;
end
end