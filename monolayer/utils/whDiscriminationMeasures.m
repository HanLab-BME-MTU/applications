function [DaviesBouldinIndex,DunnIndex,SilhouetteCoefficient] = whDiscriminationMeasures(featsGEF,featsControl)
meanGene = mean(featsGEF,2)';
meanControl = mean(featsControl,2)';

% distance between clusters centers
d = pdist2(meanGene,meanControl);

% intra cluster distance to center
genesDist = pdist2(featsGEF',meanGene);
controlDist = pdist2(featsControl',meanControl);

% genesStd = std(genesDist);
% controlStd = std(controlDist);
c1 = mean(genesDist); % genes mean distance to cluster center
c2 = mean(controlDist); % control mean distance to cluster center

% (c1+c2) / d
DaviesBouldinIndex = d / (c1+c2);

% max(c1,c2) / d
DunnIndex = d / max(max(genesDist),max(controlDist));

% Silhouette Coefficient
SilhouetteCoefficient = Silhouette(featsGEF,featsControl);

end

function [SilhouetteCoefficient] = Silhouette(featsGEF,featsControl)
sGEF = Silhouette1(featsGEF,featsControl); % calculate s for all GEFs
sControl = Silhouette1(featsControl,featsGEF); % calculate s for all Controls

SilhouetteCoefficient = mean([sGEF,sControl]);
end

function [s] = Silhouette1(feats1,feats2)
nFeats1 = size(feats1,2);
s = nan(1,nFeats1);

for i = 1 : nFeats1
    a = mean(pdist2([feats1(:,1:i-1), feats1(:,i+1:end)]',feats1(:,i)'));
    b = mean(pdist2(feats2',feats1(:,i)'));
    s(i) = (b - a) / max(a,b);
end
end