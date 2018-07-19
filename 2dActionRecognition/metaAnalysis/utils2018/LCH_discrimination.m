%%
function [confMat,confMatNorm,normAccuracy] = LCH_discrimination(feats,labels)

uniqueLabels = unique(labels);
nUniqueLabels = length(uniqueLabels);

ns = [];
for i = 1 : nUniqueLabels
    ns = [ns, sum(strcmp(labels,uniqueLabels{i}))];
end

costVal = LCH_getCostValMultLabels(ns);


sortLabels = [];
sortFeats = [];
for i = 1 : nUniqueLabels
    sortLabels = [sortLabels, ones(1,ns(i))*i];
    inds = strcmp(labels,uniqueLabels{i});
    sortFeats = [sortFeats; feats(inds,:)];
end

curCls = fitcdiscr(sortFeats,sortLabels,'ClassNames',1:nUniqueLabels,'Cost',costVal,'CrossVal','on');
confMat = confusionmat(curCls.Y,kfoldPredict(curCls));
confMatNorm = confMat;

for i = 1 : nUniqueLabels
    confMatNorm(i,:) = confMatNorm(i,:) ./ sum(confMatNorm(i,:));
end

normAccuracy = sum(diag(confMatNorm)) /  sum(confMatNorm(:));
% quadisc = fitcdiscr(feats,labels,'DiscrimType','quadratic','ClassNames',[1,2],'Cost',costVal);
% qerror = resubLoss(quadisc)
% cvmodel = crossval(quadisc,'kfold',5);
% cverror = kfoldLoss(cvmodel)
end