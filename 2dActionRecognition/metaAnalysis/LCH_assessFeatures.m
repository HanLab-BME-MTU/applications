
%% 
% Input: features + experiment information
% Output: visualizations :-)
% Dependency: TBD
% Notes: developed from LCH_clsFeaturesVisualization_Train, excluding the
% annotations
function [] = LCH_assessFeatures(featsDname,featsStrID,tsneParams,always)

close all; clc;

baseDname = [featsDname filesep featsStrID];
assert(logical(exist(baseDname,'dir')));

inFname = [baseDname filesep featsStrID '.mat'];

assert(logical(exist(inFname,'file')));

%% load
load(inFname);

cellTypeStr = 'CellTypes';
datesStr = 'Dates';
sourcesStr = 'Sources';
metEffStr = 'MetEff';
annotationStr = 'Annotations';

if(size(featsNormAll,1) == 1)
    return;
end

%% tSNE
allTSNE = tsne(featsNormAll', 'NumPCAComponents', min(tsneParams.init_dims,size(featsNormAll,1)), 'Perplexity',tsneParams.perplexity);

%% Annotations


%% Source: Tumor vs. Cell Lines vs. Melanocytes
% [uniqueSources,indsSources] = selectByUniqueID(sources);
% LCH_visualize2D(allScore,'','',sources,sourcesStr,baseDname,'PCA');
% LCH_visualize2D(allScore(:,3:end),'','',sources,sourcesStr,baseDname,'PCA_34');
% LCH_visualize2D(allScore(:,4:end),'','',sources,sourcesStr,baseDname,'PCA_45');
% LCH_visualize2D(allTSNE,'','',sources,sourcesStr,baseDname,'TSNE');

%% Inds 
% indsTumors1 = indsSources{strcmp('Tumors',uniqueSources)};
% indsCellLines1 = indsSources{strcmp('CellLines',uniqueSources)};
% indsMelanocytes1 = indsSources{strcmp('Melanocytes',uniqueSources)};



%% Cell Type and days: compare cell types in different days... this is the critical test. directory per cell type.
[uniqueCellType,indsCellTYpes] = selectByUniqueID(cellTypes);
% LCH_visualize2D(allScore,uniqueCellType,indsCellTYpes,dates,cellTypeStr,baseDname,'all_PCA',always);
LCH_visualize2D(allTSNE,uniqueCellType,indsCellTYpes,dates,datesStr,baseDname,'all_TSNE',always);
[confMatDates,confMatNormDates,normAccuracyDates] = LCH_discrimination(featsNormAll',dates);
h = figure; 
imagesc(confMatNormDates); 
caxis([0,1]); 
colorbar; 
hold on; 
title(sprintf('norm accuracy = %.2f',normAccuracyDates)); 
hold off; 
saveas(h,[baseDname filesep datesStr filesep 'confMatNorm.jpg']);

%% Day: each cell type within every day (directory per day)
[uniqueDates,indsDates] = selectByUniqueID(dates);
% LCH_visualize2D(allScore,uniqueDates,indsDates,cellTypes,datesStr,baseDname,'PCA',always);
LCH_visualize2D(allTSNE,uniqueDates,indsDates,cellTypes,cellTypeStr,baseDname,'TSNE',always);
[confMatCellTypes,confMatNormCellTypes,normAccuracyCellTypes] = LCH_discrimination(featsNormAll',cellTypes);
h = figure; 
imagesc(confMatNormCellTypes); 
caxis([0,1]); 
colorbar; 
hold on; 
title(sprintf('norm accuracy = %.2f',normAccuracyCellTypes)); 
hold off; 
saveas(h,[baseDname filesep cellTypeStr filesep 'confMatNorm.jpg']);

% LCH_visualize2D(allScore,'','',dates,datesStr,baseDname,'all_PCA',always);
LCH_visualize2D(allTSNE,'','',dates,cellTypeStr,baseDname,'all_TSNE',always);

% LCH_visualize2D(allScore(indsTumors,:),'','',dates(indsTumors),datesStr,baseDname,'Tumors_PCA',always);
% LCH_visualize2D(allTSNE(indsTumors,:),'','',dates(indsTumors),datesStr,baseDname,'Tumors_TSNE',always);

%% High vs low
metEffLabels = {'Low','High'};
metEffsCell = cell(1,length(metEffs));%mat2cell(metEffs,1,ones(1,size(metEffs,2)));
for i = 1 :length(metEffs)
    tmp = metEffs(i);
    if ~isnan(tmp)        
        metEffsCell{i} = metEffLabels{tmp+1};
    end    
end
indsNoNan = find(~isnan(metEffs));

% % [uniqueMetEff,indsMetEff] = selectByUniqueID(metEffsCell);
% LCH_visualize2D(allScore(indsNoNan,:),'','',metEffsCell(indsNoNan),metEffStr,baseDname,'all_PCA_12',always);
% LCH_visualize2D(allScore(indsNoNan,3:end),'','',metEffsCell(indsNoNan),metEffStr,baseDname,'all_PCA_34',always);
LCH_visualize2D(allTSNE(indsNoNan,:),'','',metEffsCell(indsNoNan),metEffStr,baseDname,'all_TSNE_12',always);

% LCH_visualize2D(tumorsScore,'','',metEffsCell(indsTumors),metEffStr,baseDname,'tumors_PCA_12',always);
% LCH_visualize2D(tumorsScore(:,3:end),'','',metEffsCell(indsTumors),metEffStr,baseDname,'tumors_PCA_23',always);
% LCH_visualize2D(tumorsScore(:,3:end),'','',metEffsCell(indsTumors),metEffStr,baseDname,'tumors_PCA_34',always);
% LCH_visualize2D(tumorsScore(:,4:end),'','',metEffsCell(indsTumors),metEffStr,baseDname,'tumors_PCA_45',always);
% LCH_visualize2D(tumorsTSNE,'','',metEffsCell(indsTumors),metEffStr,baseDname,'tumors_TSNE_12',always);

%% temporary!
% highVsLowPCAPval(tumorsScore,metEffsCell(indsTumors));

assert(size(featsNormAll,2) == length(metEffsCell));
[confMat,confMatNorm,normAccuracy] = highVsLowPCAPval(featsNormAll',metEffsCell); %% replace with LCH_discrimination??
h = figure; imagesc(confMatNorm); caxis([0,1]); colorbar; hold on; title(sprintf('norm accuracy = %.2f',normAccuracy)); hold off; saveas(h,[baseDname filesep metEffStr filesep 'confMatNorm.jpg']);

end

%%
function [uniqueIDs,inds] = selectByUniqueID(ids)
uniqueIDs = unique(ids);
nUniqueIDs = length(uniqueIDs);

inds = cell(1,nUniqueIDs);
for i = 1 : nUniqueIDs
    curID = uniqueIDs{i};
    inds{i} = find(strcmp(ids,curID));
end
end

%% replace with LCH_discrimination??
function [confMat,confMatNorm,normAccuracy] = highVsLowPCAPval(tumorsScore,labels)
indsHigh = find(strcmp(labels,'High'));
indsLow = find(strcmp(labels,'Low'));

npcs = size(tumorsScore,2);
ps = nan(1,npcs);
for i = 1 : npcs
    ps(i) = ranksum(tumorsScore(indsHigh,i),tumorsScore(indsLow,i));
end

nSelectFeats = min(10,size(tumorsScore,2));
n1 = length(indsHigh);
n2 = length(indsLow);
costVal = LCH_getCostVal(n1,n2);
labels = [ones(1,n1)*1 ones(1,n2)*2];
% labels = labels(randperm(length(labels)));
feats = [tumorsScore(indsHigh,1:nSelectFeats); tumorsScore(indsLow,1:nSelectFeats)];
% feats = [tumorsScore(indsHigh,2:3); tumorsScore(indsLow,2:3)];
curCls = fitcdiscr(feats,labels,'ClassNames',[1,2],'Cost',costVal,'CrossVal','on');
confMat = confusionmat(curCls.Y,kfoldPredict(curCls));
confMatNorm = confMat;
confMatNorm(1,:) = confMatNorm(1,:) ./ sum(confMatNorm(1,:));
confMatNorm(2,:) = confMatNorm(2,:) ./ sum(confMatNorm(2,:));
normAccuracy = sum(diag(confMatNorm)) /  sum(confMatNorm(:));
% quadisc = fitcdiscr(feats,labels,'DiscrimType','quadratic','ClassNames',[1,2],'Cost',costVal);
% qerror = resubLoss(quadisc)
% cvmodel = crossval(quadisc,'kfold',5);
% cverror = kfoldLoss(cvmodel)
end