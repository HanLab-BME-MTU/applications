
%% Similarity and classification measures!
% Assaf Zaritsky, May. 2016
% addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition'));
function permD = pcMetaSimilarity(featsAll,D,labels,outDname,titleStr)
close all;
% Dvec = pdist(featsAll','cityblock');
% D = squareform(Dvec);

uniqueLabels = unique(labels);
uniqueLabels = uniqueLabels(~strcmp(uniqueLabels,''));

permInds = [];
nlabels = [0];
for iLabel = 1 : length(uniqueLabels)
    curLabel = uniqueLabels{iLabel};
    curInds = find(strcmp(labels,curLabel));
    permInds = [permInds curInds];
    nlabels = [nlabels, nlabels(end)+length(curInds)];
end

% cumsumLbels = cumsum(nlabels);

permD = nan(length(permInds));
permFeats = nan(size(featsAll,1),length(permInds));

for i = 1 : length(permInds)
    permFeats(:,i) = featsAll(:,permInds(i));
    for j = 1 : length(permInds)
        permD(i,j) = D(permInds(i),permInds(j));        
    end
end

%% plot 
h = figure;
imagesc(permD);
hold on;
haxes = findobj(h,'type','axes');
set(haxes,'XTick',nlabels(1:end-1)+1);
set(haxes,'YTick',nlabels(1:end-1)+1);
set(haxes,'XTickLabel',uniqueLabels);
set(get(gca,'xlabel'),'Rotation',45);
set(haxes,'YTickLabel',uniqueLabels);
set(haxes,'FontSize',8);
set(h,'Color','w');
hold off;
export_fig([outDname 'similarity_' titleStr '.eps']);

h = figure;
imagesc(permFeats);
hold on;
haxes = findobj(h,'type','axes');
set(haxes,'XTick',nlabels(1:end-1)+1);
set(haxes,'YTick',nlabels(1:end-1)+1);
set(haxes,'XTickLabel',uniqueLabels);
set(get(gca,'xlabel'),'Rotation',45);
set(haxes,'YTickLabel',uniqueLabels);
set(haxes,'FontSize',8);
set(h,'Color','w');
hold off;
export_fig([outDname 'feats_' titleStr '.eps']);

%% classification
clsLabels = [];

for i = 1 : length(uniqueLabels)
    clsLabels = [clsLabels; i*ones(nlabels(i+1)-nlabels(i),1)];    
end

[confusionMatrix,confusionMatrixNorm,nTest,nErr] = doClassify(permFeats,clsLabels);

h = figure;
imagesc(confusionMatrix);
hold on;
colormap 'jet';
colorbar;
% caxis([0,1]);
haxes = findobj(h,'type','axes');
set(haxes,'XTick',1:length(uniqueLabels));
set(haxes,'YTick',1:length(uniqueLabels));
set(haxes,'XTickLabel',uniqueLabels);
set(get(gca,'xlabel'),'Rotation',45);
set(haxes,'YTickLabel',uniqueLabels);
set(haxes,'FontSize',8);
set(h,'Color','w');
hold off;
export_fig([outDname 'cls_confMat' titleStr '.eps']);

h = figure;
imagesc(confusionMatrixNorm);
hold on;
colormap 'jet';
colorbar;
caxis([0,1]);
haxes = findobj(h,'type','axes');
set(haxes,'XTick',1:length(uniqueLabels));
set(haxes,'YTick',1:length(uniqueLabels));
set(haxes,'XTickLabel',uniqueLabels);
set(get(gca,'xlabel'),'Rotation',45);
set(haxes,'YTickLabel',uniqueLabels);
set(haxes,'FontSize',8);
set(h,'Color','w');
hold off;
export_fig([outDname 'cls_confMatNorm' titleStr '.eps']);

save([outDname 'cls_confMat' titleStr '.mat'],...
    'confusionMatrix','confusionMatrixNorm','nTest',...
    'uniqueLabels','titleStr','permFeats','clsLabels',...
    'nTest','nErr');

end



%%

function [confusionMatrix,confusionMatrixNorm,nTest,nErr] = doClassify(feats,labels)

nClass = max(labels);
confMat = zeros(nClass,nClass);

nTrain = floor(0.9*length(labels));
nTest = length(labels) - nTrain;

indsAll = 1:length(labels);

confusionMatrix = zeros(nClass);
nfail = 0;
for iter = 1 : 100
    [labelsTrain,indsTrain] = datasample(labels,nTrain,'Replace',false);
    indsTest = find(~ismember(indsAll,indsTrain));
    labelsTest = labels(indsTest);
    featsTrain = feats(:,indsTrain);
    featsTest = feats(:,indsTest);
    
    try 
        cls = fitcdiscr(featsTrain',labelsTrain);
    catch ee
        nfail = nfail + 1;
        continue;
    end
        
    [labelsPred,scoresPred] = predict(cls,featsTest');
    
    for itest = 1 : nClass        
        for ipred = 1 : nClass  
            testPred = (labelsTest == itest) & (labelsPred == ipred);
            confusionMatrix(itest,ipred) = confusionMatrix(itest,ipred) + sum(testPred);
        end
    end
end

assert (nfail < 19);
   
nTest = sum(confusionMatrix(:));
nErr = nTest-sum(diag(confusionMatrix));

confusionMatrixNorm = confusionMatrix;
for i = 1 : nClass
    confusionMatrixNorm(i,:) = confusionMatrix(i,:)./sum(confusionMatrix(i,:));
end
confusionMatrix = confusionMatrix./sum(confusionMatrix(:));
end