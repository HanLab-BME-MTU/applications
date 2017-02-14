%% At the well level: EXPLAIN HERE
% Assaf Zaritsky, May. 2016
% addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition'));
% function [pval,pvalnoMelanocytes] = pcMetaVarianceCellTypeSource(featsAll,cellTypeAll,sourceAll,outDnameScale)
function [pval,pvalnoMelanocytes] = pcMetaVarianceCellTypeSource(D,cellTypeAll,sourceAll,outDnameScale)

addpath(genpath('/home2/azaritsky/code/extern'));

close all;

% Dvec = pdist(featsAll','cityblock');
% D = squareform(Dvec);

uniqueCellType = unique(cellTypeAll);
uniqueSource = unique(sourceAll);
nUniqueCellType = length(uniqueCellType);
nUniqueSource = length(sourceAll);

allIntraSimilarity = nan(1,nUniqueCellType);
allInterSimilarity = nan(1,nUniqueCellType);
allSources = cell(1,nUniqueCellType);

for iCellType = 1 : nUniqueCellType
    curCellType = uniqueCellType{iCellType};
    curInds = find(strcmp(cellTypeAll,curCellType)); % cell type indices
    assert(length(unique(sourceAll(curInds))) == 1);
    curSource = sourceAll{curInds(1)};
    sourceInds = find(strcmp(sourceAll,curSource)); % source indices
    sourceInds = sourceInds(~ismember(sourceInds,curInds));
    
    allIntraSimilarity(iCellType) = mean(pcMetaGetIntraSimilarities(D,curInds));
    allInterSimilarity(iCellType) = mean(pcMetaGetInterSimilarities(D,curInds,sourceInds));
    allSources{iCellType} = curSource;
end

allIntraSimilarity = allIntraSimilarity(~isnan(allIntraSimilarity));
allInterSimilarity = allInterSimilarity(~isnan(allIntraSimilarity));
allSources = allSources(~isnan(allIntraSimilarity));

[pval,pvalnoMelanocytes] = plotVarianceCellTypeSource(allIntraSimilarity,allInterSimilarity,allSources,[outDnameScale filesep 'intraInterCellTypeSimilarity.eps']);

end

%%

function [pval,pvalnoMelanocytes] = plotVarianceCellTypeSource(allIntraSimilarity,allInterSimilarity,allSources,outFname)

pval = signrank(allIntraSimilarity,allInterSimilarity);

minSS = min([allIntraSimilarity,allInterSimilarity]);
maxSS = max([allIntraSimilarity,allInterSimilarity]);
xylim = [minSS,maxSS];

cmap = colormap(hsv(3));
indsCellLines = find(strcmp(allSources,'CellLines'));
indsMelanocytes = find(strcmp(allSources,'Melanocytes'));
indsTumors = find(strcmp(allSources,'Tumors'));

pvalnoMelanocytes = signrank(allIntraSimilarity([indsCellLines,indsTumors]),allInterSimilarity([indsCellLines,indsTumors]));

fontsize = 24;
h = figure;
xlabel('Intra cell type','FontSize',fontsize);
ylabel('Inter cell type','FontSize',fontsize);
hold on;
       
title(sprintf('p = %.4f',pval),'FontSize',fontsize);
plot(allIntraSimilarity(indsCellLines),allInterSimilarity(indsCellLines),'o','MarkerEdgeColor',cmap(1,:),'LineWidth',2,'MarkerSize',8);
plot(allIntraSimilarity(indsMelanocytes),allInterSimilarity(indsMelanocytes),'o','MarkerEdgeColor',cmap(2,:),'LineWidth',2,'MarkerSize',8);
plot(allIntraSimilarity(indsTumors),allInterSimilarity(indsTumors),'o','MarkerEdgeColor',cmap(3,:),'LineWidth',2,'MarkerSize',8);
% legend('Cell lines','Melanocytes','Tumors','FontSize',fontsize,'Location','East');
plot(xylim,xylim,'--k','LineWidth',2);
xlim(xylim);
ylim(xylim);
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',fontsize);
set(h,'Color','w');
axis square;
axis equal;
%         position = get(h,'position');
%         set(h,'position',[position(1:2) round(1.2*position(3:4))]);
%
export_fig(outFname);
hold off;
end