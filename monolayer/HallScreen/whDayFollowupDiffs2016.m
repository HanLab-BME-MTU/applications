function [] = whDayFollowupDiffs2016(followupDname,validateGenes,figDname,outputPrefix)

addpath(genpath('/home2/azaritsky/code/applications/monolayer/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all;

if ~exist(figDname,'dir')
    mkdir(figDname);
end

properties = {'Speed','Directionality','Coordination'};

nGenes = length(validateGenes);

ns = [];
N = 0;

for iProp = 1 : length(properties)
    curProp = properties{iProp};
    
    fprintf(sprintf('\n%s\n',curProp));
    
    curFeatPC1 = [];
    curFeatPC2 = [];
    curFeatPC3 = [];
    curFeatsHealingRate = [];
    for iGene = 1 : nGenes
        curGene = validateGenes{iGene};
        curFeat = [];
        
        load([followupDname filesep curProp filesep curGene filesep curProp '_stats.mat']);
        
        curN = length(controlMeanPC1);
        curn = length(controlRepPC1);
        
        if iProp == 1            
            N = N + curn;
            ns = [ns curn];
            
            diffHealingRate = -(controlRepHealingRate - geneRepHealingRate);
            curFeatsHealingRate = [curFeatsHealingRate,diffHealingRate];
        end
                
        diffPC1 = -(controlRepPC1 - geneRepPC1);
        diffPC2 = -(controlRepPC2 - geneRepPC2);
        diffPC3 = -(controlRepPC3 - geneRepPC3);
        
        
        curFeatPC1 = [curFeatPC1,diffPC1];
        curFeatPC2 = [curFeatPC2,diffPC2];
        curFeatPC3 = [curFeatPC3,diffPC3];
        
        fprintf(sprintf('N = %d, n = %d\n',curN,curn));
        if iProp == 1
            fprintf(sprintf('Healing rate: %.4f\n',signrank(diffHealingRate)));
        end
        fprintf(sprintf('PC1: %.4f\n',signrank(diffPC1)));
        fprintf(sprintf('PC2: %.4f\n',signrank(diffPC2)));
        fprintf(sprintf('PC3: %.4f\n',signrank(diffPC3)));
    end    
    
    if iProp == 1
        plotDiffs(curFeatsHealingRate,ns,validateGenes,figDname,[outputPrefix curProp '_HealingRate'],[curProp ' HealingRate']);
    end
    plotDiffs(curFeatPC1,ns,validateGenes,figDname,[outputPrefix curProp '_PC1'],[curProp ' PC1']);
    plotDiffs(curFeatPC2,ns,validateGenes,figDname,[outputPrefix curProp '_PC2'],[curProp ' PC2']);
    plotDiffs(curFeatPC3,ns,validateGenes,figDname,[outputPrefix curProp '_PC3'],[curProp ' PC3']);
    
end
end

function [] = plotDiffs(curFeats,ns,validateGenes,figDname,outputPrefix,titleStr)
nGenes = length(validateGenes);
cumsumNs = [0 cumsum(ns)];

N = floor(cumsumNs(end)+4);
if N < 100
    step = 50;
else 
    step = 100;
end

maxDiff = 1.1 * (max(abs(curFeats)));

cmap = colormap(hsv(nGenes));

% [left bottom width height]
FPosition = [0 0 900 300]; % FPosition = [0 0 1000 300];
APosition = [0.1 0.2 0.7 0.75]; 

fontsize = 10;

h = figure;
xlabel('Experiment','FontSize',fontsize);
ylabel('Control - KD','FontSize',fontsize);
hold on;
title(titleStr,'FontSize',fontsize);
for igene = 1 : nGenes
    inds = cumsumNs(igene)+1:cumsumNs(igene+1);
    curPerm = randperm(ns(igene));
    indsPerm = inds(curPerm);
    plot(inds,curFeats(indsPerm),'o','MarkerEdgeColor',cmap(igene,:),'LineWidth',2,'MarkerSize',8);
end
legend([validateGenes(:)],'Location','eastoutside');
plot([1,cumsumNs(end)],[0,0],'--k','LineWidth',2);
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[-3,(cumsumNs(end)+4)]);
set(haxes,'YLim',[-maxDiff,maxDiff]);
set(haxes,'XTick',0:step:N);
set(haxes,'XTickLabel',0:step:N);
set(haxes,'YTick',[-maxDiff,0,maxDiff]);
set(haxes,'YTickLabel',[-maxDiff,0,maxDiff]);
set(haxes,'FontSize',fontsize);
set(h,'Color','w');
set(h,'Position',FPosition,'PaperPositionMode','auto');
axisHandle= findobj(h,'type','axes');
set(axisHandle,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
set(get(axisHandle,'XLabel'),'FontSize',fontsize); set(get(axisHandle,'YLabel'),'FontSize',fontsize);

legendFname = [figDname 'legend.eps'];
if ~exist(legendFname,'file')
    export_fig(legendFname);    
end
legend off;

FPosition = [0 0 700 300];
APosition = [0.1 0.2 0.85 0.75]; 
set(h,'Position',FPosition,'PaperPositionMode','auto');
axisHandle= findobj(h,'type','axes');
set(axisHandle,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
set(get(axisHandle,'XLabel'),'FontSize',fontsize); set(get(axisHandle,'YLabel'),'FontSize',fontsize);
export_fig([figDname outputPrefix '.eps']);
hold off;
end