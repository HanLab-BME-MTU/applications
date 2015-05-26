function [] = whControlInterdayAssessment(geneDayDiff,mainDirname,propertyStr,metaData,targetGenesStr,healingRateControl,healingRateGene)
nDayGeneSeq = length(geneDayDiff);

params.timePerFrame = metaData.timePerFrame;
params.patchSize = 15;

allMeanControls = [];
allMeanGenes = [];
for iDayGeneSeq = 1 : nDayGeneSeq    
    allMeanControls = [allMeanControls, geneDayDiff{iDayGeneSeq}.meanControl];
    allMeanGenes = [allMeanGenes, geneDayDiff{iDayGeneSeq}.meanGene];
end

metaMeanControl = mean(allMeanControls,2);
allDistsControl = pdist2(metaMeanControl',allMeanControls');
allDistsGene = pdist2(metaMeanControl',allMeanGenes');

[sortedDists, sortedInds] = sort(allDistsControl);

remapMetaMeanControl = repmat(metaMeanControl,1,size(allMeanControls,2));
allControlDiffs = allMeanControls - remapMetaMeanControl;
allGeneDiffs = allMeanGenes - remapMetaMeanControl;

L1DistControls =  sum(allControlDiffs,1);
L1DistGene =  sum(allGeneDiffs,1);
% absDistControls =  sum(abs(allControlDiffs),1);
% absDistdistGene =  sum(abs(allGeneDiffs),1);

loggerFname = [mainDirname 'controlInterdayAssessment/log_controlAssessment_' propertyStr '.txt'];
logger = fopen(loggerFname,'w');


fprintf(logger,'    Distance from mean control (L1 control, L1 KD)    \n');
fprintf(logger,'*************\n');
for i = 1 : iDayGeneSeq
    curInd = sortedInds(i);        
    fprintf(logger,sprintf('%s: %.1f (%.1f, %.1f)\n',geneDayDiff{curInd}.dayGeneSeqStr,sortedDists(i),L1DistControls(curInd),L1DistGene(curInd)));
end
fprintf(logger,'*************\n');

fclose(logger);

metaMeanHealingControl = mean(healingRateControl);
controlHealingDists = healingRateControl - metaMeanHealingControl;
geneHealingDists = healingRateGene - metaMeanHealingControl;

% Distribution of mean control healing rate
[nelements,healingRateControlCenters] = hist(controlHealingDists,20); 
healingRateControlDistribution = nelements./sum(nelements);
h = figure;
hold on;
bar(healingRateControlCenters,healingRateControlDistribution,'r');
xlabel('Relative healing rate (\mum hour{-1})','FontSize',22);
ylabel('Percent','FontSize',22);
haxes = get(h,'CurrentAxes');
% set(haxes,'XLim',[0,90]);
% set(haxes,'XTick',0:45:90);
% set(haxes,'XTickLabel',0:45:90);
% set(haxes,'YLim',[0,0.4]);
% set(haxes,'YTick',0:0.1:0.4);
% set(haxes,'YTickLabel',0:0.1:0.4);
set(haxes,'FontSize',22);
set(h,'Color','none');
hold off;
export_fig([mainDirname 'controlInterdayAssessment/controlDiffDistribution.eps']);


%% visualize: distance from mean control: control vs. gene
doVisualize(L1DistControls,L1DistGene,geneDayDiff,targetGenesStr,[mainDirname 'controlInterdayAssessment/' propertyStr '_relativeToMeanControlL1'],propertyStr);
% doVisualize(absDistControls,absDistdistGene,geneDayDiff,targetGenesStr,[mainDirname 'controlInterdayAssessment/' propertyStr '_relativeToMeanControlAbs'],propertyStr);
doVisualize(allDistsControl,allDistsGene,geneDayDiff,targetGenesStr,[mainDirname 'controlInterdayAssessment/' propertyStr '_relativeToMeanControlL2'],propertyStr);
doVisualize(controlHealingDists,geneHealingDists,geneDayDiff,targetGenesStr,[mainDirname 'controlInterdayAssessment/' propertyStr '_relativeHealing'],propertyStr);

end
%%
function [] = doVisualize(controlx,geney,geneDayDiff,targetGenesStr,fnamePrefix,propertyStr)
nTargets = length(targetGenesStr);
nConditions = nTargets + 3; % 0% KD, CDC42/RAC1/beta-PIX, > 50% KD (rest)

cmap = colormap(hsv(nConditions));
fontsize = 24;
markerSize = 6;
LineWidth = 3;


[negCtrlInds,restInds,targetsInds,posCntrl] = whGetTargetInds(geneDayDiff,targetGenesStr);

h = figure;
xlabel('Control','FontSize',fontsize);
ylabel('Gene','FontSize',fontsize);
hold on;

plot(controlx(negCtrlInds),geney(negCtrlInds),'o','MarkerEdgeColor',cmap(1,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','0% KD');   
plot(controlx(restInds),geney(restInds),'o','MarkerEdgeColor',cmap(2,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','> 50% KD');   
for t = 1 : nTargets
    plot(controlx(targetsInds{t}),geney(targetsInds{t}),'o','MarkerEdgeColor',cmap(t+2,:),'LineWidth',LineWidth,'MarkerSize',markerSize+3,'DisplayName',targetGenesStr{t});
end
plot(controlx(posCntrl.inds),geney(posCntrl.inds),'o','MarkerEdgeColor',cmap(nConditions,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','Pos Ctrl');

haxes = get(h,'CurrentAxes');

xlim = get(haxes,'xlim');

set(haxes,'FontSize',fontsize);

legend('show','Location','NorthEastOutside');

plot([xlim(1),xlim(2)],[xlim(1),xlim(2)],'-k','LineWidth',2);
plot(0,0,'*k','LineWidth',2,'MarkerSize',10); 

set(h,'Color','none');

outFname = [fnamePrefix '_legend.eps'];
export_fig(outFname);

legend off;
outFname = [fnamePrefix '.eps'];
export_fig(outFname);

hold off;
end