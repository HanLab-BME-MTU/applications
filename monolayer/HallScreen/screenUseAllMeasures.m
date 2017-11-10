function [] = screenUseAllMeasures()

clc;

dname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ScreenFinal/kd50time200_RHOA/screen/';

propertyStr = {'Speed','Directionality'};
measureStr = {'PC1','PC2','PC3'};

posContMeasures = [];
offtargetContMeasures = [];
screenMeasures = [];
for ip = 1 : length(propertyStr)
    for im = 1 : length(measureStr)
        load([dname propertyStr{ip} '_' measureStr{im}]);
        posContMeasures = [posContMeasures, sumZScorePosControl'];
        offtargetContMeasures = [offtargetContMeasures, sumZScoreKD0'];
        screenMeasures = [screenMeasures, sumZScoreScreen'];
    end
end

% sumPosCont = sum(abs(posContMeasures),2);
% sumOffCont = sum(abs(offtargetContMeasures),2);
% sumScreen = sum(abs(screenMeasures),2);

sumPosCont = sum(posContMeasures,2);
sumOffCont = sum(offtargetContMeasures,2);
sumScreen = sum(screenMeasures,2);

loggerFname = [dname 'allMeasures/logScreenMeasures.txt'];
logger = fopen(loggerFname,'w');

%% Selecting a threshold based on positive and off-target controls
nBootstrap = 1000;
allThresholds = estimateZScoreThresholdDistribution(sumZScorePosControl,sumZScoreKD0,nBootstrap);

mean(allThresholds)
std(allThresholds)
max(allThresholds)
prctile(allThresholds,99)
prctile(allThresholds,99.99)


% zRange = 1 : 0.1 : 10;
% nZRange = length(zRange);
% fMeasure = nan(1,nZRange);
% 
% for zi = 1 : nZRange
%     z = zRange(zi);
%     nHitsReal = sum(sumZScorePosControl > z);
%     nHitsDetected = nHitsReal + sum(sumZScoreKD0 > z);
%     precision = nHitsReal / nHitsDetected;
%     recall = nHitsReal / length(sumZScorePosControl);
%     fMeasure(zi) = 2*(precision * recall) / (precision + recall);
% end
% 
% maxFmeausre = max(fMeasure);
% maxInd = find(fMeasure == maxFmeausre,1,'last'); % last
% disp(zRange(maxInd));
% 
% % plot figure
% FPosition = [0 0 300 300];
% APosition = [0.2 0.2 0.75 0.75];
% fontsize = 10;
% h = figure;
% xlabel('Z-score','FontSize',fontsize);
% ylabel('F-Measure','FontSize',fontsize);
% hold on;
% 
% plot(zRange,fMeasure,'--','Color','k','LineWidth',2);
% 
% haxes = get(h,'CurrentAxes');
% axis(haxes);
% set(h,'Color','none'); % set(h,'Color','w');
% axisHandle= findobj(h,'type','axes');
% set(haxes,'XLim',[0,15]);
% set(haxes,'XTick',0:5:15);
% set(haxes,'XTickLabel',0:5:15);
% set(haxes,'YLim',[0,1]);
% set(haxes,'YTick',0:0.5:1);
% set(haxes,'YTickLabel',0:0.5:1);
% set(haxes,'FontSize',fontsize);
% set(h,'Position',FPosition,'PaperPositionMode','auto');
% set(axisHandle,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
% set(get(axisHandle,'XLabel'),'FontSize',fontsize); set(get(axisHandle,'YLabel'),'FontSize',fontsize);
% export_fig([mainDirname 'screen/' propertyStr '_' representativeStr '_sumZScoreFmeasure.eps']);
% hold off;
% close all;

%%


nPos = length(sumPosCont);
nOff = length(sumOffCont);
nScreen = length(sumScreen);

[sortedScreen,sortedIndsScreen] = sort(sumScreen);
sortedStrScreen = strScreen(sortedIndsScreen);
% screenIndsSorted = screenInds(sortedIndsZScoreScreen);

for igene = nScreen : -1 : 1
    fprintf(logger,sprintf('%s: %.1f \n',...
        sortedStrScreen{igene},sortedScreen(igene)));    
end

fclose(logger);
end
    
%     %% plot (previously detected hits)
%     notValidatedHitsStr = {'ARHGEF10','TRIO','TUBA','ARHGEF9','DOCK10'};
%     [notValidatedHitsInds,tmpInds] = getHitsInds(notValidatedHitsStr,strScreenGene,screenInds,dayGeneData);
%     notValidatedHitsInds = [notValidatedHitsInds{:}];
%     
%     hitsStr = {'SOS1','ARHGEF18','ARHGEF3','ARHGEF11','ARHGEF28'};
%     [hitsInds,restInds] = getHitsInds(hitsStr,strScreenGene,screenInds,dayGeneData);
%     outFname = [mainDirname 'screen/zScoreHits_' propertyStr '_' representativeStr '.eps'];
%     plotScreenZScore(sumZScoreKD0,sumZScorePosControl,sumZScoreScreen,hitsStr,hitsInds,restInds,notValidatedHitsInds,outFname);          

function [] = plotScreenZScore(sumZScoreKD0,sumZScorePosControl,sumZScoreScreen,hitsStr,hitsInds,restInds,notValidatedHitsInds,outFname)

restExcludeNotValidatedInds = restInds;
restExcludeNotValidatedInds(ismember(restInds,notValidatedHitsInds)) = [];

sumZScoreRest = sumZScoreScreen(restExcludeNotValidatedInds);
sumZScoreNotValidated = sumZScoreScreen(notValidatedHitsInds);

% Trancate Z-score by 15!
 sumZScorePosControl(sumZScorePosControl > 15) = 15; 
 sumZScoreKD0(sumZScoreKD0 > 15) = 15; 
 sumZScoreRest(sumZScoreRest > 15) = 15;
 sumZScoreNotValidated(sumZScoreNotValidated > 15) = 15;
 sumZScoreScreen(sumZScoreScreen > 15) = 15;  

nKD0 = length(sumZScoreKD0);
nPosControl = length(sumZScorePosControl);
nRest = length(restExcludeNotValidatedInds);
nNotValidated = length(notValidatedHitsInds);

allZScores = [sumZScorePosControl,sumZScoreKD0,sumZScoreRest,sumZScoreNotValidated];
allNs = [nPosControl,nKD0,nRest,nNotValidated];

nHits = length(hitsInds);

nsHits = zeros(1,nHits);
zScoreHits = cell(1,nHits);
for ihit = 1 : nHits
    curInds = hitsInds{ihit};
    nsHits(ihit) = length(curInds);
    allNs = [allNs nsHits(ihit)];
    curZscores = sumZScoreScreen(curInds);    
    zScoreHits{ihit} = curZscores;
    allZScores = [allZScores curZscores];
end

assert(max(allZScores) <= 15);

ns = [1 cumsum(allNs)];
maxZscore = min(15.5,max(allZScores) + 0.5);

N = ns(end);


legendStrs = [{'Pos Cntl','KD = 0','Screen'},hitsStr];

cmap = colormap(hsv(nHits+3));
cmap = cmap([1,3,6,2,4,7,5,8],:);

% [left bottom width height]
FPosition = [0 0 900 300];
APosition = [0.1 0.2 0.7 0.75]; 

fontsize = 10;

h = figure;
xlabel('Experiment','FontSize',fontsize);
ylabel('Z-score','FontSize',fontsize);
hold on;
plot(ns(1):ns(2),sumZScorePosControl,'o','MarkerEdgeColor',cmap(1,:),'LineWidth',2,'MarkerSize',6);
plot((ns(2)+1):ns(3),sumZScoreKD0,'o','MarkerEdgeColor',cmap(2,:),'LineWidth',2,'MarkerSize',6);
plot((ns(3)+1):ns(4),sumZScoreRest,'o','MarkerEdgeColor',cmap(3,:),'LineWidth',2,'MarkerSize',6);
for i = 4 : nHits+3
    %     plot((ns(i)+1):ns(i+1),zScoreHits{i-3},'o','MarkerEdgeColor',cmap(i,:),'LineWidth',2,'MarkerSize',6);
    plot((ns(i+1)+1):ns(i+2),zScoreHits{i-3},'o','MarkerEdgeColor',cmap(i,:),'LineWidth',2,'MarkerSize',6);
end
legend(legendStrs,'Location','eastoutside');
plot((ns(4)+1):ns(5),sumZScoreNotValidated,'X','MarkerEdgeColor',cmap(3,:),'LineWidth',2,'MarkerSize',6);
plot([1,N],[10,10],'--k','LineWidth',2);
plot([1,N],[-10,-10],'--k','LineWidth',2);
haxes = get(h,'CurrentAxes');
% set(haxes,'XLim',[-3,(ns(7)+4)]);
set(haxes,'XLim',[-3,(N+4)]);
set(haxes,'YLim',[-maxZscore,maxZscore]);
set(haxes,'XTick',0:50:N);
set(haxes,'XTickLabel',0:50:N);
set(haxes,'YTick',-15:5:15);
set(haxes,'YTickLabel',-15:5:15);
set(haxes,'FontSize',fontsize);
set(h,'Color','w');
set(h,'Position',FPosition,'PaperPositionMode','auto');
axisHandle= findobj(h,'type','axes');
set(axisHandle,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
set(get(axisHandle,'XLabel'),'FontSize',fontsize); set(get(axisHandle,'YLabel'),'FontSize',fontsize);
export_fig([outFname(1:end-4) '_legend.eps']);
legend off;

FPosition = [0 0 700 300];
APosition = [0.1 0.2 0.85 0.75]; 
set(h,'Position',FPosition,'PaperPositionMode','auto');
axisHandle= findobj(h,'type','axes');
set(axisHandle,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
set(get(axisHandle,'XLabel'),'FontSize',fontsize); set(get(axisHandle,'YLabel'),'FontSize',fontsize);
export_fig(outFname);
hold off;

end


%%
function [hitsInds,restInds] = getHitsInds(hitsStr,strScreenGene,screenInds,dayGeneData)
nHits = length(hitsStr);
hitsInds = cell(1,nHits);
rest = true(1,length(screenInds));
for ihit = 1 : nHits
    curHitStr = hitsStr{ihit};
    curHitInds = find(strcmp(curHitStr,strScreenGene));
    curInds = screenInds(curHitInds');
    rest(curHitInds') = false;
    % validation
    hitsInds{ihit} = curHitInds; % NOTE, these are not the original dayGeneData data!
    for ivalid = 1 : length(hitsInds{ihit})
        assert(strcmp(curHitStr,dayGeneData{curInds(ivalid)}.geneStr));
    end
end
restInds = find(rest);
end

%%
function [] = displayBoxPlot(data,labels,ylabelStr,outFname,dataLim)
whiskerParam = 0.7193;
colors = 'rbg';
fontsize = 22;
h = figure;
hold on;
if nargin == 4    
    boxplot(data,labels,'whisker',whiskerParam,'OutlierSize',6); % 0.7193 --> 90%
else
    boxplot(data,labels,'whisker',whiskerParam,'DataLim',dataLim,'OutlierSize',6); % 0.7193 --> 90%
end
haxes = get(h,'CurrentAxes');
set(gca,'FontSize',fontsize);
text_h = findobj(gca, 'Type', 'text');
for cnt = 1:length(text_h)
    set(text_h(cnt),'FontSize', fontsize,'VerticalAlignment', 'top','HorizontalAlignment', 'center');
end
% rotateticklabel(haxes,45);
h1 = findobj(gca,'Tag','Box');
h2 = findobj(gca,'Tag','Upper Whisker');
h3 = findobj(gca,'Tag','Lower Whisker');
h4 = findobj(gca,'Tag','Median');
h5 = findobj(gca,'Tag','Upper Adjacent Value');
h6 = findobj(gca,'Tag','Lower Adjacent Value');
for j=1:length(h1)
    patch(get(h1(j),'XData'),get(h1(j),'YData'),colors(j),'FaceAlpha',.5,'LineWidth',2);
    patch(get(h2(j),'XData'),get(h2(j),'YData'),colors(j),'FaceAlpha',.5,'LineWidth',2);
    patch(get(h3(j),'XData'),get(h3(j),'YData'),colors(j),'FaceAlpha',.5,'LineWidth',2);
    patch(get(h4(j),'XData'),get(h4(j),'YData'),colors(j),'FaceAlpha',.5,'LineWidth',2);
    patch(get(h5(j),'XData'),get(h5(j),'YData'),colors(j),'FaceAlpha',.5,'LineWidth',2);
    patch(get(h6(j),'XData'),get(h6(j),'YData'),colors(j),'FaceAlpha',.5,'LineWidth',2);
end

ol = findobj(gca,'Tag','Outliers');
for j = 1 : length(ol)
    set(ol(j),'LineWidth',2);
    %     patch(get(ol(j),'XData'),get(h1(j),'YData'),colors(j),'FaceAlpha',.5,'LineWidth',2);
end

ylim([dataLim(1)-1,dataLim(2)+1]);
ylabel(ylabelStr);
% oh=findobj(gca,'tag','Outliers');
% set(oh,'Visible','off');
set(h,'Color','w');
% position = get(h,'position');
% set(h,'position',[position(1:2) round(1.5*position(3:4))]);
hold off;
export_fig(outFname);
end

%%
function allThresholds = estimateZScoreThresholdDistribution(sumZScorePosControl,sumZScoreOfftargetControl,nBootstrap)

nPosControl = length(sumZScorePosControl);
nOfftargetControl = length(sumZScoreOfftargetControl);

zRange = 1 : 0.1 : 60;
nZRange = length(zRange);

allThresholds = nan(1,nBootstrap);

for iBoot = 1 : nBootstrap
    
    fMeasure = nan(1,nZRange);
    
    sumZScorePosControl_sample = datasample(sumZScorePosControl,nPosControl);
    sumZScoreOfftargetControl_sample = datasample(sumZScoreOfftargetControl,nOfftargetControl);
    
    for zi = 1 : nZRange
        z = zRange(zi);
        nHitsReal = sum(sumZScorePosControl_sample > z);
        nHitsDetected = nHitsReal + sum(sumZScoreOfftargetControl_sample > z);
        if nHitsDetected == 0
            fMeasure(zi) = 0;
        else
            precision = nHitsReal / nHitsDetected;
            recall = nHitsReal / length(sumZScorePosControl_sample);
            fMeasure(zi) = 2*(precision * recall) / (precision + recall);
        end
    end
    
    maxFmeausre = max(fMeasure);
    maxInd = find(fMeasure == maxFmeausre,1,'last'); % last
    allThresholds(iBoot) = zRange(maxInd);
end
end