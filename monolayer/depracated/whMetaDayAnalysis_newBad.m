function [] = whMetaDayAnalysis(allFeatures,strLabels,metaData,mainDirname)

warning('off','all');

% addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/algs/libsvm-3.18'),'-begin');
addpath(genpath('/apps/MATLAB/R2013a/toolbox/stats/stats'));

dayAnalysisProperty(allFeatures.speedFeats.features,strLabels,metaData,mainDirname,'Speed');
% dayAnalysisProperty(allFeatures.directionalityFeats.features,strLabels,metaData,mainDirname,'Directionality');
% % dayAnalysisProperty(allFeatures.strainRateFeats.features,strLabels,metaData,mainDirname,'StrainRate');
% dayAnalysisProperty(allFeatures.coordinationFeats.features,strLabels,metaData,mainDirname,'Coordination');
close all;
end


% TODO: update by day
function [] = dayAnalysisProperty(features,strLabels,metaData,mainDirname,propertyStr)
close all;

nFeats = size(features,1);
N = length(metaData.groupsByDays);
n = length(strLabels);

strPSup = 'pSuper';
strNT = 'NT';
indsPSup = strcmp(whToPrefix(strLabels),strPSup) | strcmp(whToPrefix(strLabels),strNT); % then filter by days for the gene

allDiffVectors = [];
allDiffToMeanControl = [];
allDiffMeanGeneToMeanControl = [];
geneDayDiff = {};

nGeneDay = 0;
for day = 1 : N
    dayStr = metaData.groupsByDays{day}.dates;
    daysInds = metaData.groupsByDays{day}.inds;
    indsDay = false(1,n);
    indsDay(daysInds) = true;
    indsPSup1 = indsPSup & indsDay;
    
    nControl = sum(indsPSup1);
    
    if nControl < 2
        continue;
    end    
    
    featsControl = features(:,indsPSup1);
    meanControl = mean(featsControl,2);
    
    % Now add the gene for that day (make sure there is only 1 of those)
    indsGenes = indsDay & ~indsPSup;
    
    if sum(indsGenes) < 4
        continue;
    end
    
    genesStr = metaData.treatment(indsGenes);
    
    geneStr = unique(whToPrefix(genesStr));        
    
    nGene = length(geneStr);
    
    for curGene = 1 : nGene
        nGeneDay = nGeneDay + 1;
        curGeneStr = geneStr{curGene};
        indsGene = strcmp(whToPrefix(strLabels),curGeneStr);
        indsGene = indsGene & indsDay;
        nCurGene = sum(indsGene);
        
        featsGene = features(:,indsGene);     
        meanGene = mean(featsGene,2);
        
        diffMeanGeneToMeanControl = meanGene - meanControl;
        
        % TODO: define distance gene-control data structure.
        % For each gene includes matrix of distance-vectors from day's
        % control
        diffGene = zeros(nFeats,nCurGene*nControl);
        diffGeneToMeanControl = zeros(nFeats,nCurGene);
        for igene = 1 : nCurGene
            diffGeneToMeanControl(:,igene) = featsGene(:,igene) - meanControl;
            for icontrol = 1 : nControl
                diffGene(:,(igene-1)*nControl+icontrol) = featsGene(:,igene) - featsControl(:,icontrol);
            end
        end
        allDiffVectors = [allDiffVectors diffGene];
        allDiffToMeanControl = [allDiffToMeanControl diffGeneToMeanControl];
        allDiffMeanGeneToMeanControl = [allDiffMeanGeneToMeanControl diffMeanGeneToMeanControl];
        %         allDiffLabels{nGeneDay} = curGeneStr;
        
        geneDayDiff{nGeneDay}.geneStr = curGeneStr;
        geneDayDiff{nGeneDay}.dayStr = dayStr;
        geneDayDiff{nGeneDay}.diff = diffGene;
        geneDayDiff{nGeneDay}.diffInds = (size(allDiffVectors,2)-size(diffGene,2)+1):size(allDiffVectors,2);
        geneDayDiff{nGeneDay}.geneFeatures = featsGene;
        geneDayDiff{nGeneDay}.controlFeatures = featsControl;
        geneDayDiff{nGeneDay}.geneInds = indsGene;
        geneDayDiff{nGeneDay}.controlInds = indsPSup1;
        geneDayDiff{nGeneDay}.nGeneFeatures = nCurGene;
        geneDayDiff{nGeneDay}.nControlFeatures = nControl;
        geneDayDiff{nGeneDay}.meanControl = meanControl;
        geneDayDiff{nGeneDay}.diffToMeanControl = diffGeneToMeanControl;
        geneDayDiff{nGeneDay}.diffToMeanControlInds = (size(allDiffToMeanControl,2)-size(diffGeneToMeanControl,2)+1):size(allDiffToMeanControl,2);
        geneDayDiff{nGeneDay}.diffMeanGeneToMeanControl = diffGeneToMeanControl;
        geneDayDiff{nGeneDay}.diffMeanGeneToMeanControlInds = size(allDiffMeanGeneToMeanControl,2);
        fprintf(sprintf('%s: gene: %d, control: %d\n',geneDayDiff{nGeneDay}.geneStr,geneDayDiff{nGeneDay}.nGeneFeatures,geneDayDiff{nGeneDay}.nControlFeatures))
    end
    
    close all;        
end

[coeff,score,latent] = pca(allDiffVectors');
accVariance = cumsum(latent)./sum(latent);

% TODO: now take the geneDayDiff per gene & day --> use the PCA to
% transform. Check how the PCA looks.

% % Get the 7 independent controls
% ngroups = length(metaData.groupsByTreatments);
% for i = 1 : ngroups
%     curGroup = metaData.groupsByTreatments{i};
%     if strcmp(curGroup.treatment,'Control')
%         interDayControlInds = curGroup.inds;
%     end
% end

% Calculate day variance of control + GEF
% % dayVariance(features,geneDayDiff,interDayControlInds,indsPSup,mainDirname,propertyStr);
% dayVariance(features,geneDayDiff,indsPSup,mainDirname,propertyStr);
% dayAnalysis(allDiffVectors,allDiffToMeanControl,allDiffMeanGeneToMeanControl,geneDayDiff,mainDirname,propertyStr);
% dayClassification(features,geneDayDiff,mainDirname,propertyStr);
dayClustering(geneDayDiff,allDiffVectors,allDiffMeanGeneToMeanControl,mainDirname,propertyStr);
end


%%
% function [] = dayVariance(features,geneDayDiff,interDayControlInds,indsPSup,mainDirname,propertyStr)
function [] = dayVariance(features,geneDayDiff,indsPSup,mainDirname,propertyStr)
[coeff,score,latent] = pca(features');
accVariance = cumsum(latent)./sum(latent);

nGeneDay = length(geneDayDiff);

varControl = zeros(1,nGeneDay);
varGene = zeros(1,nGeneDay);


cmap = colormap(hsv(nGeneDay));
fontsize = 24;

h = figure;
xlabel('Control variance','FontSize',fontsize);
ylabel('Gene variance','FontSize',fontsize);
hold on;

varInterDayPSup = var(score(indsPSup,1));
% varInterDayControl = var(score(interDayControlInds,1));
generalVariance = var(score(:,1));

for iGeneDay = 1 : nGeneDay
    varControl(iGeneDay) = var(score(geneDayDiff{iGeneDay}.controlInds,1));
    varGene(iGeneDay) = var(score(geneDayDiff{iGeneDay}.geneInds,1));
%     plot(varControl(iGeneDay),varGene(iGeneDay),sprintf('%s',markersPerm(ceil(iGeneDay/nColors))),'MarkerEdgeColor','k','MarkerFaceColor',sprintf('%s',colorsPerm(mod(iGeneDay,nColors)+1)),'MarkerSize',10,...
%         'DisplayName',geneDayDiff{iGeneDay}.geneStr);   
    plot(varControl(iGeneDay),varGene(iGeneDay),'o','MarkerEdgeColor',cmap(iGeneDay,:),'LineWidth',2,'MarkerSize',7,'DisplayName',[geneDayDiff{iGeneDay}.geneStr '_' geneDayDiff{1}.dayStr]);   
end

legend('Location','EastOutside');

haxes = get(h,'CurrentAxes');
if strcmp(propertyStr,'Speed')
    maxVar = 2100;
    assert(max([varControl varGene]) < maxVar);
    set(haxes,'XLim',[0,maxVar]);
    set(haxes,'XTick',0:500:maxVar);
    set(haxes,'XTickLabel',0:500:maxVar);
    set(haxes,'YLim',[0,maxVar]);
    set(haxes,'YTick',0:500:maxVar);
    set(haxes,'YTickLabel',0:500:maxVar);
else if strcmp(propertyStr,'Directionality')
        maxVar = 13;
        assert(max([varControl varGene]) < maxVar);
        set(haxes,'XLim',[0,maxVar]);
        set(haxes,'XTick',0:4:maxVar);
        set(haxes,'XTickLabel',0:4:maxVar);
        set(haxes,'YLim',[0,maxVar]);
        set(haxes,'YTick',0:4:maxVar);
        set(haxes,'YTickLabel',0:4:maxVar);
    else if strcmp(propertyStr,'Coordination')
            maxVar = 0.3;
            assert(max([varControl varGene]) < maxVar);
            set(haxes,'XLim',[0,maxVar]);
            set(haxes,'XTick',0:0.1:maxVar);
            set(haxes,'XTickLabel',0:0.1:maxVar);
            set(haxes,'YLim',[0,maxVar]);
            set(haxes,'YTick',0:0.1:maxVar);
            set(haxes,'YTickLabel',0:0.1:maxVar);
        end
    end
end

set(haxes,'FontSize',fontsize);

set(h,'Color','none');


plot([0,maxVar],[0,maxVar],'--k','LineWidth',3);
plot([varInterDayPSup,varInterDayPSup],[0,maxVar],'--g','LineWidth',3);
% plot([varInterDayControl,varInterDayControl],[0,maxVar],'--b','LineWidth',3);
plot([0,maxVar],[generalVariance,generalVariance],'--r','LineWidth',3);
plot([generalVariance,generalVariance],[0,maxVar],'--r','LineWidth',3);

controlGeneVarFname = [mainDirname 'variance/controlGeneDayVariance' propertyStr '_legend.eps'];
export_fig(controlGeneVarFname);

legend off;
controlGeneVarFname = [mainDirname 'variance/controlGeneDayVariance' propertyStr '.eps'];
export_fig(controlGeneVarFname);

hold off;

%% Daily control PCs
h = figure;
xlabel('PC 1','FontSize',fontsize);
ylabel('PC 2','FontSize',fontsize);
hold on;

% plot(score(interDayControlInds,1),score(interDayControlInds,2),'o','MarkerEdgeColor','k','LineWidth',2,'MarkerSize',7,'DisplayName','Control');

for iGeneDay = 1 : nGeneDay
    plot(score(geneDayDiff{iGeneDay}.controlInds,1),score(geneDayDiff{iGeneDay}.controlInds,2),'o','MarkerEdgeColor',cmap(iGeneDay,:),'LineWidth',2,'MarkerSize',7,'DisplayName',geneDayDiff{iGeneDay}.dayStr);
end

haxes = get(h,'CurrentAxes');
if strcmp(propertyStr,'Speed')
    assert(max(score(:,1)) < 100 && min(score(:,1)) > -100);
    assert(max(score(:,2)) < 50 && min(score(:,2)) > -50);
    set(haxes,'XLim',[-100,100]);
    set(haxes,'XTick',-100:50:100);
    set(haxes,'XTickLabel',-100:50:100);
    set(haxes,'YLim',[-50,50]);
    set(haxes,'YTick',-50:25:50);
    set(haxes,'YTickLabel',-50:25:50);
else if strcmp(propertyStr,'Directionality')
        assert(max(score(:,1)) < 8 && min(score(:,1)) > -8);
        assert(max(score(:,2)) < 7.5 && min(score(:,2)) > -3.2);
        set(haxes,'XLim',[-8,8]);
        set(haxes,'XTick',-8:4:8);
        set(haxes,'XTickLabel',-8:4:8);
        set(haxes,'YLim',[-3.2,7.5]);
        set(haxes,'YTick',-3:3:6);
        set(haxes,'YTickLabel',-3:3:6);
    else if strcmp(propertyStr,'Coordination')
            assert(max(score(:,1)) < 1.25 && min(score(:,1)) > -1);
            assert(max(score(:,2)) < 0.65 && min(score(:,2)) > -0.5);
            set(haxes,'XLim',[-1,1.25]);
            set(haxes,'XTick',-1:0.5:1);
            set(haxes,'XTickLabel',-1:0.5:1);
            set(haxes,'YLim',[-0.5,0.65]);
            set(haxes,'YTick',-0.5:0.25:0.5);
            set(haxes,'YTickLabel',-0.5:0.25:0.5);
        end
    end
end

set(haxes,'FontSize',fontsize);

set(h,'Color','none');

legend('Location','EastOutside');

controlGeneVarFname = [mainDirname 'variance/pcaPSup_' propertyStr '_legend.eps'];
export_fig(controlGeneVarFname);

legend off;
controlGeneVarFname = [mainDirname 'variance/pcaPSup_' propertyStr '.eps'];
export_fig(controlGeneVarFname);

hold off;
end
%% 
function [] = dayAnalysis(allDiffVectors,allDiffToMeanControl,allDiffMeanGeneToMeanControl,geneDayDiff,mainDirname,propertyStr)
[coeff,score,latent] = pca(allDiffVectors');
accVariance = cumsum(latent)./sum(latent);

nGeneDay = length(geneDayDiff);

cmap = colormap(hsv(nGeneDay));

% colorsPerm = 'ymcrgbk'; 
% nColors = length(colorsPerm);
% markersPerm ='os^Vph><+*X';
fontsize = 24;

%% PCA: all gene x pSup
h = figure;
xlabel('PC 1','FontSize',fontsize);
ylabel('PC 2','FontSize',fontsize);
hold on;

for iGeneDay = 1 : nGeneDay     
    %     plot(score(geneDayDiff{iGeneDay}.diffInds,1),score(geneDayDiff{iGeneDay}.diffInds,2),sprintf('%s',markersPerm(ceil(iGeneDay/nColors))),'MarkerEdgeColor','k','MarkerFaceColor',sprintf('%s',colorsPerm(mod(iGeneDay,nColors)+1)),'MarkerSize',10,...
    %         'DisplayName',geneDayDiff{iGeneDay}.geneStr);
    plot(score(geneDayDiff{iGeneDay}.diffInds,1),score(geneDayDiff{iGeneDay}.diffInds,2),'o','MarkerEdgeColor',cmap(iGeneDay,:),'LineWidth',2,'MarkerSize',7,'DisplayName',geneDayDiff{iGeneDay}.geneStr);   
end

haxes = get(h,'CurrentAxes');
if strcmp(propertyStr,'Speed')
    assert(max(abs(score(:,1))) < 150);
    assert(max(abs(score(:,2))) < 100);
    set(haxes,'XLim',[-150,150]);
    set(haxes,'XTick',-150:75:150);
    set(haxes,'XTickLabel',-150:75:150);
    set(haxes,'YLim',[-100,100]);
    set(haxes,'YTick',-100:50:100);
    set(haxes,'YTickLabel',-100:50:100);
else if strcmp(propertyStr,'Directionality')
        assert(max(abs(score(:,1))) < 12);
        assert(max(abs(score(:,2))) < 5);
        set(haxes,'XLim',[-12,12]);
        set(haxes,'XTick',-12:6:12);
        set(haxes,'XTickLabel',-12:6:12);
        set(haxes,'YLim',[-5,5]);
        set(haxes,'YTick',-5:5:5);
        set(haxes,'YTickLabel',-5:5:5);
    else if strcmp(propertyStr,'Coordination')
            assert(max(abs(score(:,1))) < 1.6);
            assert(max(abs(score(:,2))) < 0.62);
            set(haxes,'XLim',[-1.6,1.6]);
            set(haxes,'XTick',-1.5:1.5:1.5);
            set(haxes,'XTickLabel',-1.5:1.5:1.5);
            set(haxes,'YLim',[-0.62,0.62]);
            set(haxes,'YTick',-0.6:0.3:0.6);
            set(haxes,'YTickLabel',-0.6:0.3:0.6);
        end
    end
end

set(haxes,'FontSize',fontsize);

set(h,'Color','none');

legend('Location','EastOutside');

plot(0,0,'*k','LineWidth',4,'MarkerSize',20);

controlGeneVarFname = [mainDirname 'day/pca_' propertyStr '_legend.eps'];
export_fig(controlGeneVarFname);

legend off;
controlGeneVarFname = [mainDirname 'day/pca_' propertyStr '.eps'];
export_fig(controlGeneVarFname);

hold off;


%% PCA: each gene x pSup

for iGeneDay = 1 : nGeneDay
    h = figure;
    title(geneDayDiff{iGeneDay}.geneStr,'FontSize',fontsize);
    xlabel('PC 1','FontSize',fontsize);
    ylabel('PC 2','FontSize',fontsize);
    hold on;
    plot(score(geneDayDiff{iGeneDay}.diffInds,1),score(geneDayDiff{iGeneDay}.diffInds,2),'o','MarkerEdgeColor',cmap(iGeneDay,:),'LineWidth',2,'MarkerSize',8);
    
    haxes = get(h,'CurrentAxes');
    if strcmp(propertyStr,'Speed')
        assert(max(abs(score(:,1))) < 150);
        assert(max(abs(score(:,2))) < 100);
        set(haxes,'XLim',[-150,150]);
        set(haxes,'XTick',-150:75:150);
        set(haxes,'XTickLabel',-150:75:150);
        set(haxes,'YLim',[-100,100]);
        set(haxes,'YTick',-100:50:100);
        set(haxes,'YTickLabel',-100:50:100);
    else if strcmp(propertyStr,'Directionality')
            assert(max(abs(score(:,1))) < 12);
            assert(max(abs(score(:,2))) < 5);
            set(haxes,'XLim',[-12,12]);
            set(haxes,'XTick',-12:6:12);
            set(haxes,'XTickLabel',-12:6:12);
            set(haxes,'YLim',[-5,5]);
            set(haxes,'YTick',-5:5:5);
            set(haxes,'YTickLabel',-5:5:5);
        else if strcmp(propertyStr,'Coordination')
                assert(max(abs(score(:,1))) < 1.6);
                assert(max(abs(score(:,2))) < 0.62);
                set(haxes,'XLim',[-1.6,1.6]);
                set(haxes,'XTick',-1.5:1.5:1.5);
                set(haxes,'XTickLabel',-1.5:1.5:1.5);
                set(haxes,'YLim',[-0.62,0.62]);
                set(haxes,'YTick',-0.6:0.3:0.6);
                set(haxes,'YTickLabel',-0.6:0.3:0.6);
            end
        end
    end
    
    set(haxes,'FontSize',fontsize);
    
    set(h,'Color','none');
    
    legend();
    
    plot(0,0,'*k','LineWidth',4,'MarkerSize',20);
    
    controlGeneVarFname = [mainDirname 'day/gene/' geneDayDiff{iGeneDay}.geneStr '_' propertyStr '.eps'];
    export_fig(controlGeneVarFname);
    
    hold off;
    close all;
end

%% mean(gene) - mean(pSup) - averaging the PCA of all data!
h = figure;
xlabel('PC 1','FontSize',fontsize);
ylabel('PC 2','FontSize',fontsize);
hold on;

for iGeneDay = 1 : nGeneDay     
    plot(mean(score(geneDayDiff{iGeneDay}.diffInds,1)),mean(score(geneDayDiff{iGeneDay}.diffInds,2)),'o','MarkerEdgeColor',cmap(iGeneDay,:),'LineWidth',2,'MarkerSize',7,'DisplayName',geneDayDiff{iGeneDay}.geneStr);
end

haxes = get(h,'CurrentAxes');
if strcmp(propertyStr,'Speed')
    assert(max(abs(score(:,1))) < 150);
    assert(max(abs(score(:,2))) < 100);
    set(haxes,'XLim',[-150,150]);
    set(haxes,'XTick',-150:75:150);
    set(haxes,'XTickLabel',-150:75:150);
    set(haxes,'YLim',[-100,100]);
    set(haxes,'YTick',-100:50:100);
    set(haxes,'YTickLabel',-100:50:100);
else if strcmp(propertyStr,'Directionality')
        assert(max(abs(score(:,1))) < 12);
        assert(max(abs(score(:,2))) < 5);
        set(haxes,'XLim',[-12,12]);
        set(haxes,'XTick',-12:6:12);
        set(haxes,'XTickLabel',-12:6:12);
        set(haxes,'YLim',[-5,5]);
        set(haxes,'YTick',-5:5:5);
        set(haxes,'YTickLabel',-5:5:5);
    else if strcmp(propertyStr,'Coordination')
            assert(max(abs(score(:,1))) < 1.6);
            assert(max(abs(score(:,2))) < 0.62);
            set(haxes,'XLim',[-1.6,1.6]);
            set(haxes,'XTick',-1.5:1.5:1.5);
            set(haxes,'XTickLabel',-1.5:1.5:1.5);
            set(haxes,'YLim',[-0.62,0.62]);
            set(haxes,'YTick',-0.6:0.3:0.6);
            set(haxes,'YTickLabel',-0.6:0.3:0.6);
        end
    end
end

set(haxes,'FontSize',fontsize);

set(h,'Color','none');

legend('Location','EastOutside');

plot(0,0,'*k','LineWidth',4,'MarkerSize',20);

controlGeneVarFname = [mainDirname 'day/pcaMean_' propertyStr '_legend.eps'];
export_fig(controlGeneVarFname);

legend off;
controlGeneVarFname = [mainDirname 'day/pcaMean_' propertyStr '.eps'];
export_fig(controlGeneVarFname);

hold off;

%% PC Plot
pcPlot(coeff,latent,propertyStr,[mainDirname 'day/']);
close all;

%% PCA: gene - mean(pSup)

[coeffMeanControl,scoreMeanControl,latentMeanControl] = pca(allDiffToMeanControl');
accVarianceMeanControl = cumsum(latentMeanControl)./sum(latentMeanControl);

h = figure;
xlabel('PC 1','FontSize',fontsize);
ylabel('PC 2','FontSize',fontsize);
hold on;

for iGeneDay = 1 : nGeneDay     
    %     plot(scoreMeanControl(geneDayDiff{iGeneDay}.diffToMeanControlInds,1),scoreMeanControl(geneDayDiff{iGeneDay}.diffToMeanControlInds,2),sprintf('%s',markersPerm(ceil(iGeneDay/nColors))),'MarkerEdgeColor','k','MarkerFaceColor',sprintf('%s',colorsPerm(mod(iGeneDay,nColors)+1)),'MarkerSize',10,...
    %         'DisplayName',geneDayDiff{iGeneDay}.geneStr);
    plot(scoreMeanControl(geneDayDiff{iGeneDay}.diffToMeanControlInds,1),scoreMeanControl(geneDayDiff{iGeneDay}.diffToMeanControlInds,2),'o','MarkerEdgeColor',cmap(iGeneDay,:),'LineWidth',2,'MarkerSize',7,'DisplayName',geneDayDiff{iGeneDay}.geneStr);   
end

haxes = get(h,'CurrentAxes');
if strcmp(propertyStr,'Speed')
    assert(max(abs(score(:,1))) < 150);
    assert(max(abs(score(:,2))) < 100);
    set(haxes,'XLim',[-150,150]);
    set(haxes,'XTick',-150:75:150);
    set(haxes,'XTickLabel',-150:75:150);
    set(haxes,'YLim',[-100,100]);
    set(haxes,'YTick',-100:50:100);
    set(haxes,'YTickLabel',-100:50:100);
else if strcmp(propertyStr,'Directionality')
        assert(max(abs(score(:,1))) < 12);
        assert(max(abs(score(:,2))) < 5);
        set(haxes,'XLim',[-12,12]);
        set(haxes,'XTick',-12:6:12);
        set(haxes,'XTickLabel',-12:6:12);
        set(haxes,'YLim',[-5,5]);
        set(haxes,'YTick',-5:5:5);
        set(haxes,'YTickLabel',-5:5:5);
    else if strcmp(propertyStr,'Coordination')
            assert(max(abs(score(:,2))) < 1.6);
            assert(max(abs(score(:,2))) < 0.62);
            set(haxes,'XLim',[-1.6,1.6]);
            set(haxes,'XTick',-1.5:1.5:1.5);
            set(haxes,'XTickLabel',-1.5:1.5:1.5);
            set(haxes,'YLim',[-0.62,0.62]);
            set(haxes,'YTick',-0.6:0.3:0.6);
            set(haxes,'YTickLabel',-0.6:0.3:0.6);
        end
    end
end

set(haxes,'FontSize',fontsize);

set(h,'Color','none');

legend('Location','EastOutside');

plot(0,0,'*k','LineWidth',4,'MarkerSize',20);

controlGeneVarFname = [mainDirname 'day/pcaMeanControl_' propertyStr '_legend.eps'];
export_fig(controlGeneVarFname);

legend off;
controlGeneVarFname = [mainDirname 'day/pcaMeanControl_' propertyStr '.eps'];
export_fig(controlGeneVarFname);

hold off;

%% PCA: mean(gene) - mean(pSup)

[coeffMeanGeneMeanControl,scoreMeanGeneMeanControl,latentMeanGeneMeanControl] = pca(allDiffMeanGeneToMeanControl');
accVarianceMeanGeneMeanControl = cumsum(latentMeanGeneMeanControl)./sum(latentMeanGeneMeanControl);

h = figure;
xlabel('PC 1','FontSize',fontsize);
ylabel('PC 2','FontSize',fontsize);
hold on;

for iGeneDay = 1 : nGeneDay    
    plot(scoreMeanGeneMeanControl(geneDayDiff{iGeneDay}.diffMeanGeneToMeanControlInds,1),scoreMeanGeneMeanControl(geneDayDiff{iGeneDay}.diffMeanGeneToMeanControlInds,2),'o','MarkerEdgeColor',cmap(iGeneDay,:),'LineWidth',2,'MarkerSize',7,'DisplayName',geneDayDiff{iGeneDay}.geneStr);   
end

haxes = get(h,'CurrentAxes');
if strcmp(propertyStr,'Speed')
    assert(max(abs(score(:,1))) < 150);
    assert(max(abs(score(:,2))) < 100);
    set(haxes,'XLim',[-150,150]);
    set(haxes,'XTick',-150:75:150);
    set(haxes,'XTickLabel',-150:75:150);
    set(haxes,'YLim',[-100,100]);
    set(haxes,'YTick',-100:50:100);
    set(haxes,'YTickLabel',-100:50:100);
else if strcmp(propertyStr,'Directionality')
        assert(max(abs(score(:,1))) < 12);
        assert(max(abs(score(:,2))) < 5);
        set(haxes,'XLim',[-12,12]);
        set(haxes,'XTick',-12:6:12);
        set(haxes,'XTickLabel',-12:6:12);
        set(haxes,'YLim',[-5,5]);
        set(haxes,'YTick',-5:5:5);
        set(haxes,'YTickLabel',-5:5:5);
    else if strcmp(propertyStr,'Coordination')
            assert(max(abs(score(:,2))) < 1.6);
            assert(max(abs(score(:,2))) < 0.62);
            set(haxes,'XLim',[-1.6,1.6]);
            set(haxes,'XTick',-1.5:1.5:1.5);
            set(haxes,'XTickLabel',-1.5:1.5:1.5);
            set(haxes,'YLim',[-0.62,0.62]);
            set(haxes,'YTick',-0.6:0.3:0.6);
            set(haxes,'YTickLabel',-0.6:0.3:0.6);
        end
    end
end

set(haxes,'FontSize',fontsize);

set(h,'Color','none');

legend('Location','EastOutside');

plot(0,0,'*k','LineWidth',4,'MarkerSize',20); 

controlGeneVarFname = [mainDirname 'day/pcaMeanGeneMeanControl_' propertyStr '_legend.eps'];
export_fig(controlGeneVarFname);

legend off;
controlGeneVarFname = [mainDirname 'day/pcaMeanGeneMeanControl_' propertyStr '.eps'];
export_fig(controlGeneVarFname);

hold off;
end

%% SVM classification and statistical significance!
function [] = dayClassification(features,geneDayDiff,mainDirname,propertyStr)
addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/algs/libsvm-3.18'),'-begin');
nGeneDay = length(geneDayDiff);

doNormalize = true;

% Gene vs. pSup
for iGeneDay = 1 : nGeneDay    
    rmpath('/apps/MATLAB/R2013a/toolbox/stats/stats');
    curGeneStr = geneDayDiff{iGeneDay}.geneStr;
    geneFeats = features(:,geneDayDiff{iGeneDay}.geneInds);
    controlFeats = features(:,geneDayDiff{iGeneDay}.controlInds);        
    
    nGene = sum(geneDayDiff{iGeneDay}.geneInds);
    nControl = sum(geneDayDiff{iGeneDay}.controlInds);
    feats = [geneFeats controlFeats];
    
    if doNormalize
        for i = 1 : size(feats,1)
            minf = min(feats(i,:));
            maxf = max(feats(i,:));
            feats(i,:) = (feats(i,:)-minf)./(maxf-minf);
        end
    end
    
    labels = [ones(1,nGene) (-1)*ones(1,nControl)];
    %     labels = [ones(1,nGene) zeros(1,nControl)];
    nsamples = length(labels);
    
    
    %     % Similarity
    %     DIST = pdist(feats');
    %     DIST =  squareform(DIST');
    %     figure; imagesc(DIST');
    
    clsWeights = zeros(1,nsamples);
    clsLabels = zeros(1,nsamples);
    for i = 1 : nsamples
        trainFeats = [feats(:,1:i-1) feats(:,i+1:nsamples)];
        trainLabels = [labels(:,1:i-1) labels(:,i+1:nsamples)];
        testFeat = feats(:,i);
        testLabel = labels(i);
        
        %         params.c = 1;
        model = svmtrain(trainLabels', trainFeats', '-c 1');
        %         [model] = cgi_svm_train(trainFeats,trainLabels,params);
        %         [result] = cgi_svm_predict(testFeat,model);
        [predicted_label, accuracy, decision_values] = svmpredict(testLabel,testFeat',model);
        
        %% trying        
        [tmp_predicted_label, tmp_accuracy, tmp_decision_values] = svmpredict(labels',feats',model);
        tmpTrainFeats = [tmp_decision_values(1:i-1); tmp_decision_values(i+1:nsamples)];
        tmpTestFeat = tmp_decision_values(i);
        tmpModel = svmtrain(trainLabels', tmpTrainFeats, '-c 1');
        [predicted_label, accuracy, decision_values] = svmpredict(testLabel,tmpTestFeat',tmpModel);
        
        %%
        
        clsWeights(i) = decision_values;
        clsLabels(i) = predicted_label;
    end
    
    %     % another 2 rounds to make things right
    %     feats = clsWeights;
    %     model = svmtrain(feats', labels', '-c 1');
    %     [predicted_label, accuracy, decision_values] = svmpredict(labels',feats',model);
    %     clsWeights = decision_values;
    %     clsLabels = predicted_label;
    
    
    %     clsWeights = zeros(1,nsamples);
    %     clsLabels = zeros(1,nsamples);
    %     for i = 1 : nsamples
    %         trainFeats = [feats(:,1:i-1) feats(:,i+1:nsamples)];
    %         trainLabels = [labels(:,1:i-1) labels(:,i+1:nsamples)];
    %         testFeat = feats(:,i);
    %         testLabel = labels(i);
    %
    %         %         params.c = 1;
    %         model = svmtrain(trainLabels', trainFeats', '-c 1');
    %         %         [model] = cgi_svm_train(trainFeats,trainLabels,params);
    %         %         [result] = cgi_svm_predict(testFeat,model);
    %         [predicted_label, accuracy, decision_values] = svmpredict(testLabel,testFeat',model);
    %         clsWeights(i) = decision_values;
    %         clsLabels(i) = predicted_label;
    %     end
    %
    %     feats = clsWeights;
    %     clsWeights = zeros(1,nsamples);
    %     clsLabels = zeros(1,nsamples);
    %     for i = 1 : nsamples
    %         trainFeats = [feats(:,1:i-1) feats(:,i+1:nsamples)];
    %         trainLabels = [labels(:,1:i-1) labels(:,i+1:nsamples)];
    %         testFeat = feats(:,i);
    %         testLabel = labels(i);
    %
    %         %         params.c = 1;
    %         model = svmtrain(trainLabels', trainFeats', '-c 1');
    %         %         [model] = cgi_svm_train(trainFeats,trainLabels,params);
    %         %         [result] = cgi_svm_predict(testFeat,model);
    %         [predicted_label, accuracy, decision_values] = svmpredict(testLabel,testFeat',model);
    %         clsWeights(i) = decision_values;
    %         clsLabels(i) = predicted_label;
    %     end
    
    addpath(genpath('/apps/MATLAB/R2013a/toolbox/stats/stats'));
    pval = ranksum(clsWeights(1:nGene),clsWeights(nGene+1:end));
    pval1 = ranksum(clsWeights(nGene+1:end),clsWeights(1:nGene));
    geneDayDiff{iGeneDay}.pval = pval;
    geneDayDiff{iGeneDay}.pval1 = pval1;
    
    % success and ROC
    %     successNum = sum(clsLabels == labels);
    %     nGeneCls = sum(clsLabels == 1 & labels == 1);
    %     nControlCls = sum(clsLabels == -1 & labels == -1);
    %     successPercent = successNum/nsamples;
    %     fprintf(sprintf('%s (%s) = %d / %d, Control = %d / %d (success rate = %g, pval = %g)\n',curGeneStr,propertyStr,nGeneCls,nGene,nControlCls,nControl),successPercent,pval);
    %
    [sortedWeights,inds] = sort(clsWeights);
    sortedAllLabels = (labels(inds) == 1);
    % sortedLabels = (clsLabels(inds) == 1);
    hitRate = zeros(1,nsamples);
    falseRate = zeros(1,nsamples);

    for i = 1 : nsamples
        TH = sortedWeights(i);
        clsLabelsTh = ~(sortedWeights <= TH);
        hits = (clsLabelsTh & sortedAllLabels);
        nhits = sum(hits);
        falses = (clsLabelsTh & ~sortedAllLabels);
        nfalses = sum(falses);
        hitRate(i) = double(nhits)/nGene;
        falseRate(i) = double(nfalses)/nControl;
    end
    
    fontsize = 24;
    
    h = figure;
    plot(falseRate,hitRate,'ok','LineWidth',4,'MarkerSize',15);
    title(sprintf('ROC: %s vs. Control',curGeneStr));
    xlabel('false positive %','FontSize',fontsize);
    ylabel('true positive %','FontSize',fontsize);    
    hold on;
    
    haxes = get(h,'CurrentAxes');
    set(haxes,'XLim',[0,1]);
    set(haxes,'XTick',0:0.5:1);
    set(haxes,'XTickLabel',0:0.5:1);
    set(haxes,'YLim',[0,1]);
    set(haxes,'YTick',0:0.5:1);
    set(haxes,'YTickLabel',0:0.5:1);
    
    plot([0 1],[0 1],'-g','LineWidth',2);
    set(haxes,'FontSize',fontsize);
    
    set(h,'Color','none');
    outFname = [mainDirname 'svm/cls_' curGeneStr '_' propertyStr '.eps'];
    export_fig(outFname);
    hold off;
    close all;
end

fprintf(sprintf('\n ****** Confidence ********\n'));
for iGeneDay = 1 : nGeneDay    
    nGene = sum(geneDayDiff{iGeneDay}.geneInds);
    nControl = sum(geneDayDiff{iGeneDay}.controlInds);
    fprintf(sprintf('%s (gene = %d, control = %d) = %f.4 (%f.4))\n',geneDayDiff{iGeneDay}.geneStr,nGene,nControl,geneDayDiff{iGeneDay}.pval,geneDayDiff{iGeneDay}.pval1));
end

end

%% Clustering
function [] = dayClustering(geneDayDiff,allDiffVectors,allDiffMeanGeneToMeanControl,mainDirname,propertyStr)

nClusters = 5;

fontsize = 24;

nGeneDay = length(geneDayDiff);

labels = [];
genesStr = cell(1,nGeneDay);
for iGeneDay = 1 : nGeneDay 
    genesStr{iGeneDay} = geneDayDiff{iGeneDay}.geneStr;
    labels = [labels iGeneDay*ones(1,size(geneDayDiff{iGeneDay}.diff,2))];
end

clustersInds = kmeans(allDiffVectors',nClusters);

distributionGeneClusters = zeros(nClusters,nGeneDay);
maxGeneCluster = zeros(1,nGeneDay); % encodes for each gene its most probable cluster

for iGeneDay = 1 : nGeneDay 
    indsGene = labels == iGeneDay;
    for iCluster = 1 : nClusters
        distributionGeneClusters(iCluster,iGeneDay) = sum(clustersInds(indsGene) == iCluster);
    end
    distributionGeneClusters(:,iGeneDay) = distributionGeneClusters(:,iGeneDay) ./ sum(distributionGeneClusters(:,iGeneDay));
    [tmp, maxGeneCluster(iGeneDay)] = max(distributionGeneClusters(:,iGeneDay));
end

[sortedClusters,sortedOrder] = sort(maxGeneCluster);
distributionGeneClustersSorted = distributionGeneClusters(:,sortedOrder);

genesStrSorted = genesStr(sortedOrder);

DIST = pdist(distributionGeneClustersSorted');
DIST =  squareform(DIST');
figure; imagesc(DIST');

[coeff,score,latent] = pca(allDiffVectors');

% cmap = colormap(hsv(nClusters));
% fontsize = 24;
% 
% %% PC1 vs. PC2
% h = figure;
% xlabel('PC 1','FontSize',fontsize);
% ylabel('PC 2','FontSize',fontsize);
% hold on;
% 
% for iCluster = 1 : nClusters
%     genesInCluster = find(sortedClusters == iCluster);
%     geneInClusterInds = [];
%     for i = 1 : length(genesInCluster)
%         geneInClusterInds = [geneInClusterInds geneDayDiff{genesInCluster(i)}.diffInds];
%     end
%     plot(score(geneInClusterInds,1),score(geneInClusterInds,2),'o','MarkerEdgeColor',cmap(iCluster,:),'LineWidth',2,'MarkerSize',5);
% end
% 
% haxes = get(h,'CurrentAxes');
% if strcmp(propertyStr,'Speed')
%     assert(max(abs(score(:,1))) < 100);
%     assert(max(abs(score(:,2))) < 50);
%     set(haxes,'XLim',[-100,100]);
%     set(haxes,'XTick',-100:50:100);
%     set(haxes,'XTickLabel',-100:50:100);
%     set(haxes,'YLim',[-50,50]);
%     set(haxes,'YTick',-50:50:50);
%     set(haxes,'YTickLabel',-50:50:50);
% else if strcmp(propertyStr,'Directionality')
%         assert(max(abs(score(:,1))) < 12);
%         assert(max(abs(score(:,2))) < 5);
%         set(haxes,'XLim',[-12,12]);
%         set(haxes,'XTick',-12:6:12);
%         set(haxes,'XTickLabel',-12:6:12);
%         set(haxes,'YLim',[-5,5]);
%         set(haxes,'YTick',-5:5:5);
%         set(haxes,'YTickLabel',-5:5:5);
%     else if strcmp(propertyStr,'Coordination')
%             assert(max(abs(score(:,1))) < 1.6);
%             assert(max(abs(score(:,2))) < 0.62);
%             set(haxes,'XLim',[-1.6,1.6]);
%             set(haxes,'XTick',-1.5:1.5:1.5);
%             set(haxes,'XTickLabel',-1.5:1.5:1.5);
%             set(haxes,'YLim',[-0.62,0.62]);
%             set(haxes,'YTick',-0.6:0.3:0.6);
%             set(haxes,'YTickLabel',-0.6:0.3:0.6);
%         end
%     end
% end
% 
% set(haxes,'FontSize',fontsize);
% 
% set(h,'Color','none');
% 
% plot(0,0,'*k','LineWidth',4,'MarkerSize',20);
% 
% outFname = [mainDirname 'clusters/pca12_' propertyStr '.eps'];
% export_fig(outFname);
% 
% %% PC3 vs. PC4
% h = figure;
% xlabel('PC 3','FontSize',fontsize);
% ylabel('PC 4','FontSize',fontsize);
% hold on;
% 
% for iCluster = 1 : nClusters
%     genesInCluster = find(sortedClusters == iCluster);
%     geneInClusterInds = [];
%     for i = 1 : length(genesInCluster)
%         geneInClusterInds = [geneInClusterInds geneDayDiff{genesInCluster(i)}.diffInds];
%     end
%     plot(score(geneInClusterInds,3),score(geneInClusterInds,4),'o','MarkerEdgeColor',cmap(iCluster,:),'LineWidth',2,'MarkerSize',5);
% end
% 
% haxes = get(h,'CurrentAxes');
% if strcmp(propertyStr,'Speed')
%     assert(max(abs(score(:,3))) < 30);
%     assert(max(abs(score(:,4))) < 25);
%     set(haxes,'XLim',[-30,30]);
%     set(haxes,'XTick',-30:15:30);
%     set(haxes,'XTickLabel',-30:15:30);
%     set(haxes,'YLim',[-20,20]);
%     set(haxes,'YTick',-20:10:20);
%     set(haxes,'YTickLabel',-20:10:20);
% else if strcmp(propertyStr,'Directionality')
%         assert(max(abs(score(:,3))) < 5);
%         assert(max(abs(score(:,4))) < 4);
%         set(haxes,'XLim',[-4,4]);
%         set(haxes,'XTick',-4:2:4);
%         set(haxes,'XTickLabel',-4:2:4);
%         set(haxes,'YLim',[-3,3]);
%         set(haxes,'YTick',-3:1.5:3);
%         set(haxes,'YTickLabel',-3:1.5:3);
%     else if strcmp(propertyStr,'Coordination')
%             assert(max(abs(score(:,3))) < 0.6);
%             assert(max(abs(score(:,4))) < 0.4);
%             set(haxes,'XLim',[-0.4,0.4]);
%             set(haxes,'XTick',-0.4:0.2:0.4);
%             set(haxes,'XTickLabel',-0.4:0.2:0.4);
%             set(haxes,'YLim',[-0.4,0.4]);
%             set(haxes,'YTick',-0.4:0.2:0.4);
%             set(haxes,'YTickLabel',-0.4:0.2:0.4);
%         end
%     end
% end
% 
% set(haxes,'FontSize',fontsize);
% 
% set(h,'Color','none');
% 
% plot(0,0,'*k','LineWidth',4,'MarkerSize',20);
% 
% outFname = [mainDirname 'clusters/pca34_' propertyStr '.eps'];
% export_fig(outFname);

allDiffMeanGeneToMeanControlSorted = allDiffMeanGeneToMeanControl(:,sortedOrder);
TreeMeans = linkage(allDiffMeanGeneToMeanControlSorted','average');
% TreeMeans = linkage(allDiffMeanGeneToMeanControlSorted','average');
DMeans = pdist(allDiffMeanGeneToMeanControlSorted');
LeafOrderMeans = optimalleaforder(TreeMeans,DMeans);
h = figure();
dendrogram(TreeMeans,'Reorder',LeafOrderMeans,'Labels',genesStrSorted);
hold on;
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',fontsize);
childrenAxes = (get(h,'Children'));
grandChildrenAxes = (get(childrenAxes,'Children'));
for i = 1 : length(grandChildrenAxes)
    set(grandChildrenAxes(i),'LineWidth',3);
end
set(h,'Color','none');
position = get(h,'position');
set(h,'position',[position(1:2) round(3*position(3)) position(4)]);
outFname = [mainDirname 'clusters/' propertyStr '_ClusterMean.eps'];
export_fig(outFname);
hold off;
close all;

TreeDistributions = linkage(distributionGeneClustersSorted');
% TreeDistributions = linkage(distributionGeneClustersSorted','average');
DDistributions = pdist(distributionGeneClustersSorted');
LeafOrderDistributions = optimalleaforder(TreeDistributions,DDistributions);
h = figure;
dendrogram(TreeDistributions,'Reorder',LeafOrderDistributions,'Labels',genesStrSorted);
hold on;
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',fontsize);
childrenAxes = (get(h,'Children'));
grandChildrenAxes = (get(childrenAxes,'Children'));
for i = 1 : length(grandChildrenAxes)
    set(grandChildrenAxes(i),'LineWidth',3);
end
set(h,'Color','none');
position = get(h,'position');
set(h,'position',[position(1:2) round(3*position(3)) position(4)]);
outFname = [mainDirname 'clusters/' propertyStr '_ClusterDistribution.eps'];
export_fig(outFname);
hold off;
close all;

genesStrSortedZero = genesStrSorted;
genesStrSortedZero{length(genesStrSortedZero)+1} = 'pSuper';
% with all zeros (control)
allDiffMeanGeneToMeanControlSortedZero = allDiffMeanGeneToMeanControl(:,sortedOrder);
allDiffMeanGeneToMeanControlSortedZero = [allDiffMeanGeneToMeanControlSortedZero,zeros(size(allDiffMeanGeneToMeanControl,1),1)];
TreeMeansZero = linkage(allDiffMeanGeneToMeanControlSortedZero');
% TreeMeansZero = linkage(allDiffMeanGeneToMeanControlSortedZero','average');
DMeansZero = pdist(allDiffMeanGeneToMeanControlSortedZero');
LeafOrderMeansZero = optimalleaforder(TreeMeansZero,DMeansZero);
figure;
dendrogram(TreeMeansZero,'Reorder',LeafOrderMeansZero,'Labels',genesStrSortedZero);
hold on;
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',fontsize);
childrenAxes = (get(h,'Children'));
grandChildrenAxes = (get(childrenAxes,'Children'));
for i = 1 : length(grandChildrenAxes)
    set(grandChildrenAxes(i),'LineWidth',3);
end
set(h,'Color','none');
position = get(h,'position');
set(h,'position',[position(1:2) round(3*position(3)) position(4)]);
outFname = [mainDirname 'clusters/' propertyStr '_ClusterZero.eps'];
export_fig(outFname);
hold off;
close all;
end