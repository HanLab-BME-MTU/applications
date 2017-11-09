%% Delete?

function [allDiffMeanGeneToMeanControl,healingRateOut,geneDayDiff] = whMetaDayAnalysisProperty2016(data,targetGenesStr,healingRate,strLabels,metaData,mainDirname,propertyStr)
close all;

flags.variance = 0;
flags.dayAnalysisPCA = 0;
flags.dayClassification = 0;
flags.dayClustering = 0; % do not activate!
flags.dayVisKymographs = 1;
flags.dayVisTargets = 0;
flags.controlInterdayAssessment = 0;


features = data.features;
kymographs = data.kymographs;

% uniqueGenes = whGetUniqueGeneLabels(strLabels);

dayGeneDaya = whCollectDayGeneData2016();


nFeats = size(features,1);
N = length(metaData.groupsByDays);
n = length(strLabels);

strPSup = 'pSuper';
strNT = 'NT';
strDMSO = 'DMSO';
indsPSup = strcmp(whToPrefix(strLabels),strPSup) | strcmp(whToPrefix(strLabels),strNT) | strcmp(whToPrefix(strLabels),strDMSO); % then filter by days for the gene

allDiffVectors = [];
allDiffToMeanControl = [];
allDiffMeanGeneToMeanControl = [];

healingRateOut = [];
healingRateControl = [];
healingRateGene = [];

geneDayDiff = {};

nDayGeneSeq = 0;
for day = 1 : N
    dayStr = metaData.groupsByDays{day}.dates;
    daysInds = metaData.groupsByDays{day}.inds;
    indsDay = false(1,n);
    indsDay(daysInds) = true;
    indsPSup1 = indsPSup & indsDay;
    
    nControl = sum(indsPSup1);
    
    if nControl < 3
        warning('%s: %d < 3 controls',dayStr,nControl);
        continue;
    end    
    
    featsControl = features(:,indsPSup1);
    meanControl = mean(featsControl,2);
    
    % Now add the gene for that day (make sure there is only 1 of those)
    indsGenes = indsDay & ~indsPSup;        
    
    genesStr = metaData.treatment(indsGenes);
    
    geneStr = unique(whToPrefix(genesStr));        
    
    nGene = length(geneStr);
    
    for curGene = 1 : nGene
        curGeneStr = geneStr{curGene};
        indsGene = strcmp(whToPrefix(strLabels),curGeneStr);
        indsGene = indsGene & indsDay;        
        
        [seqStr, seqInds] = whGetSequencesStr(strLabels,curGeneStr);
        
        
         for seq = 1 : length(seqStr)
            seqIndsDay = seqInds{seq};
            seqIndsDay = seqIndsDay & indsGene;
            
            nCurDayGeneSeq = sum(seqIndsDay);
            
            if nCurDayGeneSeq < 3
                continue;
            end
            
            %             assert(nCurDayGeneSeq <= 6);
            
            depletionRate = metaData.KD{find(seqIndsDay,1)};
            seqStrStr = seqStr{seq}; seqStrStr = seqStrStr(2:end-1);
            dayGeneSeqStr = [curGeneStr '_' seqStrStr '_' dayStr ' (' num2str(depletionRate) ')'];
            
        
            nDayGeneSeq = nDayGeneSeq + 1;
                        
            featsGene = features(:,seqIndsDay);
            meanGene = mean(featsGene,2);                    
            
            
            %% Kymographs
            kymographsControl = kymographs(indsPSup1);
            kymographsGene = kymographs(seqIndsDay);
            
            meanKDKymograph = whGetMeanKymograph(kymographsGene);
            meanKControlKymograph = whGetMeanKymograph(kymographsControl);
            diffMeanKymograph = meanKDKymograph - meanKControlKymograph;
            
            printKymographs = true;
            if printKymographs
                params.timePerFrame = metaData.timePerFrame;
                params.patchSize = 15;
                
                if strcmp(propertyStr,'Speed')
                    params.caxis = [0 60];
                else if strcmp(propertyStr,'Directionality')
                        params.caxis = [0 8];
                    else
                        if strcmp(propertyStr,'Coordination')
                            params.caxis = [0 1];
                        end
                    end
                end
                
                for curReplicate = 1 : length(kymographsGene)
                    curKymograph = kymographsGene{curReplicate};
                    params.fname = [mainDirname 'dayWellReplicatesKymographs/' dayGeneSeqStr '_' propertyStr '_KD_' num2str(curReplicate) '.eps'];
                    plotKymograph(curKymograph,params);
                end
                
                for curReplicate = 1 : length(kymographsControl)
                    curKymograph = kymographsControl{curReplicate};
                    params.fname = [mainDirname 'dayWellReplicatesKymographs/' dayGeneSeqStr '_' propertyStr '_Ctrl_' num2str(curReplicate) '.eps'];
                    plotKymograph(curKymograph,params);
                end
                close all;
            end            
            %%           
            
            diffMeanGeneToMeanControl = meanGeneflags.variance = 0;
flags.dayAnalysisPCA = 0;
flags.dayClassification = 0;
flags.dayClustering = 0; % do not activate!
flags.dayVisKymographs = 1;
flags.dayVisTargets = 0;
flags.controlInterdayAssessment = 0; - meanControl;
            pdist2MeanGeneToMeanControl = pdist2(meanGene',meanControl');
        
            
            % For each gene-sh includes matrix of distance-vectors from day's
            % control
            diffGene = zeros(nFeats,nCurDayGeneSeq*nControl);
            diffGeneToMeanControl = zeros(nFeats,nCurDayGeneSeq);
            for igene = 1 : nCurDayGeneSeq
                diffGeneToMeanControl(:,igene) = featsGene(:,igene) - meanControl;
                for icontrol = 1 : nControl
                    diffGene(:,(igene-1)*nControl+icontrol) = featsGene(:,igene) - featsControl(:,icontrol);
                end
            end
            allDiffVectors = [allDiffVectors diffGene];
            allDiffToMeanControl = [allDiffToMeanControl diffGeneToMeanControl];
            allDiffMeanGeneToMeanControl = [allDiffMeanGeneToMeanControl diffMeanGeneToMeanControl];
            %         allDiffLabels{nDayGeneSeq} = curGeneStr;
            
            healingRateOut(nDayGeneSeq) = mean(healingRate(seqIndsDay)) - mean(healingRate(indsPSup1));
            healingRateControl(nDayGeneSeq) = mean(healingRate(indsPSup1));
            healingRateGene(nDayGeneSeq) = mean(healingRate(seqIndsDay));
            
            geneDayDiff{nDayGeneSeq}.pixelSize = metaData.pixelSize{find(seqIndsDay,1)};
            geneDayDiff{nDayGeneSeq}.geneStr = curGeneStr;
            geneDayDiff{nDayGeneSeq}.dayGeneSeqStr = dayGeneSeqStr;
            geneDayDiff{nDayGeneSeq}.SeqStr = seqStrStr;
            geneDayDiff{nDayGeneSeq}.KD = depletionRate;
            geneDayDiff{nDayGeneSeq}.dayStr = dayStr;
            geneDayDiff{nDayGeneSeq}.diff = diffGene;
            geneDayDiff{nDayGeneSeq}.diffInds = (size(allDiffVectors,2)-size(diffGene,2)+1):size(allDiffVectors,2);
            geneDayDiff{nDayGeneSeq}.geneFeatures = featsGene;
            geneDayDiff{nDayGeneSeq}.controlFeatures = featsControl;
            geneDayDiff{nDayGeneSeq}.geneInds = seqIndsDay;
            geneDayDiff{nDayGeneSeq}.controlInds = indsPSup1;
            geneDayDiff{nDayGeneSeq}.nGeneFeatures = nCurDayGeneSeq;
            geneDayDiff{nDayGeneSeq}.nControlFeatures = nControl;
            geneDayDiff{nDayGeneSeq}.meanControl = meanControl;
            geneDayDiff{nDayGeneSeq}.meanGene = meanGene;
            geneDayDiff{nDayGeneSeq}.diffToMeanControl = diffGeneToMeanControl;
            geneDayDiff{nDayGeneSeq}.diffToMeanControlInds = (size(allDiffToMeanControl,2)-size(diffGeneToMeanControl,2)+1):size(allDiffToMeanControl,2);
            geneDayDiff{nDayGeneSeq}.diffMeanGeneToMeanControl = diffMeanGeneToMeanControl;
            geneDayDiff{nDayGeneSeq}.diffMeanGeneToMeanControlInds = size(allDiffMeanGeneToMeanControl,2);
            geneDayDiff{nDayGeneSeq}.pdist2MeanGeneToMeanControl = pdist2MeanGeneToMeanControl;
            
            geneDayDiff{nDayGeneSeq}.meanKDKymograph = meanKDKymograph;
            geneDayDiff{nDayGeneSeq}.meanKControlKymograph = meanKControlKymograph;
            geneDayDiff{nDayGeneSeq}.diffKymograph = diffMeanKymograph;
            fprintf(sprintf('%s: gene: %d, control: %d\n',geneDayDiff{nDayGeneSeq}.geneStr,geneDayDiff{nDayGeneSeq}.nGeneFeatures,geneDayDiff{nDayGeneSeq}.nControlFeatures))
            
            printSpaceTime = false;
            if printSpaceTime  
                yLabelStr = propertyStr;
                yLimits = nan;
                if strcmp(propertyStr,'Speed')
                    yLimits = [0 60];
                else if strcmp(propertyStr,'Directionality')
                        yLimits = [0 8];                        
                    else
                        if strcmp(propertyStr,'Coordination')
                            yLimits = [0 1];                            
                        end
                    end
                end
                
                % Space
                xSpace = [30,90,150];
                xLabelStr = 'Distance from edge (\mum)';               
                outFname = [mainDirname 'dayGenePlots/' dayGeneSeqStr '_' propertyStr '_Space.eps'];
                titleStr = [strjoin(strsplit(dayGeneSeqStr,'_'),' ') ': spatial ' lower(propertyStr)];
                
                ysControl = [mean(featsControl(1:4,:));mean(featsControl(5:8,:));mean(featsControl(9:12,:))];
                ysKD = [mean(featsGene(1:4,:));mean(featsGene(5:8,:));mean(featsGene(9:12,:))];                
                
                plotTimeOrSpacePlot(xSpace,ysControl,ysKD,xLabelStr,yLabelStr,yLimits,titleStr,outFname);
                
                % Time
                xTime = [28,78,128,178];
                xLabelStr = 'Time (minutes)';                                
                outFname = [mainDirname 'dayGenePlots/' dayGeneSeqStr '_' propertyStr '_Time.eps'];
                titleStr = [strjoin(strsplit(dayGeneSeqStr,'_'),' ') ': temporal ' lower(propertyStr)];
                
                ysControl = [mean(featsControl(1:4:12,:));mean(featsControl(2:4:12,:));mean(featsControl(3:4:12,:));mean(featsControl(4:4:12,:))];
                ysKD = [mean(featsGene(1:4:12,:));mean(featsGene(2:4:12,:));mean(featsGene(3:4:12,:));mean(featsGene(4:4:12,:))];      
                
                plotTimeOrSpacePlot(xTime,ysControl,ysKD,xLabelStr,yLabelStr,yLimits,titleStr,outFname);
                
                close all;
            end
         end
    end    
    close all;        
end

save([mainDirname '/plithotaxisOut/geneDay' propertyStr '.mat'],'geneDayDiff');

%% 2. PCA of (gene KD - control)

[coeff,score,latent] = pca(allDiffVectors');
accVariance = cumsum(latent)./sum(latent);

%% Unused code
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
% % dayVariance(features,geneDayDiff,indsPSup,mainDirname,propertyStr);

%% 3. Assess variance
if flags.variance
    geneDayDiff = whVariance(features,geneDayDiff,indsPSup,mainDirname,propertyStr);
    save([mainDirname '/plithotaxisOut/geneDay' propertyStr '.mat'],'geneDayDiff');
end

%% 4. dayAnalysis: PCA analysis: printing PCs, mean PCA per day, and for each condition KD - Control
if flags.dayAnalysisPCA
    dayAnalysis(allDiffVectors,allDiffToMeanControl,allDiffMeanGeneToMeanControl,geneDayDiff,mainDirname,propertyStr);
end

%% 5. dayClassification - hit detection via cluster seperation measures (in svm directory)
if flags.dayClassification
    dayClassification(features,geneDayDiff,mainDirname,propertyStr);
end

%% 6. dayClustering: clustering GEFs by phenotype (not active)
if flags.dayClustering
    dayClustering(geneDayDiff,allDiffVectors,allDiffMeanGeneToMeanControl,mainDirname,propertyStr);
end

%% 7. dayVisualizeKymographs: mean kymographs of control, KD and the subtraction KD - control (in dayGeneControlKymograph directory)
if flags.dayVisKymographs
    dayVisualizeKymographs(geneDayDiff,mainDirname,propertyStr,metaData);
end

%% 8. whVisualizeTargets: PCA of (KD - control), highlighting targets. Output in targets directory.
if flags.dayVisTargets
    whVisualizeTargets(geneDayDiff,allDiffMeanGeneToMeanControl,targetGenesStr,mainDirname,propertyStr);
end

%% 9. whControlInterdayAssessment: 
if flags.controlInterdayAssessment
    whControlInterdayAssessment(geneDayDiff,mainDirname,propertyStr,metaData,targetGenesStr,healingRateControl,healingRateGene);
end
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
        assert(max(abs(score(:,1))) < 17);
        assert(max(abs(score(:,2))) < 10);
        set(haxes,'XLim',[-17,17]);
        set(haxes,'XTick',-16:8:16);
        set(haxes,'XTickLabel',-16:8:16);
        set(haxes,'YLim',[-10,10]);
        set(haxes,'YTick',-10:5:10);
        set(haxes,'YTickLabel',-10:5:10);
    else if strcmp(propertyStr,'Coordination')
            assert(max(abs(score(:,1))) < 2);
            assert(max(abs(score(:,2))) < 1.01);
            set(haxes,'XLim',[-2,2]);
            set(haxes,'XTick',-2:1:2);
            set(haxes,'XTickLabel',-2:1:2);
            set(haxes,'YLim',[-1,1.01]);
            set(haxes,'YTick',-1:0.5:1);
            set(haxes,'YTickLabel',-1:0.5:1);
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
            assert(max(abs(score(:,1))) < 17);
            assert(max(abs(score(:,2))) < 10);
            set(haxes,'XLim',[-17,17]);
            set(haxes,'XTick',-16:8:16);
            set(haxes,'XTickLabel',-16:8:16);
            set(haxes,'YLim',[-10,10]);
            set(haxes,'YTick',-10:5:10);
            set(haxes,'YTickLabel',-10:5:10);
        else if strcmp(propertyStr,'Coordination')
                assert(max(abs(score(:,1))) < 2);
                assert(max(abs(score(:,2))) < 1.01);
                set(haxes,'XLim',[-2,2]);
                set(haxes,'XTick',-2:1:2);
                set(haxes,'XTickLabel',-2:1:2);
                set(haxes,'YLim',[-1,1.01]);
                set(haxes,'YTick',-1:0.5:1);
                set(haxes,'YTickLabel',-1:0.5:1);
            end
        end
    end
    
    set(haxes,'FontSize',fontsize);
    
    set(h,'Color','none');
    
    legend();
    
    plot(0,0,'*k','LineWidth',4,'MarkerSize',20);
    
    controlGeneVarFname = [mainDirname 'day/gene/' geneDayDiff{iGeneDay}.geneStr '_' geneDayDiff{iGeneDay}.SeqStr '_' geneDayDiff{iGeneDay}.dayStr '_' propertyStr '.eps'];
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
        assert(max(abs(score(:,1))) < 17);
        assert(max(abs(score(:,2))) < 10);
        set(haxes,'XLim',[-17,17]);
        set(haxes,'XTick',-16:8:16);
        set(haxes,'XTickLabel',-16:8:16);
        set(haxes,'YLim',[-10,10]);
        set(haxes,'YTick',-10:5:10);
        set(haxes,'YTickLabel',-10:5:10);
    else if strcmp(propertyStr,'Coordination')
            assert(max(abs(score(:,1))) < 2);
            assert(max(abs(score(:,2))) < 1.01);
            set(haxes,'XLim',[-2,2]);
            set(haxes,'XTick',-2:1:2);
            set(haxes,'XTickLabel',-2:1:2);
            set(haxes,'YLim',[-1,1.01]);
            set(haxes,'YTick',-1:0.5:1);
            set(haxes,'YTickLabel',-1:0.5:1);
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
        assert(max(abs(score(:,1))) < 17);
        assert(max(abs(score(:,2))) < 10);
        set(haxes,'XLim',[-17,17]);
        set(haxes,'XTick',-16:8:16);
        set(haxes,'XTickLabel',-16:8:16);
        set(haxes,'YLim',[-10,10]);
        set(haxes,'YTick',-10:5:10);
        set(haxes,'YTickLabel',-10:5:10);
    else if strcmp(propertyStr,'Coordination')
            assert(max(abs(score(:,1))) < 2);
            assert(max(abs(score(:,2))) < 1.01);
            set(haxes,'XLim',[-2,2]);
            set(haxes,'XTick',-2:1:2);
            set(haxes,'XTickLabel',-2:1:2);
            set(haxes,'YLim',[-1,1.01]);
            set(haxes,'YTick',-1:0.5:1);
            set(haxes,'YTickLabel',-1:0.5:1);
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
        assert(max(abs(score(:,1))) < 17);
        assert(max(abs(score(:,2))) < 10);
        set(haxes,'XLim',[-17,17]);
        set(haxes,'XTick',-16:8:16);
        set(haxes,'XTickLabel',-16:8:16);
        set(haxes,'YLim',[-10,10]);
        set(haxes,'YTick',-10:5:10);
        set(haxes,'YTickLabel',-10:5:10);
    else if strcmp(propertyStr,'Coordination')
            assert(max(abs(score(:,1))) < 2);
            assert(max(abs(score(:,2))) < 1.01);
            set(haxes,'XLim',[-2,2]);
            set(haxes,'XTick',-2:1:2);
            set(haxes,'XTickLabel',-2:1:2);
            set(haxes,'YLim',[-1,1.01]);
            set(haxes,'YTick',-1:0.5:1);
            set(haxes,'YTickLabel',-1:0.5:1);
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

%% Cluster measures for detection
% TODO: use Z-score for ranking!
function [] = dayClassification(features,geneDayDiff,mainDirname,propertyStr)
addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/algs/libsvm-3.18'),'-begin');
nGeneDay = length(geneDayDiff);
% 
% doNormalize = false;

[coeff,score,latent] = pca(features');
PC1 = score(:,1); % magnitude
PC2 = score(:,2); % temporal
PC3 = score(:,3); % spatial

loggerFname = [mainDirname 'svm/log_seperation_' propertyStr '.txt'];
logger = fopen(loggerFname,'w');


% Gene vs. pSup
for iGeneDay = 1 : nGeneDay        
    curGeneStr = geneDayDiff{iGeneDay}.geneStr;
    dayGeneSeqStr = geneDayDiff{iGeneDay}.dayGeneSeqStr;
    
    geneFeats = features(:,geneDayDiff{iGeneDay}.geneInds);
    controlFeats = features(:,geneDayDiff{iGeneDay}.controlInds);        
    
    geneFeatsPC1 = PC1(geneDayDiff{iGeneDay}.geneInds)';
    controlFeatsPC1 = PC1(geneDayDiff{iGeneDay}.controlInds)';
    geneFeatsPC2 = PC2(geneDayDiff{iGeneDay}.geneInds)';
    controlFeatsPC2 = PC2(geneDayDiff{iGeneDay}.controlInds)';
    geneFeatsPC3 = PC3(geneDayDiff{iGeneDay}.geneInds)';
    controlFeatsPC3 = PC3(geneDayDiff{iGeneDay}.controlInds)';
    
    nGene = sum(geneDayDiff{iGeneDay}.geneInds);
    nControl = sum(geneDayDiff{iGeneDay}.controlInds);
    
    if nGene < 3 || nControl < 3
        fprintf(logger,sprintf('Excluding %s (gene = %d, control = %d)\n',geneDayDiff{iGeneDay}.geneStr,nGene,nControl));    
        continue;
    end
    
    [DaviesBouldinIndex,DunnIndex,SilhouetteCoefficient] = whDiscriminationMeasures(geneFeats,controlFeats); 
    geneDayDiff{iGeneDay}.DaviesBouldinIndex = DaviesBouldinIndex;
    geneDayDiff{iGeneDay}.DunnIndex = DunnIndex;
    geneDayDiff{iGeneDay}.SilhouetteCoefficient = SilhouetteCoefficient;
    
    [DaviesBouldinIndexPC1,DunnIndexPC1,SilhouetteCoefficientPC1] = whDiscriminationMeasures(geneFeatsPC1,controlFeatsPC1); 
    geneDayDiff{iGeneDay}.DaviesBouldinIndexPC1 = DaviesBouldinIndexPC1;
    geneDayDiff{iGeneDay}.DunnIndexPC1 = DunnIndexPC1;
    geneDayDiff{iGeneDay}.SilhouetteCoefficientPC1 = SilhouetteCoefficientPC1;
    
    [DaviesBouldinIndexPC2,DunnIndexPC2,SilhouetteCoefficientPC2] = whDiscriminationMeasures(geneFeatsPC2,controlFeatsPC2); 
    geneDayDiff{iGeneDay}.DaviesBouldinIndexPC2 = DaviesBouldinIndexPC2;
    geneDayDiff{iGeneDay}.DunnIndexPC2 = DunnIndexPC2;
    geneDayDiff{iGeneDay}.SilhouetteCoefficientPC2 = SilhouetteCoefficientPC2;
    
    [DaviesBouldinIndexPC3,DunnIndexPC3,SilhouetteCoefficientPC3] = whDiscriminationMeasures(geneFeatsPC3,controlFeatsPC3); 
    geneDayDiff{iGeneDay}.DaviesBouldinIndexPC3 = DaviesBouldinIndexPC3;
    geneDayDiff{iGeneDay}.DunnIndexPC3 = DunnIndexPC3;
    geneDayDiff{iGeneDay}.SilhouetteCoefficientPC3 = SilhouetteCoefficientPC3;
end

%% calculate off-target controls mean and std for these indices (for z-score)


% --------------
printAll = true;
thresholds.speed.DBTH = 1.8;
thresholds.speed.DunnTH = 2.7;
thresholds.speed.SilTH = 0.5;

thresholds.directionality.DBTH = 0.7;
thresholds.directionality.DunnTH = 1;
thresholds.directionality.SilTH = 0.08;

thresholds.coordination.DBTH = 0.7;
thresholds.coordination.DunnTH = 1.2;
thresholds.coordination.SilTH = 0.18;
detectScreenTargets(geneDayDiff,logger,mainDirname,propertyStr,thresholds,'',printAll);

% --------------
printAll = false;
thresholdsPC1.speed.DBTH = 1.8;
thresholdsPC1.speed.DunnTH = 2.7;
thresholdsPC1.speed.SilTH = 0.5;

thresholdsPC1.directionality.DBTH = 0.7;
thresholdsPC1.directionality.DunnTH = 1;
thresholdsPC1.directionality.SilTH = 0.08;

thresholdsPC1.coordination.DBTH = 0.7;
thresholdsPC1.coordination.DunnTH = 1.2;
thresholdsPC1.coordination.SilTH = 0.18;
detectScreenTargets(geneDayDiff,logger,mainDirname,propertyStr,thresholdsPC1,'PC1',printAll);

% --------------
printAll = false;
thresholdsPC2.speed.DBTH = 1.8;
thresholdsPC2.speed.DunnTH = 2.7;
thresholdsPC2.speed.SilTH = 0.5;

thresholdsPC2.directionality.DBTH = 0.7;
thresholdsPC2.directionality.DunnTH = 1;
thresholdsPC2.directionality.SilTH = 0.08;

thresholdsPC2.coordination.DBTH = 0.7;
thresholdsPC2.coordination.DunnTH = 1.2;
thresholdsPC2.coordination.SilTH = 0.18;
detectScreenTargets(geneDayDiff,logger,mainDirname,propertyStr,thresholdsPC2,'PC2',printAll);

% --------------
printAll = false;
thresholdsPC3.speed.DBTH = 1.8;
thresholdsPC3.speed.DunnTH = 2.7;
thresholdsPC3.speed.SilTH = 0.5;

thresholdsPC3.directionality.DBTH = 0.7;
thresholdsPC3.directionality.DunnTH = 1;
thresholdsPC3.directionality.SilTH = 0.08;

thresholdsPC3.coordination.DBTH = 0.7;
thresholdsPC3.coordination.DunnTH = 1.2;
thresholdsPC3.coordination.SilTH = 0.18;
detectScreenTargets(geneDayDiff,logger,mainDirname,propertyStr,thresholdsPC3,'PC3',printAll);

% --------------

fclose(logger);

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
% %          assert(max(abs(score(:,1))) < 16);
%         assert(max(abs(score(:,2))) < 10);
%         set(haxes,'XLim',[-16,16]);
%         set(haxes,'XTick',-16:8:16);
%         set(haxes,'XTickLabel',-16:8:16);
%         set(haxes,'YLim',[-10,10]);
%         set(haxes,'YTick',-10:5:10);
%         set(haxes,'YTickLabel',-10:5:10);
% %     else if strcmp(propertyStr,'Coordination')
%             assert(max(abs(score(:,1))) < 2);
%             assert(max(abs(score(:,2))) < 1.01);
%             set(haxes,'XLim',[-2,2]);
%             set(haxes,'XTick',-2:1:2);
%             set(haxes,'XTickLabel',-2:1:2);
%             set(haxes,'YLim',[-1,1.01]);
%             set(haxes,'YTick',-1:0.5:1);
%             set(haxes,'YTickLabel',-1:0.5:1);
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

%%
function [] = displayBoxPlot(data,labels,ylabelStr,outFname)
colors = 'rbg';
fontsize = 22;
h = figure;
hold on;
boxplot(data,labels,'whisker',0.7193); % 0.7193 --> 90%
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
ylabel(ylabelStr);
% oh=findobj(gca,'tag','Outliers');
% set(oh,'Visible','off');
set(h,'Color','none');
% position = get(h,'position');
% set(h,'position',[position(1:2) round(1.5*position(3:4))]);
hold off;
export_fig(outFname);
end

%%

function [] = detectScreenTargets(geneDayDiff,logger,mainDirname,propertyStr,thresholds,representativeStr,printAll)

nGeneDay = length(geneDayDiff);

DaviesBouldinIndices = [];
DunnIndices = [];
SilhouetteCoefficients = [];
dayGenesSeqStr = {};

controlIntrawellVariance = [];
geneIntrawellVariance = [];
kd0InterwellDist = [];

DaviesBouldinIndicesKD0 = [];
DunnIndicesKD0 = [];
SilhouetteCoefficientsKD0 = [];

DaviesBouldinIndicesCdc42Rac1 = [];
DunnIndicesCdc42Rac1 = [];
SilhouetteCoefficientsCdc42Rac1 = [];

DaviesBouldinIndicesRest = [];
DunnIndicesRest = [];
SilhouetteCoefficientsRest = [];

strKD0 = {}; ikd0 = 0;
strCdc42Rac1 = {}; icdc42rac1 = 0;
strRest = {}; irest = 0;

restInds = [];
kd0Inds = [];
Cdc42Rac1Inds = [];

curExp = 0;
for iGeneDay = 1 : nGeneDay   
    if isfield(geneDayDiff{iGeneDay},'DaviesBouldinIndex')
        curExp = curExp + 1;       
        
        controlIntrawellVariance = [controlIntrawellVariance geneDayDiff{iGeneDay}.varControl];
        geneIntrawellVariance = [geneIntrawellVariance geneDayDiff{iGeneDay}.varGene];
        
        if strcmp(representativeStr,'')
            DaviesBouldinIndex = geneDayDiff{iGeneDay}.DaviesBouldinIndex;
            DunnIndex = geneDayDiff{iGeneDay}.DunnIndex;
            SilhouetteCoefficient = geneDayDiff{iGeneDay}.SilhouetteCoefficient;            
        else if strcmp(representativeStr,'PC1')
                DaviesBouldinIndex = geneDayDiff{iGeneDay}.DaviesBouldinIndexPC1;
                DunnIndex = geneDayDiff{iGeneDay}.DunnIndexPC1;
                SilhouetteCoefficient = geneDayDiff{iGeneDay}.SilhouetteCoefficientPC1;
            else if strcmp(representativeStr,'PC2')
                    DaviesBouldinIndex = geneDayDiff{iGeneDay}.DaviesBouldinIndexPC2;
                    DunnIndex = geneDayDiff{iGeneDay}.DunnIndexPC2;
                    SilhouetteCoefficient = geneDayDiff{iGeneDay}.SilhouetteCoefficientPC2;
                else if strcmp(representativeStr,'PC3')
                        DaviesBouldinIndex = geneDayDiff{iGeneDay}.DaviesBouldinIndexPC3;
                        DunnIndex = geneDayDiff{iGeneDay}.DunnIndexPC3;
                        SilhouetteCoefficient = geneDayDiff{iGeneDay}.SilhouetteCoefficientPC3;
                    end
                end
            end
        end
        
        
        DaviesBouldinIndices = [DaviesBouldinIndices DaviesBouldinIndex];
        DunnIndices = [DunnIndices DunnIndex];
        SilhouetteCoefficients = [SilhouetteCoefficients SilhouetteCoefficient];
        dayGenesSeqStr{curExp} = geneDayDiff{iGeneDay}.dayGeneSeqStr;
        
        
        if geneDayDiff{iGeneDay}.KD == 0
            ikd0 = ikd0 + 1;
            strKD0{ikd0} = dayGenesSeqStr{curExp};
            DaviesBouldinIndicesKD0 = [DaviesBouldinIndicesKD0 DaviesBouldinIndex];
            DunnIndicesKD0 = [DunnIndicesKD0 DunnIndex];
            SilhouetteCoefficientsKD0 = [SilhouetteCoefficientsKD0 SilhouetteCoefficient];                                    
                        
            kd0Inds = [kd0Inds iGeneDay];

            kd0InterwellDist = [kd0InterwellDist geneDayDiff{iGeneDay}.pdist2MeanGeneToMeanControl];
        else if strcmp(geneDayDiff{iGeneDay}.geneStr,'RAC1') || strcmp(geneDayDiff{iGeneDay}.geneStr,'CDC42') || ...
                    strcmp(geneDayDiff{iGeneDay}.geneStr,'beta-PIX') %|| strcmp(geneDayDiff{iGeneDay}.geneStr,'RHOA')
                icdc42rac1 = icdc42rac1 + 1;
                strCdc42Rac1{icdc42rac1} = dayGenesSeqStr{curExp};
                DaviesBouldinIndicesCdc42Rac1 = [DaviesBouldinIndicesCdc42Rac1 DaviesBouldinIndex];
                DunnIndicesCdc42Rac1 = [DunnIndicesCdc42Rac1 DunnIndex];
                SilhouetteCoefficientsCdc42Rac1 = [SilhouetteCoefficientsCdc42Rac1 SilhouetteCoefficient];
                
                Cdc42Rac1Inds = [Cdc42Rac1Inds iGeneDay];
            else
                irest = irest + 1;
                strRest{irest} = dayGenesSeqStr{curExp};
                DaviesBouldinIndicesRest = [DaviesBouldinIndicesRest DaviesBouldinIndex];
                DunnIndicesRest = [DunnIndicesRest DunnIndex];
                SilhouetteCoefficientsRest = [SilhouetteCoefficientsRest SilhouetteCoefficient];
                
                restInds = [restInds iGeneDay];
            end
            
        end
    end
end

nZero = length(DaviesBouldinIndicesKD0);
nKnown = length(DaviesBouldinIndicesCdc42Rac1);
nRest = length(DaviesBouldinIndicesRest);

%% off-target controls for z-score calculations
meanDB = mean(DaviesBouldinIndicesKD0); stdDB = std(DaviesBouldinIndicesKD0);
meanDunn = mean(DunnIndicesKD0); stdDunn = std(DunnIndicesKD0);
meanSilhouette = mean(SilhouetteCoefficientsKD0); stdSilhouette = std(SilhouetteCoefficientsKD0);

zScoreDaviesBouldinKD0 = (DaviesBouldinIndicesKD0 - meanDB) ./ stdDB;
zScoreDunnKD0 = (DunnIndicesKD0 - meanDunn) ./ stdDunn;
zScoreSilhouetteKD0 = (SilhouetteCoefficientsKD0 - meanSilhouette) ./ stdSilhouette;

zScoreDaviesBouldinCdc42Rac1 = (DaviesBouldinIndicesCdc42Rac1 - meanDB) ./ stdDB;
zScoreDunnCdc42Rac1 = (DunnIndicesCdc42Rac1 - meanDunn) ./ stdDunn;
zScoreSilhouetteCdc42Rac1 = (SilhouetteCoefficientsCdc42Rac1 - meanSilhouette) ./ stdSilhouette;

zScoreDaviesBouldinRest = (DaviesBouldinIndicesRest - meanDB) ./ stdDB;
zScoreDunnRest = (DunnIndicesRest - meanDunn) ./ stdDunn;
zScoreSilhouetteRest = (SilhouetteCoefficientsRest - meanSilhouette) ./ stdSilhouette;


meanControlIntrawellVar = mean(controlIntrawellVariance); stdControlIntrawellVar = std(controlIntrawellVariance);
meanGeneIntrawellVar = mean(geneIntrawellVariance); stdGeneIntrawellVar = std(geneIntrawellVariance);
meanInterwellDist = mean(kd0InterwellDist); stdInterwellVar = std(kd0InterwellDist);

medianIntrawellVar = median([controlIntrawellVariance geneIntrawellVariance]);
medianInterwellDist = median(kd0InterwellDist);
%%
labels = [repmat(' KD = 0%',nZero,1);repmat('Positive',nKnown,1);repmat('KD > 50%',nRest,1)];
dataDaviesBouldin = [DaviesBouldinIndicesKD0';DaviesBouldinIndicesCdc42Rac1';DaviesBouldinIndicesRest'];
dataDunn = [DunnIndicesKD0';DunnIndicesCdc42Rac1';DunnIndicesRest'];
dataSilhouette = [SilhouetteCoefficientsKD0';SilhouetteCoefficientsCdc42Rac1';SilhouetteCoefficientsRest'];

displayBoxPlot(dataDaviesBouldin,labels,'Davies Bouldin Index',[mainDirname 'svm/' propertyStr representativeStr '_DaviesBouldin.eps']);
displayBoxPlot(dataDunn,labels,'Dunn Index',[mainDirname 'svm/' propertyStr representativeStr '_Dunn.eps']);
displayBoxPlot(dataSilhouette,labels,'Silhouette Coefficients',[mainDirname 'svm/' propertyStr representativeStr '_Silhouette.eps']);

%% boxplots z-score
zScoreTH = 2;
zScoreSumTH = 5;

zScoreDaviesBouldin = [zScoreDaviesBouldinKD0';zScoreDaviesBouldinCdc42Rac1';zScoreDaviesBouldinRest'];
zScoreDunn = [zScoreDunnKD0';zScoreDunnCdc42Rac1';zScoreDunnRest'];
zScoreSilhouette = [zScoreSilhouetteKD0';zScoreSilhouetteCdc42Rac1';zScoreSilhouetteRest'];
displayBoxPlot(zScoreDaviesBouldin,labels,'Z-score',[mainDirname 'svm/' propertyStr representativeStr '_DaviesBouldinZScore.eps']);
displayBoxPlot(zScoreDunn,labels,'Z-score',[mainDirname 'svm/' propertyStr representativeStr '_DunnZScore.eps']);
displayBoxPlot(zScoreSilhouette,labels,'Z-score',[mainDirname 'svm/' propertyStr representativeStr '_SilhouetteZScore.eps']);


sumZScoreRest = zScoreDaviesBouldinRest + zScoreDunnRest + zScoreSilhouetteRest;
[sortedZScoreRest,sortedIndsZScoreRest] = sort(sumZScoreRest);
sortedStrRest = strRest(sortedIndsZScoreRest);
restIndsSorted = restInds(sortedIndsZScoreRest);

if printAll
    fprintf(logger,sprintf('\n ****** Z-score ********\n'));
    for iGeneDay = 1 : length(sortedZScoreRest)
        curOrigInd = restIndsSorted(iGeneDay);
        curSortedInd = sortedIndsZScoreRest(iGeneDay);
        % THIS IS NOT Z-SCORE OF COURSE, JUST # OF MEDIANS
        varControlZscore = geneDayDiff{curOrigInd}.varControl / medianIntrawellVar;%(geneDayDiff{curOrigInd}.varControl - meanControlIntrawellVar) ./ stdControlIntrawellVar;
        varGeneZscore = geneDayDiff{curOrigInd}.varGene / medianIntrawellVar;%(geneDayDiff{curOrigInd}.varGene - meanGeneIntrawellVar) ./ stdGeneIntrawellVar;
        distInterwell = geneDayDiff{curOrigInd}.pdist2MeanGeneToMeanControl / medianInterwellDist;%(geneDayDiff{curOrigInd}.pdist2MeanGeneToMeanControl - meanInterwellDist) ./ stdInterwellVar;
        fprintf(logger,sprintf('%s: %.1f (%.1f,%.1f,%.1f)\n                    dist: %.1f, variance (cnt,kd): (%.1f,%.1f) \n',...
            sortedStrRest{iGeneDay},sortedZScoreRest(iGeneDay),...
            zScoreDaviesBouldinRest(curSortedInd),zScoreDunnRest(curSortedInd),...
            zScoreSilhouetteRest(curSortedInd),...
            distInterwell,varControlZscore,varGeneZscore));
    end
end

% nPosControl = sum((abs(zScoreDaviesBouldinCdc42Rac1) > zScoreTH & abs(zScoreDunnCdc42Rac1) > zScoreTH & abs(zScoreSilhouetteCdc42Rac1) > zScoreTH));
% nOffTargetControl = sum((abs(zScoreDaviesBouldinKD0) > zScoreTH & abs(zScoreDunnKD0) > zScoreTH & abs(zScoreSilhouetteKD0) > zScoreTH));
% nHits = sum((abs(zScoreDaviesBouldinRest) > zScoreTH & abs(zScoreDunnRest) > zScoreTH & abs(zScoreSilhouetteRest) > zScoreTH));

nPosControl = sum((zScoreDaviesBouldinCdc42Rac1 +  zScoreDunnCdc42Rac1 + zScoreSilhouetteCdc42Rac1) >= zScoreSumTH);
nOffTargetControl = sum((zScoreDaviesBouldinKD0 + zScoreDunnKD0 + zScoreSilhouetteKD0) >= zScoreSumTH);
nHits = sum((zScoreDaviesBouldinRest+ zScoreDunnRest+ zScoreSilhouetteRest) >= zScoreSumTH);

fprintf(logger,sprintf('\n --------------------------------\n'));
fprintf(logger,sprintf('\n ****** Targets %s ********\n',representativeStr));
fprintf(logger,sprintf('Targets: %d/%d\n',nHits,nRest));
fprintf(logger,sprintf('KD == 0: %d/%d\n',nOffTargetControl,nZero));
fprintf(logger,sprintf('CDC42/RAC1/beta-PIX: %d/%d\n',nPosControl,nKnown));
% fprintf(logger,sprintf('Threshholds: %.2f, %.2f, %.2f\n\n',zScoreTH,zScoreTH,zScoreTH));
fprintf(logger,sprintf('Threshholds: sum of 3 z-scores > %.1f\n\n',zScoreSumTH));

targetsInds = find((zScoreDaviesBouldinRest + zScoreDunnRest + zScoreSilhouetteRest) >= zScoreSumTH);
nTargets = length(targetsInds);

sumZScoreHits = sumZScoreRest(targetsInds);
[sortedSumZScoreHits,sortedIndsSumZScoreHits] = sort(sumZScoreHits);

for iTarget = 1 : nTargets
    curSortedInd = sortedIndsSumZScoreHits(iTarget);
    curRestInd = targetsInds(curSortedInd);
    curOriginalInd = restInds(curRestInd);
    if strcmp(representativeStr,'')
        % THIS IS NOT Z-SCORE OF COURSE, JUST # OF MEDIANS
        %         varControlZscore = (geneDayDiff{curOriginalInd}.varControl - meanControlIntrawellVar) ./ stdControlIntrawellVar;
        %         varGeneZscore = (geneDayDiff{curOriginalInd}.varGene - meanGeneIntrawellVar) ./ stdGeneIntrawellVar;
        %         distInterwell = (geneDayDiff{curOriginalInd}.pdist2MeanGeneToMeanControl - meanInterwellDist) ./ stdInterwellVar;
        varControlZscore = geneDayDiff{curOriginalInd}.varControl / medianIntrawellVar;
        varGeneZscore = geneDayDiff{curOriginalInd}.varGene / medianIntrawellVar;
        distInterwell = geneDayDiff{curOriginalInd}.pdist2MeanGeneToMeanControl / medianInterwellDist;

        fprintf(logger,sprintf('%s: %.1f (%.1f, %.1f, %.1f)\n                   dist: %.1f, variance (cnt,kd): (%.1f,%.1f) \n',...
            geneDayDiff{curOriginalInd}.dayGeneSeqStr,...
            sortedSumZScoreHits(iTarget),...
            zScoreDaviesBouldinRest(curRestInd),...
            zScoreDunnRest(curRestInd),...
            zScoreSilhouetteRest(curRestInd),...
            distInterwell,varControlZscore,varGeneZscore));
    else
        fprintf(logger,sprintf('%s: %.1f (%.1f, %.1f, %.1f)\n',...
            geneDayDiff{curOriginalInd}.dayGeneSeqStr,...
            sortedSumZScoreHits(iTarget),...
            zScoreDaviesBouldinRest(curRestInd),...
            zScoreDunnRest(curRestInd),...
            zScoreSilhouetteRest(curRestInd)));
    end
end

end



function [] = plotPCsDirectionAcrossProperties(geneDayDiff,speedDiffsScores,directionalityDiffsScores,coordinationDiffsScores,xlimScore,ylimScore,mainDirname,fontsize,iGeneDay,pcNumber)
h = figure;
hold on;
quiver(0,0,speedDiffsScores(iGeneDay,1),speedDiffsScores(iGeneDay,pcNumber),'-k','LineWidth',6);
quiver(0,0,directionalityDiffsScores(iGeneDay,1),directionalityDiffsScores(iGeneDay,pcNumber),'-r','LineWidth',6);
quiver(0,0,coordinationDiffsScores(iGeneDay,1),coordinationDiffsScores(iGeneDay,pcNumber),'-b','LineWidth',6);
xlabel('PC 1','FontSize',fontsize);
ylabel(['PC ' num2str(pcNumber)],'FontSize',fontsize);
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',xlimScore);
set(haxes,'YLim',ylimScore);
set(haxes,'FontSize',fontsize);
set(h,'Color','none');
hold off;
pcaFname = [mainDirname 'combinedProperties/anglesPC/' geneDayDiff{iGeneDay}.dayGeneSeqStr '_PC1PC' num2str(pcNumber) '.eps'];
export_fig(pcaFname);
close all;
end


function [] = plotTimeOrSpacePlot(x,ysControl,ysKD,xLabelStr,yLabelStr,yLimits,titleStr,outFname)
fontsize = 24;
lineWidth = 2;
markerSize = 8;
h = figure;
hold on;
title(titleStr,'FontSize',fontsize);
plot(x,ysControl,'ok--','MarkerEdgeColor','k','LineWidth',lineWidth,'MarkerSize',markerSize,'DisplayName','Control');
plot(x,ysKD,'or--','MarkerEdgeColor','r','LineWidth',lineWidth,'MarkerSize',markerSize,'DisplayName','KD');
ylim(yLimits);
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',32);
xlabel(xLabelStr,'FontSize',fontsize); ylabel(yLabelStr,'FontSize',fontsize);
set(haxes,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
set(h,'Color','w');
hold off;
export_fig(outFname);
end