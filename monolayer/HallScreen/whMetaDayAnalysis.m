function [] = whMetaDayAnalysis(allFeatures,healingRate,strLabels,metaData,mainDirname)

warning('off','all');

% addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/algs/libsvm-3.18'),'-begin');
addpath(genpath('/apps/MATLAB/R2013a/toolbox/stats/stats'));

% targetGenesStr = {'TRIO','TUBA','FGD6','ARHGEF3','VAV1','SOS1','ARHGEF11'};
targetGenesStr = {'TRIO','ARHGEF3','ARHGEF11'};

[speedDiffs,healingRateSpeed,geneDayDiffSpeed] = dayAnalysisProperty(allFeatures.speedFeats,targetGenesStr,healingRate,strLabels,metaData,mainDirname,'Speed');
[directionalityDiffs,healingRateDirectionality,geneDayDiffDirectionality] = dayAnalysisProperty(allFeatures.directionalityFeats,targetGenesStr,healingRate,strLabels,metaData,mainDirname,'Directionality');
% % dayAnalysisProperty(allFeatures.strainRateFeats.features,healingRate,strLabels,metaData,mainDirname,'StrainRate');
[coordinationDiffs,healingRateCoordination,geneDayDiffCoordination] = dayAnalysisProperty(allFeatures.coordinationFeats,targetGenesStr,healingRate,strLabels,metaData,mainDirname,'Coordination');

whAssociatePropertyDirection(geneDayDiffSpeed,healingRateSpeed,speedDiffs,directionalityDiffs,coordinationDiffs,targetGenesStr,mainDirname);

close all;
end


% TODO: update by day
function [allDiffMeanGeneToMeanControl,healingRateOut,geneDayDiff] = dayAnalysisProperty(data,targetGenesStr,healingRate,strLabels,metaData,mainDirname,propertyStr)
close all;

features = data.features;
kymographs = data.kymographs;

% uniqueGenes = whGetUniqueGeneLabels(strLabels);

nFeats = size(features,1);
N = length(metaData.groupsByDays);
n = length(strLabels);

strPSup = 'pSuper';
strNT = 'NT';
indsPSup = strcmp(whToPrefix(strLabels),strPSup) | strcmp(whToPrefix(strLabels),strNT); % then filter by days for the gene

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
            
            assert(nCurDayGeneSeq <= 6);                        
            
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
            
            printKymographs = false;
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
            
            diffMeanGeneToMeanControl = meanGene - meanControl;
        
            
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
            
            geneDayDiff{nDayGeneSeq}.meanKDKymograph = meanKDKymograph;
            geneDayDiff{nDayGeneSeq}.meanKControlKymograph = meanKControlKymograph;
            geneDayDiff{nDayGeneSeq}.diffKymograph = diffMeanKymograph;
            fprintf(sprintf('%s: gene: %d, control: %d\n',geneDayDiff{nDayGeneSeq}.geneStr,geneDayDiff{nDayGeneSeq}.nGeneFeatures,geneDayDiff{nDayGeneSeq}.nControlFeatures))
         end
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
% % dayVariance(features,geneDayDiff,indsPSup,mainDirname,propertyStr);
% whVariance(features,geneDayDiff,indsPSup,mainDirname,propertyStr);
% dayAnalysis(allDiffVectors,allDiffToMeanControl,allDiffMeanGeneToMeanControl,geneDayDiff,mainDirname,propertyStr);
% dayClassification(features,geneDayDiff,mainDirname,propertyStr);
% dayClustering(geneDayDiff,allDiffVectors,allDiffMeanGeneToMeanControl,mainDirname,propertyStr);
% dayVisualizeKymographs(geneDayDiff,mainDirname,propertyStr,metaData);
% 
% whVisualizeTargets(geneDayDiff,allDiffMeanGeneToMeanControl,targetGenesStr,mainDirname,propertyStr);
whControlInterdayAssessment(geneDayDiff,mainDirname,propertyStr,metaData,targetGenesStr,healingRateControl,healingRateGene);
end


%%
% function [] = dayVariance(features,geneDayDiff,interDayControlInds,indsPSup,mainDirname,propertyStr)
% function [] = dayVariance(features,geneDayDiff,indsPSup,mainDirname,propertyStr)
% % Calculates variance of 1st PC for each well
% [coeff,score,latent] = pca(features');
% accVariance = cumsum(latent)./sum(latent);
% 
% nGeneDay = length(geneDayDiff);
% 
% varControl = zeros(1,nGeneDay);
% varGene = zeros(1,nGeneDay);
% 
% 
% cmap = colormap(hsv(nGeneDay));
% fontsize = 24;
% 
% h = figure;
% xlabel('Control variance','FontSize',fontsize);
% ylabel('Gene variance','FontSize',fontsize);
% hold on;
% 
% varInterDayPSup = var(score(indsPSup,1));
% % varInterDayControl = var(score(interDayControlInds,1));
% generalVariance = var(score(:,1));
% 
% for iGeneDay = 1 : nGeneDay
%     varControl(iGeneDay) = var(score(geneDayDiff{iGeneDay}.controlInds,1));
%     varGene(iGeneDay) = var(score(geneDayDiff{iGeneDay}.geneInds,1));
% %     plot(varControl(iGeneDay),varGene(iGeneDay),sprintf('%s',markersPerm(ceil(iGeneDay/nColors))),'MarkerEdgeColor','k','MarkerFaceColor',sprintf('%s',colorsPerm(mod(iGeneDay,nColors)+1)),'MarkerSize',10,...
% %         'DisplayName',geneDayDiff{iGeneDay}.geneStr);   
%     plot(varControl(iGeneDay),varGene(iGeneDay),'o','MarkerEdgeColor',cmap(iGeneDay,:),'LineWidth',2,'MarkerSize',7,'DisplayName',[geneDayDiff{iGeneDay}.geneStr '_' geneDayDiff{1}.dayStr]);   
% end
% 
% legend('Location','EastOutside');
% 
% haxes = get(h,'CurrentAxes');
% if strcmp(propertyStr,'Speed')
%     maxVar = 2100;
%     assert(max([varControl varGene]) < maxVar);
%     set(haxes,'XLim',[0,maxVar]);
%     set(haxes,'XTick',0:500:maxVar);
%     set(haxes,'XTickLabel',0:500:maxVar);
%     set(haxes,'YLim',[0,maxVar]);
%     set(haxes,'YTick',0:500:maxVar);
%     set(haxes,'YTickLabel',0:500:maxVar);
% else if strcmp(propertyStr,'Directionality')
%         maxVar = 13;
%         assert(max([varControl varGene]) < maxVar);
%         set(haxes,'XLim',[0,maxVar]);
%         set(haxes,'XTick',0:4:maxVar);
%         set(haxes,'XTickLabel',0:4:maxVar);
%         set(haxes,'YLim',[0,maxVar]);
%         set(haxes,'YTick',0:4:maxVar);
%         set(haxes,'YTickLabel',0:4:maxVar);
%     else if strcmp(propertyStr,'Coordination')
%             maxVar = 0.3;
%             assert(max([varControl varGene]) < maxVar);
%             set(haxes,'XLim',[0,maxVar]);
%             set(haxes,'XTick',0:0.1:maxVar);
%             set(haxes,'XTickLabel',0:0.1:maxVar);
%             set(haxes,'YLim',[0,maxVar]);
%             set(haxes,'YTick',0:0.1:maxVar);
%             set(haxes,'YTickLabel',0:0.1:maxVar);
%         end
%     end
% end
% 
% set(haxes,'FontSize',fontsize);
% 
% set(h,'Color','none');
% 
% 
% plot([0,maxVar],[0,maxVar],'--k','LineWidth',3);
% plot([varInterDayPSup,varInterDayPSup],[0,maxVar],'--g','LineWidth',3);
% % plot([varInterDayControl,varInterDayControl],[0,maxVar],'--b','LineWidth',3);
% plot([0,maxVar],[generalVariance,generalVariance],'--r','LineWidth',3);
% plot([generalVariance,generalVariance],[0,maxVar],'--r','LineWidth',3);
% 
% controlGeneVarFname = [mainDirname 'variance/controlGeneDayVariance' propertyStr '_legend.eps'];
% export_fig(controlGeneVarFname);
% 
% legend off;
% controlGeneVarFname = [mainDirname 'variance/controlGeneDayVariance' propertyStr '.eps'];
% export_fig(controlGeneVarFname);
% 
% hold off;
% 
% %% Daily control PCs
% h = figure;
% xlabel('PC 1','FontSize',fontsize);
% ylabel('PC 2','FontSize',fontsize);
% hold on;
% 
% % plot(score(interDayControlInds,1),score(interDayControlInds,2),'o','MarkerEdgeColor','k','LineWidth',2,'MarkerSize',7,'DisplayName','Control');
% 
% for iGeneDay = 1 : nGeneDay
%     plot(score(geneDayDiff{iGeneDay}.controlInds,1),score(geneDayDiff{iGeneDay}.controlInds,2),'o','MarkerEdgeColor',cmap(iGeneDay,:),'LineWidth',2,'MarkerSize',7,'DisplayName',geneDayDiff{iGeneDay}.dayStr);
% end
% 
% haxes = get(h,'CurrentAxes');
% if strcmp(propertyStr,'Speed')
%     assert(max(score(:,1)) < 100 && min(score(:,1)) > -100);
%     assert(max(score(:,2)) < 50 && min(score(:,2)) > -50);
%     set(haxes,'XLim',[-100,100]);
%     set(haxes,'XTick',-100:50:100);
%     set(haxes,'XTickLabel',-100:50:100);
%     set(haxes,'YLim',[-50,50]);
%     set(haxes,'YTick',-50:25:50);
%     set(haxes,'YTickLabel',-50:25:50);
% else if strcmp(propertyStr,'Directionality')
%         assert(max(score(:,1)) < 8 && min(score(:,1)) > -8);
%         assert(max(score(:,2)) < 7.5 && min(score(:,2)) > -3.2);
%         set(haxes,'XLim',[-8,8]);
%         set(haxes,'XTick',-8:4:8);
%         set(haxes,'XTickLabel',-8:4:8);
%         set(haxes,'YLim',[-3.2,7.5]);
%         set(haxes,'YTick',-3:3:6);
%         set(haxes,'YTickLabel',-3:3:6);
%     else if strcmp(propertyStr,'Coordination')
%             assert(max(score(:,1)) < 1.25 && min(score(:,1)) > -1);
%             assert(max(score(:,2)) < 0.65 && min(score(:,2)) > -0.5);
%             set(haxes,'XLim',[-1,1.25]);
%             set(haxes,'XTick',-1:0.5:1);
%             set(haxes,'XTickLabel',-1:0.5:1);
%             set(haxes,'YLim',[-0.5,0.65]);
%             set(haxes,'YTick',-0.5:0.25:0.5);
%             set(haxes,'YTickLabel',-0.5:0.25:0.5);
%         end
%     end
% end
% 
% set(haxes,'FontSize',fontsize);
% 
% set(h,'Color','none');
% 
% legend('Location','EastOutside');
% 
% controlGeneVarFname = [mainDirname 'variance/pcaPSup_' propertyStr '_legend.eps'];
% export_fig(controlGeneVarFname);
% 
% legend off;
% controlGeneVarFname = [mainDirname 'variance/pcaPSup_' propertyStr '.eps'];
% export_fig(controlGeneVarFname);
% 
% hold off;
% end
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
        assert(max(abs(score(:,1))) < 16);
        assert(max(abs(score(:,2))) < 10);
        set(haxes,'XLim',[-16,16]);
        set(haxes,'XTick',-16:8:16);
        set(haxes,'XTickLabel',-16:8:16);
        set(haxes,'YLim',[-10,10]);
        set(haxes,'YTick',-10:5:10);
        set(haxes,'YTickLabel',-10:5:10);
    else if strcmp(propertyStr,'Coordination')
            assert(max(abs(score(:,1))) < 2);
            assert(max(abs(score(:,2))) < 1);
            set(haxes,'XLim',[-2,2]);
            set(haxes,'XTick',-2:1:2);
            set(haxes,'XTickLabel',-2:1:2);
            set(haxes,'YLim',[-1,1]);
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
            assert(max(abs(score(:,1))) < 16);
            assert(max(abs(score(:,2))) < 10);
            set(haxes,'XLim',[-16,16]);
            set(haxes,'XTick',-16:8:16);
            set(haxes,'XTickLabel',-16:8:16);
            set(haxes,'YLim',[-10,10]);
            set(haxes,'YTick',-10:5:10);
            set(haxes,'YTickLabel',-10:5:10);
        else if strcmp(propertyStr,'Coordination')
                assert(max(abs(score(:,1))) < 2);
                assert(max(abs(score(:,2))) < 1);
                set(haxes,'XLim',[-2,2]);
                set(haxes,'XTick',-2:1:2);
                set(haxes,'XTickLabel',-2:1:2);
                set(haxes,'YLim',[-1,1]);
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
        assert(max(abs(score(:,1))) < 16);
        assert(max(abs(score(:,2))) < 10);
        set(haxes,'XLim',[-16,16]);
        set(haxes,'XTick',-16:8:16);
        set(haxes,'XTickLabel',-16:8:16);
        set(haxes,'YLim',[-10,10]);
        set(haxes,'YTick',-10:5:10);
        set(haxes,'YTickLabel',-10:5:10);
    else if strcmp(propertyStr,'Coordination')
            assert(max(abs(score(:,1))) < 2);
            assert(max(abs(score(:,2))) < 1);
            set(haxes,'XLim',[-2,2]);
            set(haxes,'XTick',-2:1:2);
            set(haxes,'XTickLabel',-2:1:2);
            set(haxes,'YLim',[-1,1]);
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
        assert(max(abs(score(:,1))) < 16);
        assert(max(abs(score(:,2))) < 10);
        set(haxes,'XLim',[-16,16]);
        set(haxes,'XTick',-16:8:16);
        set(haxes,'XTickLabel',-16:8:16);
        set(haxes,'YLim',[-10,10]);
        set(haxes,'YTick',-10:5:10);
        set(haxes,'YTickLabel',-10:5:10);
    else if strcmp(propertyStr,'Coordination')
            assert(max(abs(score(:,1))) < 2);
            assert(max(abs(score(:,2))) < 1);
            set(haxes,'XLim',[-2,2]);
            set(haxes,'XTick',-2:1:2);
            set(haxes,'XTickLabel',-2:1:2);
            set(haxes,'YLim',[-1,1]);
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
        assert(max(abs(score(:,1))) < 16);
        assert(max(abs(score(:,2))) < 10);
        set(haxes,'XLim',[-16,16]);
        set(haxes,'XTick',-16:8:16);
        set(haxes,'XTickLabel',-16:8:16);
        set(haxes,'YLim',[-10,10]);
        set(haxes,'YTick',-10:5:10);
        set(haxes,'YTickLabel',-10:5:10);
    else if strcmp(propertyStr,'Coordination')
            assert(max(abs(score(:,1))) < 2);
            assert(max(abs(score(:,2))) < 1);
            set(haxes,'XLim',[-2,2]);
            set(haxes,'XTick',-2:1:2);
            set(haxes,'XTickLabel',-2:1:2);
            set(haxes,'YLim',[-1,1]);
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

%% SVM classification and statistical significance!
function [] = dayClassification(features,geneDayDiff,mainDirname,propertyStr)
addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/algs/libsvm-3.18'),'-begin');
nGeneDay = length(geneDayDiff);
% 
% doNormalize = false;

[coeff,score,latent] = pca(features');
PC1 = score(:,1); % magnitude
PC2 = score(:,2); % temporal
PC3 = score(:,3); % spatial

loggerFname = [mainDirname 'svm/log_' propertyStr '.txt'];
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
%             assert(max(abs(score(:,2))) < 1);
%             set(haxes,'XLim',[-2,2]);
%             set(haxes,'XTick',-2:1:2);
%             set(haxes,'XTickLabel',-2:1:2);
%             set(haxes,'YLim',[-1,1]);
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

DaviesBouldinIndicesKD0 = [];
DunnIndicesKD0 = [];
SilhouetteCoefficientsKD0 = [];

DaviesBouldinIndicesCdc42Rac1 = [];
DunnIndicesCdc42Rac1 = [];
SilhouetteCoefficientsCdc42Rac1 = [];

DaviesBouldinIndicesRest = [];
DunnIndicesRest = [];
SilhouetteCoefficientsRest = [];

curExp = 0;
for iGeneDay = 1 : nGeneDay   
    if isfield(geneDayDiff{iGeneDay},'DaviesBouldinIndex')
        curExp = curExp + 1;       
        
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
            DaviesBouldinIndicesKD0 = [DaviesBouldinIndicesKD0 DaviesBouldinIndex];
            DunnIndicesKD0 = [DunnIndicesKD0 DunnIndex];
            SilhouetteCoefficientsKD0 = [SilhouetteCoefficientsKD0 SilhouetteCoefficient];            
        else if strcmp(geneDayDiff{iGeneDay}.geneStr,'RAC1') || strcmp(geneDayDiff{iGeneDay}.geneStr,'CDC42') || strcmp(geneDayDiff{iGeneDay}.geneStr,'beta-PIX')
                DaviesBouldinIndicesCdc42Rac1 = [DaviesBouldinIndicesCdc42Rac1 DaviesBouldinIndex];
                DunnIndicesCdc42Rac1 = [DunnIndicesCdc42Rac1 DunnIndex];
                SilhouetteCoefficientsCdc42Rac1 = [SilhouetteCoefficientsCdc42Rac1 SilhouetteCoefficient];
            else
                DaviesBouldinIndicesRest = [DaviesBouldinIndicesRest DaviesBouldinIndex];
                DunnIndicesRest = [DunnIndicesRest DunnIndex];
                SilhouetteCoefficientsRest = [SilhouetteCoefficientsRest SilhouetteCoefficient];
            end
            
        end
    end
end

nZero = length(DaviesBouldinIndicesKD0);
nKnown = length(DaviesBouldinIndicesCdc42Rac1);
nRest = length(DaviesBouldinIndicesRest);

labels = [repmat(' KD = 0%',nZero,1);repmat('Positive',nKnown,1);repmat('KD > 50%',nRest,1)];
dataDaviesBouldin = [DaviesBouldinIndicesKD0';DaviesBouldinIndicesCdc42Rac1';DaviesBouldinIndicesRest'];
dataDunn = [DunnIndicesKD0';DunnIndicesCdc42Rac1';DunnIndicesRest'];
dataSilhouette = [SilhouetteCoefficientsKD0';SilhouetteCoefficientsCdc42Rac1';SilhouetteCoefficientsRest'];

displayBoxPlot(dataDaviesBouldin,labels,'Davies Bouldin Index',[mainDirname 'svm/' propertyStr representativeStr '_DaviesBouldinZero-effectors.eps']);
displayBoxPlot(dataDunn,labels,'Dunn Index',[mainDirname 'svm/' propertyStr representativeStr '_DunnZero-effectors.eps']);
displayBoxPlot(dataSilhouette,labels,'Silhouette Coefficients',[mainDirname 'svm/' propertyStr representativeStr '_SilhouetteZero-effectors.eps']);

%
[sortedDB,sortedIndsDB] = sort(DaviesBouldinIndices);
sortedDayGeneSeqsStrDB = dayGenesSeqStr(sortedIndsDB);

[sortedDunn,sortedIndsDunn] = sort(DunnIndices);
sortedDayGeneSeqsStrDunn = dayGenesSeqStr(sortedIndsDunn);

[sortedSC,sortedIndsSC] = sort(SilhouetteCoefficients);
sortedDayGeneSeqsStrSC = dayGenesSeqStr(sortedIndsSC);


if printAll
    fprintf(logger,sprintf('\n ****** Davies Bouldin ********\n'));
    for iGeneDay = 1 : length(sortedDayGeneSeqsStrDB)
        fprintf(logger,sprintf('%s: %.2f\n',sortedDayGeneSeqsStrDB{iGeneDay},sortedDB(iGeneDay)));
    end
    
    fprintf(logger,sprintf('\n ****** Dunn ********\n'));
    for iGeneDay = 1 : length(sortedDayGeneSeqsStrDB)
        fprintf(logger,sprintf('%s: %.2f\n',sortedDayGeneSeqsStrDunn{iGeneDay},sortedDunn(iGeneDay)));
    end
    
    fprintf(logger,sprintf('\n ****** Silhouette ********\n'));
    for iGeneDay = 1 : length(sortedDayGeneSeqsStrDB)
        fprintf(logger,sprintf('%s: %.2f\n',sortedDayGeneSeqsStrSC{iGeneDay},sortedSC(iGeneDay)));
    end
end


if strcmp(propertyStr,'Speed')
    DBTH = thresholds.speed.DBTH;
    DunnTH = thresholds.speed.DunnTH;
    SilTH = thresholds.speed.SilTH;
    inds1 = (DaviesBouldinIndicesCdc42Rac1 > DBTH & DunnIndicesCdc42Rac1 > DunnTH & SilhouetteCoefficientsCdc42Rac1 > SilTH);
    inds0 = (DaviesBouldinIndicesKD0 > DBTH & DunnIndicesKD0 > DunnTH & SilhouetteCoefficientsKD0 > SilTH);
else if strcmp(propertyStr,'Directionality')
        DBTH = thresholds.directionality.DBTH;
        DunnTH = thresholds.directionality.DunnTH;
        SilTH = thresholds.directionality.SilTH;
        inds1 = (DaviesBouldinIndicesCdc42Rac1 > DBTH & DunnIndicesCdc42Rac1 > DunnTH & SilhouetteCoefficientsCdc42Rac1 > SilTH);
        inds0 = (DaviesBouldinIndicesKD0 > DBTH & DunnIndicesKD0 > DunnTH & SilhouetteCoefficientsKD0 > SilTH);
    else if strcmp(propertyStr,'Coordination')
            DBTH = thresholds.coordination.DBTH;
            DunnTH = thresholds.coordination.DunnTH;
            SilTH = thresholds.coordination.SilTH;
            inds1 = (DaviesBouldinIndicesCdc42Rac1 > DBTH & DunnIndicesCdc42Rac1 > DunnTH & SilhouetteCoefficientsCdc42Rac1 > SilTH);
            inds0 = (DaviesBouldinIndicesKD0 > DBTH & DunnIndicesKD0 > DunnTH & SilhouetteCoefficientsKD0 > SilTH);
        end
    end
end
fprintf(logger,sprintf('\n --------------------------------\n'));
fprintf(logger,sprintf('\n ****** Targets %s ********\n',representativeStr));
fprintf(logger,sprintf('KD == 0: %d/%d\n',sum(inds0),nZero));
fprintf(logger,sprintf('CDC42/RAC1/beta-PIX: %d/%d\n',sum(inds1),nKnown));
fprintf(logger,sprintf('Threshholds: %.2f, %.2f, %.2f\n\n',DBTH,DunnTH,SilTH));

for iGeneDay = 1 : length(DaviesBouldinIndices)   
    if DaviesBouldinIndices(iGeneDay) > DBTH && DunnIndices(iGeneDay) > DunnTH && SilhouetteCoefficients(iGeneDay) > SilTH
        fprintf(logger,sprintf('%s: %.2f, %.2f, %.2f \n',dayGenesSeqStr{iGeneDay},DaviesBouldinIndices(iGeneDay),DunnIndices(iGeneDay),SilhouetteCoefficients(iGeneDay)));
    end
end
end


function [] = whAssociatePropertyDirection(geneDayDiff,healingRate,speedDiffs,directionalityDiffs,coordinationDiffs,targetGenesStr,mainDirname)
n = length(speedDiffs);
assert(n == length(directionalityDiffs) && n == length(coordinationDiffs));

speed5 = prctile(speedDiffs(:),5);
speed95 = prctile(speedDiffs(:),95);
directionality5 = prctile(directionalityDiffs(:),5);
directionality95 = prctile(directionalityDiffs(:),95);
coordination5 = prctile(coordinationDiffs(:),5);
coordination95 = prctile(coordinationDiffs(:),95);

speedDiffs1 = speedDiffs;
directionalityDiffs1 = directionalityDiffs;
coordinationDiffs1 = coordinationDiffs;


% > 0 
speedDiffs1(speedDiffs < speed5) = speed5;
speedDiffs1(speedDiffs > speed95) = speed95;
% speedDiffs = speedDiffs - speed5;
directionalityDiffs1(directionalityDiffs < directionality5) = directionality5;
directionalityDiffs1(directionalityDiffs > directionality95) = directionality95;
% directionalityDiffs = directionalityDiffs - directionality5;
coordinationDiffs1(coordinationDiffs < coordination5) = coordination5;
coordinationDiffs1(coordinationDiffs > coordination95) = coordination95;
% coordinationDiffs = coordinationDiffs - coordination5;

% normalize
speedDiffsNorm = speedDiffs1/ max(abs(speed95),abs(speed5));
directionalityDiffsNorm = directionalityDiffs1/ max(abs(directionality95),abs(directionality5));
coordinationDiffsNorm = coordinationDiffs1/ max(abs(coordination95),abs(coordination5));

allDiffs = [speedDiffsNorm,directionalityDiffsNorm,coordinationDiffsNorm];

[coeff,score,latent] = pca(allDiffs');
accVariance = cumsum(latent)./sum(latent);

% plos PCs
pcPlot(coeff,latent,'CombinedPropoerties',[mainDirname 'combinedProperties/']);
close all;

speedDiffsScores = score(1:n,:);
directionalityDiffsScores = score(n+1:2*n,:);
coordinationDiffsScores = score(2*n+1:end,:);

nGeneDay = length(geneDayDiff);
maxpc1 = max(abs(min(score(:,1))),abs(max(score(:,1))));
maxpc2 = max(abs(min(score(:,2))),abs(max(score(:,2))));
maxpc3 = max(abs(min(score(:,3))),abs(max(score(:,3))));
maxpc4 = max(abs(min(score(:,4))),abs(max(score(:,4))));
maxpc5 = max(abs(min(score(:,5))),abs(max(score(:,5))));

corrSpeedDirectionality = nan(1,nGeneDay);
pvalSpeedDirectionality = nan(1,nGeneDay);

corrSpeedCoordination = nan(1,nGeneDay);
pvalSpeedCoordination = nan(1,nGeneDay);

fontsize = 24;
for iGeneDay = 1 : nGeneDay        
    % PC 1 vs. PC 2
    plotPCsDirectionAcrossProperties(geneDayDiff,speedDiffsScores,directionalityDiffsScores,coordinationDiffsScores,[-maxpc1 maxpc1],[-maxpc2 maxpc2],mainDirname,fontsize,iGeneDay,2);
    % PC 1 vs. PC 3
    plotPCsDirectionAcrossProperties(geneDayDiff,speedDiffsScores,directionalityDiffsScores,coordinationDiffsScores,[-maxpc1 maxpc1],[-maxpc3 maxpc3],mainDirname,fontsize,iGeneDay,3);    
    % PC 1 vs. PC 4
    plotPCsDirectionAcrossProperties(geneDayDiff,speedDiffsScores,directionalityDiffsScores,coordinationDiffsScores,[-maxpc1 maxpc1],[-maxpc4 maxpc4],mainDirname,fontsize,iGeneDay,4);    
    % PC 1 vs. PC 5
    plotPCsDirectionAcrossProperties(geneDayDiff,speedDiffsScores,directionalityDiffsScores,coordinationDiffsScores,[-maxpc1 maxpc1],[-maxpc5 maxpc5],mainDirname,fontsize,iGeneDay,5);    
    
    % Correlations between speed and directionality
    [corrSpeedDirectionality(iGeneDay), pvalSpeedDirectionality(iGeneDay)] = corr(speedDiffs(:,iGeneDay),directionalityDiffs(:,iGeneDay));
    
    % Correlations between speed and coordination
    [corrSpeedCoordination(iGeneDay), pvalSpeedCoordination(iGeneDay)] = corr(speedDiffs(:,iGeneDay),coordinationDiffs(:,iGeneDay));    
end

[negCtrlInds,restInds,targetsInds,posCntrl] = whGetTargetInds(geneDayDiff,targetGenesStr);

nTargets = length(targetGenesStr);
nConditions = nTargets + 3;
cmap = colormap(hsv(nConditions));
markerSize = 7;
LineWidth = 2;

% Rho vs. pval - correlations
h = figure;
xlabel('Rho','FontSize',fontsize);
ylabel('pval','FontSize',fontsize);
hold on;
% plot(corrSpeedDirectionality,pvalSpeedDirectionality,'o','MarkerEdgeColor','k','LineWidth',2,'MarkerSize',7);       
plot(corrSpeedDirectionality(negCtrlInds),pvalSpeedDirectionality(negCtrlInds),'o','MarkerEdgeColor',cmap(1,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','0% KD');   
plot(corrSpeedDirectionality(restInds),pvalSpeedDirectionality(restInds),'o','MarkerEdgeColor',cmap(2,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','> 50% KD');   
for t = 1 : nTargets
    plot(corrSpeedDirectionality(targetsInds{t}),pvalSpeedDirectionality(targetsInds{t}),'o','MarkerEdgeColor',cmap(t+2,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName',targetGenesStr{t});                
end
plot(corrSpeedDirectionality(posCntrl.inds),pvalSpeedDirectionality(posCntrl.inds),'o','MarkerEdgeColor',cmap(nConditions,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','CDC42/RAC1/beta-PIX');

haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[-1,1]);
set(haxes,'FontSize',fontsize);
set(h,'Color','none');
export_fig([mainDirname 'combinedProperties/SpeedDirectionalityCorr.eps']);

h = figure;
xlabel('Rho','FontSize',fontsize);
ylabel('pval','FontSize',fontsize);
hold on;
% plot(corrSpeedCoordination,pvalSpeedCoordination,'o','MarkerEdgeColor','k','LineWidth',2,'MarkerSize',7);       
plot(corrSpeedCoordination(negCtrlInds),pvalSpeedCoordination(negCtrlInds),'o','MarkerEdgeColor',cmap(1,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','0% KD');   
plot(corrSpeedCoordination(restInds),pvalSpeedCoordination(restInds),'o','MarkerEdgeColor',cmap(2,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','> 50% KD');   
for t = 1 : nTargets
    plot(corrSpeedCoordination(targetsInds{t}),pvalSpeedCoordination(targetsInds{t}),'o','MarkerEdgeColor',cmap(t+2,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName',targetGenesStr{t});                
end
plot(corrSpeedCoordination(posCntrl.inds),pvalSpeedCoordination(posCntrl.inds),'o','MarkerEdgeColor',cmap(nConditions,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','CDC42/RAC1/beta-PIX');

haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[-1,1]);
set(haxes,'FontSize',fontsize);
set(h,'Color','none');
export_fig([mainDirname 'combinedProperties/SpeedCoordinationCorr.eps']);

% correlations distributions
[nelements, xcenters] = hist(corrSpeedDirectionality,-0.8:0.2:0.8);
corrSpeedDirectionalityDistribution = nelements ./ sum(nelements);
h = figure;
hold on;
bar(xcenters,corrSpeedDirectionalityDistribution,'k');
xlabel('Pearson Rho','FontSize',22);
ylabel('Percent','FontSize',22);
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[-1,1]);
% set(haxes,'XTick',-1:0.5:1);
% set(haxes,'XTickLabel',-1:0.5:1);
set(haxes,'YLim',[0,0.5]);
% set(haxes,'YTick',0:0.05:0.1);
% set(haxes,'YTickLabel',0:0.05:0.1);
set(haxes,'FontSize',22);
set(h,'Color','none');
hold off;
outFname = [mainDirname 'combinedProperties/SpeedDirectionalityRhoDistributionCorr.eps'];
export_fig(outFname);

[nelements, xcenters] = hist(corrSpeedCoordination,-0.8:0.2:0.8);
corrSpeedCoordinationDistribution = nelements ./ sum(nelements);
h = figure;
hold on;
bar(xcenters,corrSpeedCoordinationDistribution,'k');
xlabel('Pearson Rho','FontSize',22);
ylabel('Percent','FontSize',22);
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[-1,1]);
% set(haxes,'XTick',-1:0.5:1);
% set(haxes,'XTickLabel',-1:0.5:1);
set(haxes,'YLim',[0,0.5]);
% set(haxes,'YTick',0:0.05:0.1);
% set(haxes,'YTickLabel',0:0.05:0.1);
set(haxes,'FontSize',22);
set(h,'Color','none');
hold off;
outFname = [mainDirname 'combinedProperties/SpeedCoordinationRhoDistributionCorr.eps'];
export_fig(outFname);

% Rho vs. wound healing closure rate
h = figure;
xlabel('Rho','FontSize',22);
ylabel('Healing rate (\mum hour{-1})','FontSize',22);
hold all;
% plot(corrSpeedDirectionality,healingRate,'sk','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',7);
plot(corrSpeedDirectionality(negCtrlInds),healingRate(negCtrlInds),'o','MarkerEdgeColor',cmap(1,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','0% KD');   
plot(corrSpeedDirectionality(restInds),healingRate(restInds),'o','MarkerEdgeColor',cmap(2,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','> 50% KD');   
for t = 1 : nTargets
    plot(corrSpeedDirectionality(targetsInds{t}),healingRate(targetsInds{t}),'o','MarkerEdgeColor',cmap(t+2,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName',targetGenesStr{t});                
end
plot(corrSpeedDirectionality(posCntrl.inds),healingRate(posCntrl.inds),'o','MarkerEdgeColor',cmap(nConditions,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','CDC42/RAC1/beta-PIX');
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',22);
set(h,'Color','none');
hold off;
pcaFname = [mainDirname 'combinedProperties/SpeedDirectionalityRhoHealingRate.eps'];
export_fig(pcaFname);
    
h = figure;
xlabel('Rho','FontSize',22);
ylabel('Healing rate (\mum hour{-1})','FontSize',22);
hold all;
% plot(corrSpeedCoordination,healingRate,'sk','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',7);
plot(corrSpeedCoordination(negCtrlInds),healingRate(negCtrlInds),'o','MarkerEdgeColor',cmap(1,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','0% KD');   
plot(corrSpeedCoordination(restInds),healingRate(restInds),'o','MarkerEdgeColor',cmap(2,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','> 50% KD');   
for t = 1 : nTargets
    plot(corrSpeedCoordination(targetsInds{t}),healingRate(targetsInds{t}),'o','MarkerEdgeColor',cmap(t+2,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName',targetGenesStr{t});                
end
plot(corrSpeedCoordination(posCntrl.inds),healingRate(posCntrl.inds),'o','MarkerEdgeColor',cmap(nConditions,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','CDC42/RAC1/beta-PIX');
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',22);
set(h,'Color','none');
hold off;
pcaFname = [mainDirname 'combinedProperties/SpeedCoordinationRhoHealingRate.eps'];
export_fig(pcaFname);

% Serial number vs. Wound closure rate
% Rho vs. wound healing closure rate
h = figure;
xlabel('Experiment ID','FontSize',22);
ylabel('Healing rate (\mum hour{-1})','FontSize',22);
hold all;
nn = 0;
plot((nn+1):(nn+length(negCtrlInds)),healingRate(negCtrlInds),'o','MarkerEdgeColor',cmap(1,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','0% KD');   
nn = length(negCtrlInds);
plot((nn+1):(nn+length(restInds)),healingRate(restInds),'o','MarkerEdgeColor',cmap(2,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','> 50% KD');   
nn = nn+length(restInds);
for t = 1 : nTargets
    plot((nn+1):(nn+length(targetsInds{t})),healingRate(targetsInds{t}),'o','MarkerEdgeColor',cmap(t+2,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName',targetGenesStr{t});                
    nn = nn+length(targetsInds{t});
end
plot((nn+1):(nn+length(posCntrl.inds)),healingRate(posCntrl.inds),'o','MarkerEdgeColor',cmap(nConditions,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','CDC42/RAC1/beta-PIX');
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',22);
set(h,'Color','none');
hold off;
pcaFname = [mainDirname 'combinedProperties/TargetsHealingRate.eps'];
export_fig(pcaFname);


indHealingRate = abs(healingRate) > 20;


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
