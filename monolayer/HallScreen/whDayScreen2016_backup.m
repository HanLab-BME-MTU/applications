function [] = whDayScreen2016(dayGeneData,mainDirname,propertyStr)
% output: seperation measures for the screening (not relevant for
% followups);
nGeneDay = length(dayGeneData);

loggerFname = [mainDirname 'screen/log_' propertyStr '.txt'];
logger = fopen(loggerFname,'w');

loggerAllFname = [mainDirname 'screen/logAll_' propertyStr '.txt'];
loggerAll = fopen(loggerAllFname,'w');

% Gene vs. pSup
for iGeneDay = 1 : nGeneDay        
    %     curGeneStr = dayGeneData{iGeneDay}.geneStr;
    %     dayGeneSeqStr = dayGeneData{iGeneDay}.dayGeneSeqStr;
    %
    %     geneFeats = dayGeneData{iGeneDay}.geneFeatures;%features(:,dayGeneData{iGeneDay}.geneInds);
    %     controlFeats = dayGeneData{iGeneDay}.controlFeatures;%features(:,dayGeneData{iGeneDay}.controlInds);
    
    nGene = dayGeneData{iGeneDay}.nGeneFeatures;%sum(dayGeneData{iGeneDay}.geneInds);
    nControl = dayGeneData{iGeneDay}.nControlFeatures;%sum(dayGeneData{iGeneDay}.controlInds);
    
    if nGene < 3 || nControl < 3
        fprintf(logger,sprintf('Excluding %s (gene = %d, control = %d)\n',dayGeneData{iGeneDay}.geneStr,nGene,nControl));    
        continue;
    end
    
    % PCs1-3
    for ipc = 1 : 3        
        controlFeatsPC = dayGeneData{iGeneDay}.featsControlPCA(:,ipc)';
        geneFeatsPC = dayGeneData{iGeneDay}.featsGenePCA(:,ipc)';                 
        
        [DaviesBouldinIndexPC,DunnIndexPC,SilhouetteCoefficientPC] = whDiscriminationMeasures(geneFeatsPC,controlFeatsPC); 
        measuresField = ['measuresPC' num2str(ipc)];
        dayGeneData{iGeneDay} = setfield(dayGeneData{iGeneDay},measuresField,[DaviesBouldinIndexPC,DunnIndexPC,SilhouetteCoefficientPC]);         %#ok<SFLD>
    end
        
    % Healing rate
    controlHealingRates = dayGeneData{iGeneDay}.controlHealingRates;
    geneHealingRates = dayGeneData{iGeneDay}.geneHealingRates;
    [DaviesBouldinIndexHealingRate,DunnIndexHealingRate,SilhouetteCoefficientHealingRate] = whDiscriminationMeasures(geneHealingRates,controlHealingRates);     
    measuresField = 'measure_healingRate';
    dayGeneData{iGeneDay} = setfield(dayGeneData{iGeneDay},measuresField,[DaviesBouldinIndexHealingRate,DunnIndexHealingRate,SilhouetteCoefficientHealingRate]);         %#ok<SFLD>
    
    % Magnitude, temporal and spatial derivative
    directMeasuresStr = {'mag','timeDerive','spaceDeriv'};
    for ims = 1 : 3
        controlFeatsDirect = dayGeneData{iGeneDay}.featsControlDirect(:,ims)';
        geneFeatsDirect = dayGeneData{iGeneDay}.featsGeneDirect(:,ims)';
        
        [DaviesBouldinIndexDirect,DunnIndexDirect,SilhouetteCoefficientDirect] = whDiscriminationMeasures(geneFeatsDirect,controlFeatsDirect); 
        measuresField = ['measure_' directMeasuresStr{ims}];
        dayGeneData{iGeneDay} = setfield(dayGeneData{iGeneDay},measuresField,[DaviesBouldinIndexDirect,DunnIndexDirect,SilhouetteCoefficientDirect]);         %#ok<SFLD>
    end
    
%     ipc1 = 1;
%     geneFeatsPC1 = dayGeneData{iGeneDay}.featsGenePCA(:,ipc1);%PC1(dayGeneData{iGeneDay}.geneInds)';
%     controlFeatsPC1 = dayGeneData{iGeneDay}.featsControlPCA(:,ipc1);%PC1(dayGeneData{iGeneDay}.controlInds)';
%     ipc2 = 2;
%     geneFeatsPC2 = dayGeneData{iGeneDay}.featsGenePCA(:,ipc2);%PC2(dayGeneData{iGeneDay}.geneInds)';
%     controlFeatsPC2 = dayGeneData{iGeneDay}.featsControlPCA(:,ipc2);%PC2(dayGeneData{iGeneDay}.controlInds)';
%     ipc3 = 3;
%     geneFeatsPC3 = dayGeneData{iGeneDay}.featsGenePCA(:,ipc3);%PC3(dayGeneData{iGeneDay}.geneInds)';
%     controlFeatsPC3 = dayGeneData{iGeneDay}.featsControlPCA(:,ipc3);%PC3(dayGeneData{iGeneDay}.controlInds)';
%     
%     nGene = dayGeneData{iGeneDay}.nGeneFeatures;%sum(dayGeneData{iGeneDay}.geneInds);
%     nControl = dayGeneData{iGeneDay}.nControlFeatures;%sum(dayGeneData{iGeneDay}.controlInds);
%     
%     if nGene < 3 || nControl < 3
%         fprintf(logger,sprintf('Excluding %s (gene = %d, control = %d)\n',dayGeneData{iGeneDay}.geneStr,nGene,nControl));    
%         continue;
%     end
%     
%     %% Includes all 12 features, while we focus on the PCs
%     %     [DaviesBouldinIndex,DunnIndex,SilhouetteCoefficient] = whDiscriminationMeasures(geneFeats,controlFeats);
%     %     dayGeneData{iGeneDay}.DaviesBouldinIndex = DaviesBouldinIndex;
%     %     dayGeneData{iGeneDay}.DunnIndex = DunnIndex;
%     %     dayGeneData{iGeneDay}.SilhouetteCoefficient = SilhouetteCoefficient;
%     
%     [DaviesBouldinIndexPC1,DunnIndexPC1,SilhouetteCoefficientPC1] = whDiscriminationMeasures(geneFeatsPC1,controlFeatsPC1); 
%     dayGeneData{iGeneDay}.measuresPC1 = [DaviesBouldinIndexPC1,DunnIndexPC1,SilhouetteCoefficientPC1];
%     dayGeneData{iGeneDay}.DunnIndexPC1 = DunnIndexPC1;
%     dayGeneData{iGeneDay}.SilhouetteCoefficientPC1 = SilhouetteCoefficientPC1;
%     
%     [DaviesBouldinIndexPC2,DunnIndexPC2,SilhouetteCoefficientPC2] = whDiscriminationMeasures(geneFeatsPC2,controlFeatsPC2); 
%     dayGeneData{iGeneDay}.DaviesBouldinIndexPC2 = DaviesBouldinIndexPC2;
%     dayGeneData{iGeneDay}.DunnIndexPC2 = DunnIndexPC2;
%     dayGeneData{iGeneDay}.SilhouetteCoefficientPC2 = SilhouetteCoefficientPC2;
%     
%     [DaviesBouldinIndexPC3,DunnIndexPC3,SilhouetteCoefficientPC3] = whDiscriminationMeasures(geneFeatsPC3,controlFeatsPC3); 
%     dayGeneData{iGeneDay}.DaviesBouldinIndexPC3 = DaviesBouldinIndexPC3;
%     dayGeneData{iGeneDay}.DunnIndexPC3 = DunnIndexPC3;
%     dayGeneData{iGeneDay}.SilhouetteCoefficientPC3 = SilhouetteCoefficientPC3;
end

save([mainDirname propertyStr '_screen.mat'],'dayGeneData');

%% calculate off-target controls mean and std for these indices (for z-score)


printAll = false;
detectScreenTargets(dayGeneData,logger,loggerAll,mainDirname,propertyStr,printAll);

fclose(logger);
fclose(loggerAll);

end


function [] = detectScreenTargets(dayGeneData,logger,loggerAll,mainDirname,propertyStr,printAll)

nGeneDay = length(dayGeneData);

for ipc = 1 : 7
    if ipc <= 3
        representativeStr = ['PC' num2str(ipc)];
        measuresField = ['measures' representativeStr];
    else if ipc <= 6
            directMeasuresStr = {'mag','timeDerive','spaceDeriv'};
            representativeStr = directMeasuresStr{ipc-3};
            measuresField = ['measure_' representativeStr];
        else if ipc == 7
                representativeStr = 'healingRate';
                measuresField = ['measure_' representativeStr];
            end
        end
    end
    
    
    
    
    % dayGeneData{iGeneDay} = setfield(dayGeneData{iGeneDay},measuresField,[DaviesBouldinIndexPC,DunnIndexPC,SilhouetteCoefficientPC]);         %#ok<SFLD>
    
    DaviesBouldinIndices = [];
    DunnIndices = [];
    SilhouetteCoefficients = [];
    dayGenesSeqStr = {};
    geneStr = {};
    
    %     controlIntrawellVariance = [];
    %     geneIntrawellVariance = [];
    %     kd0InterwellDist = [];
    
    DaviesBouldinIndicesKD0 = [];
    DunnIndicesKD0 = [];
    SilhouetteCoefficientsKD0 = [];
    
    DaviesBouldinIndicesPosControl = [];
    DunnIndicesPosControl = [];
    SilhouetteCoefficientsPosControl = [];
    
    DaviesBouldinIndicesScreen = [];
    DunnIndicesScreen = [];
    SilhouetteCoefficientsScreen = [];
    
    strKD0 = {}; ikd0 = 0;
    strPosControl = {}; iPosControl = 0;
    strScreen = {}; iScreen = 0;
    strScreenGene = {};
    
    screenInds = [];
    kd0Inds = [];
    PosControlInds = [];

    curExp = 0;
    for iGeneDay = 1 : nGeneDay
        if isfield(dayGeneData{iGeneDay},measuresField)
            curExp = curExp + 1;
            
            %             controlIntrawellVariance = [controlIntrawellVariance dayGeneData{iGeneDay}.varControl];
            %             geneIntrawellVariance = [geneIntrawellVariance dayGeneData{iGeneDay}.varGene];
            
            measuresPC = getfield(dayGeneData{iGeneDay},measuresField);
            DaviesBouldinIndexPC = measuresPC(1);
            DunnIndexPC = measuresPC(2);
            SilhouetteCoefficientPC = measuresPC(3);              
                        
            DaviesBouldinIndices = [DaviesBouldinIndices DaviesBouldinIndexPC];
            DunnIndices = [DunnIndices DunnIndexPC];
            SilhouetteCoefficients = [SilhouetteCoefficients SilhouetteCoefficientPC];
            dayGenesSeqStr{curExp} = dayGeneData{iGeneDay}.dayGeneSeqStr;
            geneStr{curExp} = dayGeneData{iGeneDay}.geneStr;
            
            
            %% KD == 0
            if dayGeneData{iGeneDay}.KD == 0
                % Patch for ARHGEF40 that do not have healing rate
                % calculated for!
                if ~isnan(DaviesBouldinIndexPC)
                    ikd0 = ikd0 + 1;
                    strKD0{ikd0} = dayGenesSeqStr{curExp};
                    DaviesBouldinIndicesKD0 = [DaviesBouldinIndicesKD0 DaviesBouldinIndexPC];
                    DunnIndicesKD0 = [DunnIndicesKD0 DunnIndexPC];
                    SilhouetteCoefficientsKD0 = [SilhouetteCoefficientsKD0 SilhouetteCoefficientPC];
                    
                    kd0Inds = [kd0Inds iGeneDay];
                    
                    %                 kd0InterwellDist = [kd0InterwellDist dayGeneData{iGeneDay}.pdist2MeanGeneToMeanControl];
                end
                %% Positive control
            else if strcmp(dayGeneData{iGeneDay}.geneStr,'RAC1') || strcmp(dayGeneData{iGeneDay}.geneStr,'CDC42') || ...
                        strcmp(dayGeneData{iGeneDay}.geneStr,'beta-PIX')% || strcmp(dayGeneData{iGeneDay}.geneStr,'RHOA')
                    iPosControl = iPosControl + 1;
                    strPosControl{iPosControl} = dayGenesSeqStr{curExp};
                    DaviesBouldinIndicesPosControl = [DaviesBouldinIndicesPosControl DaviesBouldinIndexPC];
                    DunnIndicesPosControl = [DunnIndicesPosControl DunnIndexPC];
                    SilhouetteCoefficientsPosControl = [SilhouetteCoefficientsPosControl SilhouetteCoefficientPC];
                    
                    PosControlInds = [PosControlInds iGeneDay];
                else
                    iScreen = iScreen + 1;
                    strScreen{iScreen} = dayGenesSeqStr{curExp};
                    strScreenGene{iScreen} = geneStr{curExp};
                    DaviesBouldinIndicesScreen = [DaviesBouldinIndicesScreen DaviesBouldinIndexPC];
                    DunnIndicesScreen = [DunnIndicesScreen DunnIndexPC];
                    SilhouetteCoefficientsScreen = [SilhouetteCoefficientsScreen SilhouetteCoefficientPC];
                    
                    screenInds = [screenInds iGeneDay];
                end
                
            end
        else
            fprintf(logger,sprintf('%s: not included in screen \n',dayGeneData{iGeneDay}.dayGeneSeqStr));            
        end        
    end
    
    nOffTarget = length(DaviesBouldinIndicesKD0);
    nPosControl = length(DaviesBouldinIndicesPosControl);
    nScreen = length(DaviesBouldinIndicesScreen);
    
    %% off-target controls for z-score calculations
    meanDB = mean(DaviesBouldinIndicesKD0); stdDB = std(DaviesBouldinIndicesKD0);
    meanDunn = mean(DunnIndicesKD0); stdDunn = std(DunnIndicesKD0);
    meanSilhouette = mean(SilhouetteCoefficientsKD0); stdSilhouette = std(SilhouetteCoefficientsKD0);
    
    zScoreDaviesBouldinKD0 = (DaviesBouldinIndicesKD0 - meanDB) ./ stdDB;
    zScoreDunnKD0 = (DunnIndicesKD0 - meanDunn) ./ stdDunn;
    zScoreSilhouetteKD0 = (SilhouetteCoefficientsKD0 - meanSilhouette) ./ stdSilhouette;
    
    zScoreDaviesBouldinPosControl = (DaviesBouldinIndicesPosControl - meanDB) ./ stdDB;
    zScoreDunnPosControl = (DunnIndicesPosControl - meanDunn) ./ stdDunn;
    zScoreSilhouettePosControl = (SilhouetteCoefficientsPosControl - meanSilhouette) ./ stdSilhouette;
    
    zScoreDaviesBouldinScreen = (DaviesBouldinIndicesScreen - meanDB) ./ stdDB;
    zScoreDunnScreen = (DunnIndicesScreen - meanDunn) ./ stdDunn;
    zScoreSilhouetteScreen = (SilhouetteCoefficientsScreen - meanSilhouette) ./ stdSilhouette;

    % meanControlIntrawellVar = mean(controlIntrawellVariance); stdControlIntrawellVar = std(controlIntrawellVariance);
    % meanGeneIntrawellVar = mean(geneIntrawellVariance); stdGeneIntrawellVar = std(geneIntrawellVariance);
    % meanInterwellDist = mean(kd0InterwellDist); stdInterwellVar = std(kd0InterwellDist);
    %
    % medianIntrawellVar = median([controlIntrawellVariance geneIntrawellVariance]);
    % medianInterwellDist = median(kd0InterwellDist);
    %% Plot box plots
    labels = [repmat(' KD = 0%',nOffTarget,1);repmat('Positive',nPosControl,1);repmat('KD > 50%',nScreen,1)];
    dataDaviesBouldin = [DaviesBouldinIndicesKD0';DaviesBouldinIndicesPosControl';DaviesBouldinIndicesScreen'];
    dataDunn = [DunnIndicesKD0';DunnIndicesPosControl';DunnIndicesScreen'];
    dataSilhouette = [SilhouetteCoefficientsKD0';SilhouetteCoefficientsPosControl';SilhouetteCoefficientsScreen'];
        
    %     displayBoxPlot(dataDaviesBouldin,labels,'Davies Bouldin Index',[mainDirname 'screen/' propertyStr '_' representativeStr '_DaviesBouldin.eps']);
    %     displayBoxPlot(dataDunn,labels,'Dunn Index',[mainDirname 'screen/' propertyStr '_' representativeStr '_Dunn.eps']);
    %     displayBoxPlot(dataSilhouette,labels,'Silhouette Coefficients',[mainDirname 'screen/' propertyStr '_' representativeStr '_Silhouette.eps']);

    %% boxplots z-score
    
    zScoreDaviesBouldin = [zScoreDaviesBouldinKD0';zScoreDaviesBouldinPosControl';zScoreDaviesBouldinScreen'];
    zScoreDunn = [zScoreDunnKD0';zScoreDunnPosControl';zScoreDunnScreen'];
    zScoreSilhouette = [zScoreSilhouetteKD0';zScoreSilhouettePosControl';zScoreSilhouetteScreen'];
    
    %     displayBoxPlot(zScoreDaviesBouldin,labels,'Z-score',[mainDirname 'screen/' propertyStr '_' representativeStr '_DaviesBouldinZScore.eps']);
    %     displayBoxPlot(zScoreDunn,labels,'Z-score',[mainDirname 'screen/' propertyStr '_' representativeStr '_DunnZScore.eps']);
    %     displayBoxPlot(zScoreSilhouette,labels,'Z-score',[mainDirname 'screen/' propertyStr '_' representativeStr '_SilhouetteZScore.eps']);

    close all;

    sumZScoreKD0 = zScoreDaviesBouldinKD0 + zScoreDunnKD0 + zScoreSilhouetteKD0;
    sumZScorePosControl = zScoreDaviesBouldinPosControl + zScoreDunnPosControl + zScoreSilhouettePosControl;
    sumZScoreScreen = zScoreDaviesBouldinScreen + zScoreDunnScreen + zScoreSilhouetteScreen;
    sumZScore =  [sumZScoreKD0';sumZScorePosControl';sumZScoreScreen'];
    displayBoxPlot(sumZScore,labels,'Z-score',[mainDirname 'screen/' propertyStr '_' representativeStr '_sumZScore.eps'],[-5,20]);
    
    [sortedZScoreScreen,sortedIndsZScoreScreen] = sort(sumZScoreScreen);
    sortedStrScreen = strScreen(sortedIndsZScoreScreen);
    screenIndsSorted = screenInds(sortedIndsZScoreScreen);

    if printAll
        fprintf(logger,sprintf('\n ****** Z-score ********\n'));
        for iGeneDay = 1 : length(sortedZScoreScreen)
            %             curOrigInd = screenIndsSorted(iGeneDay);
            curSortedInd = sortedIndsZScoreScreen(iGeneDay);
            %             % THIS IS NOT Z-SCORE OF COURSE, JUST # OF MEDIANS
            %             varControlZscore = dayGeneData{curOrigInd}.varControl / medianIntrawellVar;%(dayGeneData{curOrigInd}.varControl - meanControlIntrawellVar) ./ stdControlIntrawellVar;
            %             varGeneZscore = dayGeneData{curOrigInd}.varGene / medianIntrawellVar;%(dayGeneData{curOrigInd}.varGene - meanGeneIntrawellVar) ./ stdGeneIntrawellVar;
            %             distInterwell = dayGeneData{curOrigInd}.pdist2MeanGeneToMeanControl / medianInterwellDist;%(dayGeneData{curOrigInd}.pdist2MeanGeneToMeanControl - meanInterwellDist) ./ stdInterwellVar;
            %             fprintf(logger,sprintf('%s: %.1f (%.1f,%.1f,%.1f)\n                    dist: %.1f, variance (cnt,kd): (%.1f,%.1f) \n',...
            %                 sortedStrScreen{iGeneDay},sortedZScoreScreen(iGeneDay),...
            %                 zScoreDaviesBouldinScreen(curSortedInd),zScoreDunnScreen(curSortedInd),...
            %                 zScoreSilhouetteScreen(curSortedInd),...
            %                 distInterwell,varControlZscore,varGeneZscore));
            fprintf(logger,sprintf('%s: %.1f (%.1f,%.1f,%.1f)\n',...
                sortedStrScreen{iGeneDay},sortedZScoreScreen(iGeneDay),...
                zScoreDaviesBouldinScreen(curSortedInd),zScoreDunnScreen(curSortedInd),...
                zScoreSilhouetteScreen(curSortedInd)));
        end
    end
    
    
    %% Logger all
    fprintf(loggerAll,sprintf('\n --------------------------------\n'));
    fprintf(loggerAll,sprintf('\n ****** %s ********\n',representativeStr));
    for iGeneDay = 1 : length(sumZScoreScreen)
        curDayGeneStr = strScreen{iGeneDay};
        curDayGeneStr_split = strsplit(curDayGeneStr);
        
        curStr = curDayGeneStr_split{1};
        kdEffStr = curDayGeneStr_split{2}; kdEffStr = kdEffStr(2:end-1);
        
        fprintf(loggerAll,sprintf('%s %s %.1f \n',...
            curStr,kdEffStr,sumZScoreScreen(iGeneDay)));
    end
    
    
    %%
    
    % nPosControl = sum((abs(zScoreDaviesBouldinPosControl) > zScoreTH & abs(zScoreDunnPosControl) > zScoreTH & abs(zScoreSilhouettePosControl) > zScoreTH));
    % nOffTargetControl = sum((abs(zScoreDaviesBouldinKD0) > zScoreTH & abs(zScoreDunnKD0) > zScoreTH & abs(zScoreSilhouetteKD0) > zScoreTH));
    % nHits = sum((abs(zScoreDaviesBouldinScreen) > zScoreTH & abs(zScoreDunnScreen) > zScoreTH & abs(zScoreSilhouetteScreen) > zScoreTH));

    
    %% Save posi
    %     sumZScorePosControl = (zScoreDaviesBouldinPosControl +  zScoreDunnPosControl + zScoreSilhouettePosControl);
    %     sumZScoreKD0 = (zScoreDaviesBouldinKD0 + zScoreDunnKD0 + zScoreSilhouetteKD0);
    save([mainDirname 'screen/' propertyStr '_' representativeStr '.mat'],'sumZScoreKD0','sumZScorePosControl','sumZScoreScreen','strScreen');
        
    %% Calculate threshold - just for wound healing rate    
    if ipc == 7                
        nBootstrap = 1000;
        allThresholds = estimateZScoreThresholdDistribution(sumZScorePosControl,sumZScoreKD0,nBootstrap);
        
        
        zRange = 1 : 0.1 : 10;
        nZRange = length(zRange);
        fMeasure = nan(1,nZRange);
        
        for zi = 1 : nZRange
            z = zRange(zi);
            nHitsReal = sum(sumZScorePosControl > z);
            nHitsDetected = nHitsReal + sum(sumZScoreKD0 > z);
            precision = nHitsReal / nHitsDetected;
            recall = nHitsReal / length(sumZScorePosControl);
            fMeasure(zi) = 2*(precision * recall) / (precision + recall);
        end
        
        maxFmeausre = max(fMeasure);
        maxInd = find(fMeasure == maxFmeausre,1,'last'); % last
        disp(zRange(maxInd));
        
        % plot figure
        FPosition = [0 0 300 300];
        APosition = [0.2 0.2 0.75 0.75];
        fontsize = 10;
        h = figure;
        xlabel('Z-score','FontSize',fontsize);
        ylabel('F-Measure','FontSize',fontsize);
        hold on;
            
        plot(zRange,fMeasure,'--','Color','k','LineWidth',2);
        
        haxes = get(h,'CurrentAxes');
        axis(haxes);
        set(h,'Color','none'); % set(h,'Color','w');
        axisHandle= findobj(h,'type','axes');
        set(haxes,'XLim',[0,15]);
        set(haxes,'XTick',0:5:15);
        set(haxes,'XTickLabel',0:5:15);
        set(haxes,'YLim',[0,1]);        
        set(haxes,'YTick',0:0.5:1);
        set(haxes,'YTickLabel',0:0.5:1);
        set(haxes,'FontSize',fontsize);
        set(h,'Position',FPosition,'PaperPositionMode','auto');
        set(axisHandle,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
        set(get(axisHandle,'XLabel'),'FontSize',fontsize); set(get(axisHandle,'YLabel'),'FontSize',fontsize);
        export_fig([mainDirname 'screen/' propertyStr '_' representativeStr '_sumZScoreFmeasure.eps']);
        hold off;
        close all;
        
        
        
        %         [maxFmeausre,maxIndFmeasure] = max(fMeasure);
        %         disp(zRange(maxIndFmeasure)); % it picks the first index, the true value is 5.8
    end
    
    
    %% Detect hits
    % zScoreTH = 2;
    zScoreSumTH = 9.8;
    %     zScoreSumTH = 5.8;

    nPosControlHits = sum((zScoreDaviesBouldinPosControl +  zScoreDunnPosControl + zScoreSilhouettePosControl) >= zScoreSumTH);
    nOffTargetControlHits = sum((zScoreDaviesBouldinKD0 + zScoreDunnKD0 + zScoreSilhouetteKD0) >= zScoreSumTH);
    nHits = sum((zScoreDaviesBouldinScreen+ zScoreDunnScreen+ zScoreSilhouetteScreen) >= zScoreSumTH);

    fprintf(logger,sprintf('\n --------------------------------\n'));
    fprintf(logger,sprintf('\n ****** Targets %s ********\n',representativeStr));
    fprintf(logger,sprintf('Screen: %d/%d\n',nHits,nScreen));
    fprintf(logger,sprintf('Off targets: %d/%d\n',nOffTargetControlHits,nOffTarget));
    fprintf(logger,sprintf('Positive controls: %d/%d\n',nPosControlHits,nPosControl));
    % fprintf(logger,sprintf('Threshholds: %.2f, %.2f, %.2f\n\n',zScoreTH,zScoreTH,zScoreTH));
    fprintf(logger,sprintf('Threshholds: sum of 3 z-scores > %.1f\n\n',zScoreSumTH));

    targetsInds = find((zScoreDaviesBouldinScreen + zScoreDunnScreen + zScoreSilhouetteScreen) >= zScoreSumTH);
    nTargets = length(targetsInds);
    
    sumZScoreHits = sumZScoreScreen(targetsInds);
    [sortedSumZScoreHits,sortedIndsSumZScoreHits] = sort(sumZScoreHits);
    
    for iTarget = 1 : nTargets
        curSortedInd = sortedIndsSumZScoreHits(iTarget);
        curScreenInd = targetsInds(curSortedInd);
        curOriginalInd = screenInds(curScreenInd);
        
        %     if strcmp(representativeStr,'')
        %         % THIS IS NOT Z-SCORE OF COURSE, JUST # OF MEDIANS
        %         %         varControlZscore = (dayGeneData{curOriginalInd}.varControl - meanControlIntrawellVar) ./ stdControlIntrawellVar;
        %         %         varGeneZscore = (dayGeneData{curOriginalInd}.varGene - meanGeneIntrawellVar) ./ stdGeneIntrawellVar;
        %         %         distInterwell = (dayGeneData{curOriginalInd}.pdist2MeanGeneToMeanControl - meanInterwellDist) ./ stdInterwellVar;
        %         varControlZscore = dayGeneData{curOriginalInd}.varControl / medianIntrawellVar;
        %         varGeneZscore = dayGeneData{curOriginalInd}.varGene / medianIntrawellVar;
        %         distInterwell = dayGeneData{curOriginalInd}.pdist2MeanGeneToMeanControl / medianInterwellDist;
        %
        %         fprintf(logger,sprintf('%s: %.1f (%.1f, %.1f, %.1f)\n                   dist: %.1f, variance (cnt,kd): (%.1f,%.1f) \n',...
        %             dayGeneData{curOriginalInd}.dayGeneSeqStr,...
        %             sortedSumZScoreHits(iTarget),...
        %             zScoreDaviesBouldinScreen(curScreenInd),...
        %             zScoreDunnScreen(curScreenInd),...
        %             zScoreSilhouetteScreen(curScreenInd),...
        %             distInterwell,varControlZscore,varGeneZscore));
        %     else
        fprintf(logger,sprintf('%s: %.1f (%.1f, %.1f, %.1f)\n',...
            dayGeneData{curOriginalInd}.dayGeneSeqStr,...
            sortedSumZScoreHits(iTarget),...
            zScoreDaviesBouldinScreen(curScreenInd),...
            zScoreDunnScreen(curScreenInd),...
            zScoreSilhouetteScreen(curScreenInd)));
        %     end
    end
    
    %% plot (previously detected hits)
    notValidatedHitsStr = {'ARHGEF10','TRIO','TUBA','ARHGEF9','DOCK10'};
    [notValidatedHitsInds,tmpInds] = getHitsInds(notValidatedHitsStr,strScreenGene,screenInds,dayGeneData);
    notValidatedHitsInds = [notValidatedHitsInds{:}];
    
    hitsStr = {'SOS1','ARHGEF18','ARHGEF3','ARHGEF11','ARHGEF28'};
    [hitsInds,restInds] = getHitsInds(hitsStr,strScreenGene,screenInds,dayGeneData);
    outFname = [mainDirname 'screen/zScoreHits_' propertyStr '_' representativeStr '.eps'];
    plotScreenZScore(sumZScoreKD0,sumZScorePosControl,sumZScoreScreen,hitsStr,hitsInds,restInds,notValidatedHitsInds,outFname);          
end
end
%% 

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

zRange = 1 : 0.1 : 10;
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
        precision = nHitsReal / nHitsDetected;
        recall = nHitsReal / length(sumZScorePosControl_sample);
        fMeasure(zi) = 2*(precision * recall) / (precision + recall);
    end
    
    maxFmeausre = max(fMeasure);
    maxInd = find(fMeasure == maxFmeausre,1,'last'); % last
    allThresholds(iBoot) = zRange(maxInd);
end
end