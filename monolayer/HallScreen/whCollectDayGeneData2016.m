function [dayGeneData,allDiffMeanGeneToMeanControl,healingRateOut] = whCollectDayGeneData2016(data,healingRate,strLabels,metaData,mainDirname,propertyStr,plotBasicPCAFname)

features = data.features;
kymographs = data.kymographs; % kymographs are already cut by the time and space limitations see whMetaAnalysis2016-->getAllFeatureVectors-->extractFeatures-->getFeatures

dayGeneDataFname = [mainDirname 'dayGeneData_' propertyStr '.mat'];
if exist(dayGeneDataFname,'file')
    load(dayGeneDataFname);
    return;
end

doKymographStd = false;
load(plotBasicPCAFname); % outControls. Usage: score = whGetPCAPrecalc(features,outControls);
if strcmp(propertyStr,'Speed')
    basicParamsPCA = outControls.speed;
else if strcmp(propertyStr,'Directionality')
        basicParamsPCA = outControls.directional;
    else
        if strcmp(propertyStr,'Coordination')
            doKymographStd = false;
            basicParamsPCA = outControls.coordination;
        end
    end
end

if doKymographStd
    kymographsStd = data.kymographsStd; % kymographs are already cut by the time and space limitations
end

nFeats = size(features,1);

% derivatives
derivCoeff.mag = (ones(1,nFeats) .* (1.0/nFeats))';
derivCoeff.timeDerive = [0.7,0.3,-0.3,-0.7,...
    0.7,0.3,-0.3,-0.7,...
    0.7,0.3,-0.3,-0.7]';
derivCoeff.spaceDeriv = [1,1,1,1,...
    0,0,0,0,...
    -1,-1,-1,-1]';

% uniqueGenes = whGetUniqueGeneLabels(strLabels);

N = length(metaData.groupsByDays);
n = length(strLabels);

indsControl = whControlInds(strLabels);

% allDiffVectors = [];
% allDiffToMeanControl = [];
allDiffMeanGeneToMeanControl = [];
% 
healingRateOut = [];
% healingRateControl = [];
% healingRateGene = [];

dayGeneData = {};

nDayGeneSeq = 0;
for day = 1 : N
    dayStr = metaData.groupsByDays{day}.dates;
    daysInds = metaData.groupsByDays{day}.inds;
    indsDay = false(1,n);
    indsDay(daysInds) = true;
    indsControlDay = indsControl & indsDay;
    
    nControl = sum(indsControlDay);
    
    if nControl < 3
        warning('%s: %d < 3 controls',dayStr,nControl);
        continue;
    end    
    
    featsControl = features(:,indsControlDay);
    meanControl = mean(featsControl,2);
    healingRatesControl = healingRate(indsControlDay);
        
    indsGenesDay = indsDay & ~indsControlDay;        
    
    genesStr = metaData.treatment(indsGenesDay);
    
    geneStr = unique(whToPrefix(genesStr));        
    
    nGene = length(geneStr);
    
    for curGene = 1 : nGene
        curGeneStr = geneStr{curGene};
        indsGeneDay = strcmp(whToPrefix(strLabels),curGeneStr);
        indsGeneDay = indsGeneDay & indsDay;        
        
        [seqStr, seqInds] = whGetSequencesStr(strLabels,curGeneStr);
                
        for seq = 1 : length(seqStr)
            indsSeqDay = seqInds{seq};
            indsSeqDay = indsSeqDay & indsGeneDay;
            
            nCurDayGeneSeq = sum(indsSeqDay);
            
            if nCurDayGeneSeq < 3
                warning('%s: %s has %d < 3 repeats',dayStr,seqStr{seq},nControl);
                continue;
            end                        
            
            depletionRate = metaData.KD{find(indsSeqDay,1)};
            seqStrStr = seqStr{seq}; seqStrStr = seqStrStr(2:end-1);
            dayGeneSeqStr = [curGeneStr '_' seqStrStr '_' dayStr ' (' num2str(depletionRate) ')'];            
        
            nDayGeneSeq = nDayGeneSeq + 1;
                        
            featsGene = features(:,indsSeqDay);
            meanGene = mean(featsGene,2);          
            healingRatesGene = healingRate(indsSeqDay);
            
            
            %% Kymographs & printing
            kymographsControl = kymographs(indsControlDay);
            kymographsGene = kymographs(indsSeqDay);
            
            meanKDKymograph = whGetMeanKymograph(kymographsGene);
            meanKControlKymograph = whGetMeanKymograph(kymographsControl);
            diffMeanKymograph = meanKDKymograph - meanKControlKymograph;
            
            if doKymographStd
                kymographsStdControl = kymographsStd(indsControlDay);
                kymographsStdGene = kymographsStd(indsSeqDay);
                
                meanStdKDKymograph = whGetMeanKymograph(kymographsStdGene);
                meanStdKControlKymograph = whGetMeanKymograph(kymographsStdControl);
                diffMeanStdKymograph = meanStdKDKymograph - meanStdKControlKymograph;
            end
            
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
                    if exist(params.fname,'file')
                        continue;
                    end
                    plotKymograph(curKymograph,params);
                end
                
                for curReplicate = 1 : length(kymographsControl)
                    curKymograph = kymographsControl{curReplicate};
                    params.fname = [mainDirname 'dayWellReplicatesKymographs/' dayGeneSeqStr '_' propertyStr '_Ctrl_' num2str(curReplicate) '.eps'];
                    if exist(params.fname,'file')
                        continue;
                    end
                    plotKymograph(curKymograph,params);
                end
                close all;
            end            
            
            %%
            
            diffMeanGeneToMeanControl = meanGene - meanControl;            
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
            %             allDiffVectors = [allDiffVectors diffGene];
            %             allDiffToMeanControl = [allDiffToMeanControl diffGeneToMeanControl];
            allDiffMeanGeneToMeanControl = [allDiffMeanGeneToMeanControl diffMeanGeneToMeanControl];
            %             %         allDiffLabels{nDayGeneSeq} = curGeneStr;
            %
            healingRateOut(nDayGeneSeq) = mean(healingRate(indsSeqDay)) - mean(healingRate(indsControlDay));
            %             healingRateControl(nDayGeneSeq) = mean(healingRate(indsControlDay));
            %             healingRateGene(nDayGeneSeq) = mean(healingRate(indsSeqDay));

            healingRateDiff = mean(healingRate(indsSeqDay)) - mean(healingRate(indsControlDay));
            healingRateControl = mean(healingRate(indsControlDay));
            healingRateGene = mean(healingRate(indsSeqDay));

            dayGeneData{nDayGeneSeq}.pixelSize = metaData.pixelSize{find(indsSeqDay,1)};
            dayGeneData{nDayGeneSeq}.geneStr = curGeneStr;
            dayGeneData{nDayGeneSeq}.dayGeneSeqStr = dayGeneSeqStr;
            dayGeneData{nDayGeneSeq}.SeqStr = seqStrStr;
            dayGeneData{nDayGeneSeq}.KD = depletionRate;
            dayGeneData{nDayGeneSeq}.dayStr = dayStr;
            dayGeneData{nDayGeneSeq}.diff = diffGene; % Matrix of distance-vectors from the day's controls (for each hairpin)
            %             dayGeneData{nDayGeneSeq}.diffInds = (size(allDiffVectors,2)-size(diffGene,2)+1):size(allDiffVectors,2);
            dayGeneData{nDayGeneSeq}.geneFeatures = featsGene;
            dayGeneData{nDayGeneSeq}.controlFeatures = featsControl;
            dayGeneData{nDayGeneSeq}.geneHealingRates = healingRatesGene;
            dayGeneData{nDayGeneSeq}.controlHealingRates = healingRatesControl;
            dayGeneData{nDayGeneSeq}.geneInds = indsSeqDay;
            dayGeneData{nDayGeneSeq}.controlInds = indsControlDay;
            dayGeneData{nDayGeneSeq}.nGeneFeatures = nCurDayGeneSeq;
            dayGeneData{nDayGeneSeq}.nControlFeatures = nControl;
            dayGeneData{nDayGeneSeq}.meanControl = meanControl;
            dayGeneData{nDayGeneSeq}.meanGene = meanGene;
            dayGeneData{nDayGeneSeq}.diffToMeanControl = diffGeneToMeanControl;
            %             dayGeneData{nDayGeneSeq}.diffToMeanControlInds = (size(allDiffToMeanControl,2)-size(diffGeneToMeanControl,2)+1):size(allDiffToMeanControl,2);
            dayGeneData{nDayGeneSeq}.diffMeanGeneToMeanControl = diffMeanGeneToMeanControl;
            dayGeneData{nDayGeneSeq}.diffMeanGeneToMeanControlInds = size(allDiffMeanGeneToMeanControl,2);
            dayGeneData{nDayGeneSeq}.pdist2MeanGeneToMeanControl = pdist2MeanGeneToMeanControl;
            
            % PCA
            dayGeneData{nDayGeneSeq}.featsControlPCA = whGetPCAPrecalc(featsControl,basicParamsPCA);
            dayGeneData{nDayGeneSeq}.featsGenePCA = whGetPCAPrecalc(featsGene,basicParamsPCA);
            dayGeneData{nDayGeneSeq}.meanControlPCA = whGetPCAPrecalc(meanControl,basicParamsPCA);
            dayGeneData{nDayGeneSeq}.meanGenePCA = whGetPCAPrecalc(meanGene,basicParamsPCA);
            
            % magnitude, time and space derivatives
            dayGeneData{nDayGeneSeq}.featsControlDirect = whGetDerivatives(featsControl,basicParamsPCA.meanFeats,basicParamsPCA.stdFeats,derivCoeff);            
            dayGeneData{nDayGeneSeq}.featsGeneDirect = whGetDerivatives(featsGene,basicParamsPCA.meanFeats,basicParamsPCA.stdFeats,derivCoeff);            
            dayGeneData{nDayGeneSeq}.meanControlDirect = whGetDerivatives(meanControl,basicParamsPCA.meanFeats,basicParamsPCA.stdFeats,derivCoeff);            
            dayGeneData{nDayGeneSeq}.meanGeneDirect = whGetDerivatives(meanGene,basicParamsPCA.meanFeats,basicParamsPCA.stdFeats,derivCoeff);
                        
            % Healing rate
            dayGeneData{nDayGeneSeq}.healingRateDiff = healingRateDiff;
            dayGeneData{nDayGeneSeq}.healingRateControl = healingRateControl;
            dayGeneData{nDayGeneSeq}.healingRateGene = healingRateGene;
            
            % Kymograph
            dayGeneData{nDayGeneSeq}.meanKDKymograph = meanKDKymograph;
            dayGeneData{nDayGeneSeq}.meanKControlKymograph = meanKControlKymograph;
            dayGeneData{nDayGeneSeq}.diffKymograph = diffMeanKymograph;
            
            % Std kymographs
            if doKymographStd
                dayGeneData{nDayGeneSeq}.meanStdKDKymograph = meanStdKDKymograph;
                dayGeneData{nDayGeneSeq}.meanStdKControlKymograph = meanStdKControlKymograph;
                dayGeneData{nDayGeneSeq}.diffStdKymograph = diffMeanStdKymograph;
            end
            
            fprintf(sprintf('%s: gene: %d, control: %d\n',dayGeneData{nDayGeneSeq}.geneStr,dayGeneData{nDayGeneSeq}.nGeneFeatures,dayGeneData{nDayGeneSeq}.nControlFeatures))                       
         end
    end    
    close all;        
end

save(dayGeneDataFname,'dayGeneData');

save([mainDirname 'dayGeneData_' propertyStr '_aux.mat'],'allDiffMeanGeneToMeanControl','healingRateOut');

end


