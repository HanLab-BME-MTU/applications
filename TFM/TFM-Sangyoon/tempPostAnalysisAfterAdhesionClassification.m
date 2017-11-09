function []=tempPostAnalysisAfterAdhesionClassification(pathForColocalization)
    load([pathForColocalization filesep 'data' filesep 'allDataAfterClassification.mat'])
    %% average time from t_init to t_peak
    firstIncreseTimeForce_both = arrayfun(@(x) x.firstIncreseTimeForce, tracksNA(curIndices(bothTiTpIdx)));
    peakTimeForceAll_both = arrayfun(@(x) x.forcePeakFrame*tInterval, tracksNA(curIndices(bothTiTpIdx)));
    timeToPeakForceFT = peakTimeForceAll_both-firstIncreseTimeForce_both;
%     figure, histogram(timeToPeakForceFT)
%     figure, plot(firstIncreseTimeForce_both,peakTimeForceAll_both,'.')
    median(timeToPeakForceFT)

    firstIncreseTimeInt_both = arrayfun(@(x) x.firstIncreseTimeInt, tracksNA(curIndices(bothTiTpIdx)));
    peakTimeIntAll_both = arrayfun(@(x) x.intenPeakFrame*tInterval, tracksNA(curIndices(bothTiTpIdx)));
    timeToPeakIntFT = peakTimeIntAll_both-firstIncreseTimeInt_both;
%     figure, histogram(timeToPeakIntFT)
%     figure, plot(firstIncreseTimeInt_both,peakTimeIntAll_both,'.')
    median(timeToPeakIntFT)

    save([pathForColocalization filesep 'data' filesep 'timeToPeaks.mat'],'timeToPeakForceFT','timeToPeakIntFT')
    

    %% assembly rate and force growth rate in curIndices
    forceSlopeG1 =arrayfun(@(x) (x.forceSlope),tracksNA(curIndices));
    earlyAmpSlopeG1 =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(curIndices));
   
    save([pathForColocalization filesep 'data' filesep 'assemblyRateForceSlopes.mat'],'forceSlopeG1','earlyAmpSlopeG1')
    
    %% non force transmitting adhesion statistics
    % average fluorescence intensity
    avgFluoInten_nonTransmitting = arrayfun(@(x) mean(x.ampTotal(x.presence)),tracksNA(curIndices(nonTransmittingIdx)));
    mean(avgFluoInten_nonTransmitting)
    avgFluoInten_forceTransmitting = arrayfun(@(x) mean(x.ampTotal(x.presence)),tracksNA(curIndices(firstIncreseTimeIntAgainstForceAllIdxIDs)));
    mean(avgFluoInten_forceTransmitting)
    [hFI,pFI]=ttest2(avgFluoInten_nonTransmitting,avgFluoInten_forceTransmitting)
    % only for pre-detection period
    avgPreDetecFluoInten_nonTransmitting = arrayfun(@(x) mean(x.ampTotal(x.startingFrameExtraExtra:x.startingFrameExtra)),tracksNA(curIndices(nonTransmittingIdx)));
    mean(avgPreDetecFluoInten_nonTransmitting)
    avgPreDetecFluoInten_forceTransmitting = arrayfun(@(x) mean(x.ampTotal(x.startingFrameExtraExtra:x.startingFrameExtra)),tracksNA(curIndices(firstIncreseTimeIntAgainstForceAllIdxIDs)));
    mean(avgPreDetecFluoInten_forceTransmitting)
    [hPreDetecFI,pPreDetecFI]=ttest2(avgPreDetecFluoInten_nonTransmitting,avgPreDetecFluoInten_forceTransmitting)
    % earlyAmpSlope
    avgEarlyAmpSlope_nonTransmitting = arrayfun(@(x) x.earlyAmpSlope,tracksNA(curIndices(nonTransmittingIdx)));
    mean(avgEarlyAmpSlope_nonTransmitting)
    avgEarlyAmpSlope_forceTransmitting = arrayfun(@(x) x.earlyAmpSlope,tracksNA(curIndices(firstIncreseTimeIntAgainstForceAllIdxIDs)));
    mean(avgEarlyAmpSlope_forceTransmitting)
    [hEarlyAmpSlope,pEarlyAmpSlope]=ttest2(avgEarlyAmpSlope_nonTransmitting,avgEarlyAmpSlope_forceTransmitting)
    % forceMag for pre-detection period
    avgPreDetecForceMag_nonTransmitting = arrayfun(@(x) mean(x.forceMag(x.startingFrameExtraExtra:x.startingFrameExtra)),tracksNA(curIndices(nonTransmittingIdx)));
    mean(avgPreDetecForceMag_nonTransmitting)
    avgPreDetecForceMag_forceTransmitting = arrayfun(@(x) mean(x.forceMag(x.startingFrameExtraExtra:x.startingFrameExtra)),tracksNA(curIndices(firstIncreseTimeIntAgainstForceAllIdxIDs)));
    mean(avgPreDetecForceMag_forceTransmitting)
    [hPreDetecForce,pPreDetecForce]=ttest2(avgPreDetecForceMag_nonTransmitting,avgPreDetecForceMag_forceTransmitting)
    %% Save these parameters
    save([pathForColocalization filesep 'data' filesep 'nonTransmittingVsForceTransmitting.mat'],'avgFluoInten_nonTransmitting',...
        'avgFluoInten_forceTransmitting','hFI','pFI','avgPreDetecFluoInten_nonTransmitting','avgPreDetecFluoInten_forceTransmitting',...
        'hPreDetecFI','pPreDetecFI','avgEarlyAmpSlope_nonTransmitting','avgEarlyAmpSlope_forceTransmitting',...
        'hEarlyAmpSlope','pEarlyAmpSlope','avgPreDetecForceMag_nonTransmitting','avgPreDetecForceMag_forceTransmitting',...
        'hPreDetecForce','pPreDetecForce','-v7.3')
    %% Look at feature difference per each group
    distToEdge{1} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup1filtered));
    distToEdge{2} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup2filtered));
    distToEdge{3} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup3filtered));
    distToEdge{4} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup4filtered));
    distToEdge{5} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup5filtered));
    distToEdge{6} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup6filtered));
    distToEdge{7} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup7filtered));
    distToEdge{8} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup8filtered));
    distToEdge{9} =arrayfun(@(x) mean(x.distToEdge),tracksNA(idGroup9filtered));
    [lengthLongest]=max(cellfun(@(x) length(x),distToEdge));

    matrixDistToEdge = NaN(lengthLongest,9);
    for ii=1:9
        matrixDistToEdge(1:length(distToEdge{ii}),ii) = distToEdge{ii};
    end
    boxWidth=0.5;
    figure
    boxplot(matrixDistToEdge,'orientation','vertical','whisker',0.5,'notch','on',...
        'labels',{'g1','g2','g3','g4','g5','g6','g7','g8','g9'},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
%     boxplot(matrixDistToEdge,'orientation','vertical','whisker',0.5,'notch','on',...
%         'labels',{['G1' '(N=' num2str(length(distToEdge{1})) ')'],['G2 (N=' num2str(length(distToEdge{2})) ')'],...
%         ['G3 (N=' num2str(length(distToEdge{3})) ')'],['G4 (N=' num2str(length(distToEdge{4})) ')'],...
%         ['G5 (N=' num2str(length(distToEdge{5})) ')'],['G6 (N=' num2str(length(distToEdge{6})) ')'],...
%         ['G7 (N=' num2str(length(distToEdge{7})) ')'],['G8 (N=' num2str(length(distToEdge{8})) ')'],...
%         ['G9 (N=' num2str(length(distToEdge{9})) ')']},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
%     boxplot(matrixDistToEdge,'orientation','vertical','whisker',0.5,'notch','on',...
%         'labels',{sprintf('G%d\n(N=%d)',1, length(distToEdge{1})),sprintf('G%d\n(N=%d)',2, length(distToEdge{2})),...
%         sprintf('G%d\n(N=%d)',3, length(distToEdge{3})),sprintf('G%d\n(N=%d)',4, length(distToEdge{4})),...
%         sprintf('G%d\n(N=%d)',5, length(distToEdge{5})),sprintf('G%d\n(N=%d)',6, length(distToEdge{6})),...
%         sprintf('G%d\n(N=%d)',7, length(distToEdge{7})),sprintf('G%d\n(N=%d)',8, length(distToEdge{8})),...
%         sprintf('G%d\n(N=%d)',9, length(distToEdge{9}))},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
    set(findobj(gca,'LineStyle','--'),'LineStyle','-')
    set(findobj(gca,'tag','Median'),'LineWidth',2)
    ylim([-2 50])
    title('distToEdge')
    save([pathForColocalization filesep 'data' filesep 'distToEdge.mat'],'distToEdge','matrixDistToEdge','-v7.3')
%     export_fig([pathForColocalization filesep 'eps' filesep 'distToEdgeForAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    print('-depsc', '-loose', [pathForColocalization filesep 'eps' filesep 'distToEdgeForAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    
    %% Look at feature difference per each group - earlyAmpSlope
    earlyAmpSlope{1} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup1filtered));
    earlyAmpSlope{2} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup2filtered));
    earlyAmpSlope{3} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup3filtered));
    earlyAmpSlope{4} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup4filtered));
    earlyAmpSlope{5} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup5filtered));
    earlyAmpSlope{6} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup6filtered));
    earlyAmpSlope{7} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup7filtered));
    earlyAmpSlope{8} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup8filtered));
    earlyAmpSlope{9} =arrayfun(@(x) (x.earlyAmpSlope),tracksNA(idGroup9filtered));
    [lengthLongestSlope]=max(cellfun(@(x) length(x),earlyAmpSlope));

    matrixEarlyAmpSlope = NaN(lengthLongestSlope,9);
    for ii=1:9
        matrixEarlyAmpSlope(1:length(earlyAmpSlope{ii}),ii) = earlyAmpSlope{ii};
    end
    boxWidth=0.5;
    figure
    boxplot(matrixEarlyAmpSlope,'orientation','vertical','whisker',0.5,'notch','on',...
        'labels',{'g1','g2','g3','g4','g5','g6','g7','g8','g9'},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
    set(findobj(gca,'LineStyle','--'),'LineStyle','-')
    set(findobj(gca,'tag','Median'),'LineWidth',2)
%     ylim([-2 50])
    title('earlyAmpSlope')
    save([pathForColocalization filesep 'data' filesep 'earlyAmpSlope.mat'],'earlyAmpSlope','matrixEarlyAmpSlope','-v7.3')
%     export_fig([pathForColocalization filesep 'eps' filesep 'earlyAmpSlopeAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    print('-depsc', '-loose', [pathForColocalization filesep 'eps' filesep 'earlyAmpSlopeAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    %% Look at feature difference per each group - advanceDist
    advanceDist{1} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup1filtered));
    advanceDist{2} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup2filtered));
    advanceDist{3} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup3filtered));
    advanceDist{4} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup4filtered));
    advanceDist{5} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup5filtered));
    advanceDist{6} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup6filtered));
    advanceDist{7} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup7filtered));
    advanceDist{8} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup8filtered));
    advanceDist{9} =arrayfun(@(x) mean(x.advanceDist),tracksNA(idGroup9filtered));
    [lengthLongestAdvDist]=max(cellfun(@(x) length(x),advanceDist));

    matrixAdvanceDist = NaN(lengthLongestAdvDist,9);
    for ii=1:9
        matrixAdvanceDist(1:length(advanceDist{ii}),ii) = advanceDist{ii};
    end
    boxWidth=0.5;
    figure
    boxplot(matrixAdvanceDist,'orientation','vertical','whisker',0.5,'notch','on',...
        'labels',{'g1','g2','g3','g4','g5','g6','g7','g8','g9'},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
    set(findobj(gca,'LineStyle','--'),'LineStyle','-')
    set(findobj(gca,'tag','Median'),'LineWidth',2)
%     ylim([-2 50])
    title('advanceDist')
    save([pathForColocalization filesep 'data' filesep 'advanceDist.mat'],'advanceDist','matrixAdvanceDist','-v7.3')
%     export_fig([pathForColocalization filesep 'eps' filesep 'advanceDistAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    print('-depsc', '-loose', [pathForColocalization filesep 'eps' filesep 'advanceDistAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    %% Look at feature difference per each group - ampTotal
    ampTotal{1} =arrayfun(@(x) mean(x.ampTotal),tracksNA(idGroup1filtered));
    ampTotal{2} =arrayfun(@(x) mean(x.ampTotal),tracksNA(idGroup2filtered));
    ampTotal{3} =arrayfun(@(x) mean(x.ampTotal),tracksNA(idGroup3filtered));
    ampTotal{4} =arrayfun(@(x) mean(x.ampTotal),tracksNA(idGroup4filtered));
    ampTotal{5} =arrayfun(@(x) mean(x.ampTotal),tracksNA(idGroup5filtered));
    ampTotal{6} =arrayfun(@(x) mean(x.ampTotal),tracksNA(idGroup6filtered));
    ampTotal{7} =arrayfun(@(x) mean(x.ampTotal),tracksNA(idGroup7filtered));
    ampTotal{8} =arrayfun(@(x) mean(x.ampTotal),tracksNA(idGroup8filtered));
    ampTotal{9} =arrayfun(@(x) mean(x.ampTotal),tracksNA(idGroup9filtered));
    [lengthLongestAmpTotal]=max(cellfun(@(x) length(x),ampTotal));

    matrixAmpTotal = NaN(lengthLongestAmpTotal,9);
    for ii=1:9
        matrixAmpTotal(1:length(ampTotal{ii}),ii) = ampTotal{ii};
    end
    boxWidth=0.5;
    figure
    boxplot(matrixAmpTotal,'orientation','vertical','whisker',0.5,'notch','off',...
        'labels',{'g1','g2','g3','g4','g5','g6','g7','g8','g9'},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
    set(findobj(gca,'LineStyle','--'),'LineStyle','-')
    set(findobj(gca,'tag','Median'),'LineWidth',2)
%     ylim([-2 50])
    title('ampTotal')
    save([pathForColocalization filesep 'data' filesep 'ampTotal.mat'],'ampTotal','matrixAmpTotal','-v7.3')
%     export_fig([pathForColocalization filesep 'eps' filesep 'ampTotalAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    print('-depsc', '-loose', [pathForColocalization filesep 'eps' filesep 'ampTotalAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    %% Look at feature difference per each group - starting ampTotal
    startingAmpTotal{1} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup1filtered));
    startingAmpTotal{2} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup2filtered));
    startingAmpTotal{3} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup3filtered));
    startingAmpTotal{4} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup4filtered));
    startingAmpTotal{5} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup5filtered));
    startingAmpTotal{6} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup6filtered));
    startingAmpTotal{7} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup7filtered));
    startingAmpTotal{8} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup8filtered));
    startingAmpTotal{9} =arrayfun(@(x) (x.ampTotal(x.startingFrameExtra)),tracksNA(idGroup9filtered));
    [lengthLongestStartingAmpTotal]=max(cellfun(@(x) length(x),startingAmpTotal));

    matrixStartingAmpTotal = NaN(lengthLongestStartingAmpTotal,9);
    for ii=1:9
        matrixStartingAmpTotal(1:length(startingAmpTotal{ii}),ii) = startingAmpTotal{ii};
    end
    boxWidth=0.5;
    figure
    boxplot(matrixStartingAmpTotal,'orientation','vertical','whisker',0.5,'notch','off',...
        'labels',{'g1','g2','g3','g4','g5','g6','g7','g8','g9'},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
    set(findobj(gca,'LineStyle','--'),'LineStyle','-')
    set(findobj(gca,'tag','Median'),'LineWidth',2)
%     ylim([-2 50])
    title('startingAmpTotal')
    save([pathForColocalization filesep 'data' filesep 'startingAmpTotal.mat'],'startingAmpTotal','matrixStartingAmpTotal','-v7.3')
%     export_fig([pathForColocalization filesep 'eps' filesep 'startingAmpTotalAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    print('-depsc', '-loose', [pathForColocalization filesep 'eps' filesep 'startingAmpTotalAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    %% Look at feature difference per each group - edgeAdvanceDistChange
    edgeAdvanceDistChange{1} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup1filtered));
    edgeAdvanceDistChange{2} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup2filtered));
    edgeAdvanceDistChange{3} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup3filtered));
    edgeAdvanceDistChange{4} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup4filtered));
    edgeAdvanceDistChange{5} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup5filtered));
    edgeAdvanceDistChange{6} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup6filtered));
    edgeAdvanceDistChange{7} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup7filtered));
    edgeAdvanceDistChange{8} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup8filtered));
    edgeAdvanceDistChange{9} =arrayfun(@(x) (x.edgeAdvanceDistChange2min(x.endingFrameExtra)),tracksNA(idGroup9filtered));
    [lengthLongestStartingAmpTotal]=max(cellfun(@(x) length(x),edgeAdvanceDistChange));

    matrixEdgeAdvanceDistChange = NaN(lengthLongestStartingAmpTotal,9);
    for ii=1:9
        matrixEdgeAdvanceDistChange(1:length(edgeAdvanceDistChange{ii}),ii) = edgeAdvanceDistChange{ii};
    end
    boxWidth=0.5;
    figure
    boxplot(matrixEdgeAdvanceDistChange,'orientation','vertical','whisker',0.5,'notch','off',...
        'labels',{'g1','g2','g3','g4','g5','g6','g7','g8','g9'},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
    set(findobj(gca,'LineStyle','--'),'LineStyle','-')
    set(findobj(gca,'tag','Median'),'LineWidth',2)
%     ylim([-2 50])
    title('edgeAdvanceDistChange')
    save([pathForColocalization filesep 'data' filesep 'edgeAdvanceDistChange.mat'],'edgeAdvanceDistChange','matrixEdgeAdvanceDistChange','-v7.3')
%     export_fig([pathForColocalization filesep 'eps' filesep 'edgeAdvanceDistChangeAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    print('-depsc', '-loose', [pathForColocalization filesep 'eps' filesep 'edgeAdvanceDistChangeAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
      %% Look at feature difference per each group - starting forceMag
    startingForceMag{1} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup1filtered));
    startingForceMag{2} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup2filtered));
    startingForceMag{3} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup3filtered));
    startingForceMag{4} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup4filtered));
    startingForceMag{5} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup5filtered));
    startingForceMag{6} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup6filtered));
    startingForceMag{7} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup7filtered));
    startingForceMag{8} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup8filtered));
    startingForceMag{9} =arrayfun(@(x) (x.forceMag(x.startingFrameExtra)),tracksNA(idGroup9filtered));
    [lengthLongestStartingForceMag]=max(cellfun(@(x) length(x),startingForceMag));

    matrixStartingForceMag = NaN(lengthLongestStartingForceMag,9);
    for ii=1:9
        matrixStartingForceMag(1:length(startingForceMag{ii}),ii) = startingForceMag{ii};
    end
    boxWidth=0.5;
    figure
    boxplot(matrixStartingForceMag,'orientation','vertical','whisker',0.5,'notch','off',...
        'labels',{'g1','g2','g3','g4','g5','g6','g7','g8','g9'},'symbol','','widths',boxWidth,'jitter',1,'colors','k')
    set(findobj(gca,'LineStyle','--'),'LineStyle','-')
    set(findobj(gca,'tag','Median'),'LineWidth',2)
%     ylim([-2 50])
    title('startingForceMag')
    save([pathForColocalization filesep 'data' filesep 'startingForceMag.mat'],'startingForceMag','matrixStartingForceMag','-v7.3')
%     export_fig([pathForColocalization filesep 'eps' filesep 'startingForceMagAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    print('-depsc', '-loose', [pathForColocalization filesep 'eps' filesep 'startingForceMagAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    %% Look at feature difference per each group - force slope
    forceSlope{1} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup1filtered));
    forceSlope{2} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup2filtered));
    forceSlope{3} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup3filtered));
    forceSlope{4} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup4filtered));
    forceSlope{5} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup5filtered));
    forceSlope{6} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup6filtered));
    forceSlope{7} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup7filtered));
    forceSlope{8} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup8filtered));
    forceSlope{9} =arrayfun(@(x) (x.forceSlope),tracksNA(idGroup9filtered));
    [lengthForceSlope]=max(cellfun(@(x) length(x),forceSlope));

    matrixForceSlope = NaN(lengthForceSlope,9);
    for ii=1:9
        matrixForceSlope(1:length(forceSlope{ii}),ii) = forceSlope{ii};
    end
    boxWidth=0.5;
    figure
    boxplot(matrixForceSlope,'orientation','vertical','whisker',1,'notch','off',...
        'labels',{'g1','g2','g3','g4','g5','g6','g7','g8','g9'},'symbol','','widths',boxWidth,'colors','k')
    set(findobj(gca,'LineStyle','--'),'LineStyle','-')
    set(findobj(gca,'tag','Median'),'LineWidth',2)
%     ylim([-2 50])
    title('forceSlope')
    save([pathForColocalization filesep 'data' filesep 'forceSlope.mat'],'forceSlope','matrixForceSlope','-v7.3')
%     export_fig([pathForColocalization filesep 'eps' filesep 'forceSlopeAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    print('-depsc', '-loose', [pathForColocalization filesep 'eps' filesep 'forceSlopeAllGroups.eps']);% histogramPeakLagVinVsTal -transparent
    %% Save all workspace
    save([pathForColocalization filesep 'data' filesep 'allDataAfterClassification.mat'],'-v7.3')