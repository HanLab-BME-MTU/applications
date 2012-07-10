% 
% % criteria.lifeTime.min = 5;
% % indx5 = chooseTracks(tracksFinal,criteria);
% % tracksFinal = tracksFinal(indx5);
% % diffAnalysisRes = diffAnalysisRes(indx5);
% 
% [sptPropInWindow,windowTrackAssign,trackWindowAssign,windowSize] = ...
%     particleBehaviorRelToActivityOnset(tracksFinal,windowsAll,1:400:7201,...
%     protSamples,diffAnalysisRes,5,[1 1; 2 2; 3 3],[],[],[],[],diffModeAnalysisRes);
% 
% [sptPropInWindow,windowTrackAssign,trackWindowAssign,windowSize] = ...
%     particleBehaviorRelToActivityOnset(tracksFinal,windowsAll,1:400:7201,...
%     protSamples,diffAnalysisRes,5,[1 1; 2 2; 3 3],[],[],windowTrackAssign,trackWindowAssign);
% 
% save('particleBehavior','sptPropInWindow');
% 
% % save('tracksDiffLength5','tracksFinal','diffAnalysisRes','indx5');
% % save('windowsActivityTracks','protSamples','trackWindowAssign','windowSize','windowTrackAssign','windowsAll');

meanTmp = sptPropInWindow(1).modeAnalysis.fracMode3.mean;
stdTmp = sptPropInWindow(1).modeAnalysis.fracMode3.std;
numTmp = sptPropInWindow(1).modeAnalysis.fracMode3.numPoints;
meanTmp(numTmp<25) = NaN;
stdTmp(numTmp<25) = NaN;

figHandle = figure; hold on
plot(-3:6,meanTmp(1,1:10),'b','Marker','.')
myErrorbar(-3:6,meanTmp(1,1:10),stdTmp(1,1:10)./sqrt(numTmp(1,1:10)))
plot(-3:6,meanTmp(2,1:10),'c','Marker','.')
myErrorbar(-3:6,meanTmp(2,1:10),stdTmp(2,1:10)./sqrt(numTmp(2,1:10)))
plot(-3:6,meanTmp(3,1:10),'r','Marker','.')
myErrorbar(-3:6,meanTmp(3,1:10),stdTmp(3,1:10)./sqrt(numTmp(3,1:10)))
plot(-3:6,meanTmp(4,1:10),'r--','Marker','.')
myErrorbar(-3:6,meanTmp(4,1:10),stdTmp(4,1:10)./sqrt(numTmp(4,1:10)))
plot(-3:6,meanTmp(5,1:10),'r:','Marker','.')
myErrorbar(-3:6,meanTmp(5,1:10),stdTmp(5,1:10)./sqrt(numTmp(5,1:10)))
plot(-3:6,meanTmp(6,1:10),'g','Marker','.')
myErrorbar(-3:6,meanTmp(6,1:10),stdTmp(6,1:10)./sqrt(numTmp(6,1:10)))
plot(-3:6,meanTmp(7,1:10),'g--','Marker','.')
myErrorbar(-3:6,meanTmp(7,1:10),stdTmp(7,1:10)./sqrt(numTmp(7,1:10)))
plot(-3:6,meanTmp(8,1:10),'g:','Marker','.')
myErrorbar(-3:6,meanTmp(8,1:10),stdTmp(8,1:10)./sqrt(numTmp(8,1:10)))
plot(-3:6,meanTmp(9,1:10),'k','Marker','.','LineWidth',2,'MarkerSize',10)
myErrorbar(-3:6,meanTmp(9,1:10),stdTmp(9,1:10)./sqrt(numTmp(9,1:10)))
legend({'pr>re','pr>un','pr>pa(s)>pr','pr>pa(s)>re','pr>pa(s)>un',...
    'pr>pa(l)>pr','pr>pa(l)>re','pr>pa(l)>un','all'})
% saveas(figHandle,'protrusion_f2fDisp_band1','fig');

figHandle = figure; hold on
plot(-3:6,meanTmp(10,1:10),'b','Marker','.')
myErrorbar(-3:6,meanTmp(10,1:10),stdTmp(10,1:10)./sqrt(numTmp(10,1:10)))
plot(-3:6,meanTmp(11,1:10),'c','Marker','.')
myErrorbar(-3:6,meanTmp(11,1:10),stdTmp(11,1:10)./sqrt(numTmp(11,1:10)))
plot(-3:6,meanTmp(12,1:10),'r','Marker','.')
myErrorbar(-3:6,meanTmp(12,1:10),stdTmp(12,1:10)./sqrt(numTmp(12,1:10)))
plot(-3:6,meanTmp(13,1:10),'r--','Marker','.')
myErrorbar(-3:6,meanTmp(13,1:10),stdTmp(13,1:10)./sqrt(numTmp(13,1:10)))
plot(-3:6,meanTmp(14,1:10),'r:','Marker','.')
myErrorbar(-3:6,meanTmp(14,1:10),stdTmp(14,1:10)./sqrt(numTmp(14,1:10)))
plot(-3:6,meanTmp(15,1:10),'g','Marker','.')
myErrorbar(-3:6,meanTmp(15,1:10),stdTmp(15,1:10)./sqrt(numTmp(15,1:10)))
plot(-3:6,meanTmp(16,1:10),'g--','Marker','.')
myErrorbar(-3:6,meanTmp(16,1:10),stdTmp(16,1:10)./sqrt(numTmp(16,1:10)))
plot(-3:6,meanTmp(17,1:10),'g:','Marker','.')
myErrorbar(-3:6,meanTmp(17,1:10),stdTmp(17,1:10)./sqrt(numTmp(17,1:10)))
plot(-3:6,meanTmp(18,1:10),'k','Marker','.','LineWidth',2,'MarkerSize',10)
myErrorbar(-3:6,meanTmp(18,1:10),stdTmp(18,1:10)./sqrt(numTmp(18,1:10)))
legend({'re>pr','re>un','re>pa(s)>pr','re>pa(s)>re','re>pa(s)>un',...
    're>pa(l)>pr','re>pa(l)>re','re>pa(l)>un','all'})
% saveas(figHandle,'retraction_f2fDisp_band1','fig');

figHandle = figure; hold on
plot(-3:6,meanTmp(19,1:10),'b','Marker','.')
myErrorbar(-3:6,meanTmp(19,1:10),stdTmp(19,1:10)./sqrt(numTmp(19,1:10)))
plot(-3:6,meanTmp(20,1:10),'b--','Marker','.')
myErrorbar(-3:6,meanTmp(20,1:10),stdTmp(20,1:10)./sqrt(numTmp(20,1:10)))
plot(-3:6,meanTmp(21,1:10),'b:','Marker','.')
myErrorbar(-3:6,meanTmp(21,1:10),stdTmp(21,1:10)./sqrt(numTmp(21,1:10)))
plot(-3:6,meanTmp(22,1:10),'r','Marker','.')
myErrorbar(-3:6,meanTmp(22,1:10),stdTmp(22,1:10)./sqrt(numTmp(22,1:10)))
plot(-3:6,meanTmp(23,1:10),'r--','Marker','.')
myErrorbar(-3:6,meanTmp(23,1:10),stdTmp(23,1:10)./sqrt(numTmp(23,1:10)))
plot(-3:6,meanTmp(24,1:10),'r:','Marker','.')
myErrorbar(-3:6,meanTmp(24,1:10),stdTmp(24,1:10)./sqrt(numTmp(24,1:10)))
plot(-3:6,meanTmp(25,1:10),'g','Marker','.')
myErrorbar(-3:6,meanTmp(25,1:10),stdTmp(25,1:10)./sqrt(numTmp(25,1:10)))
plot(-3:6,meanTmp(26,1:10),'g--','Marker','.')
myErrorbar(-3:6,meanTmp(26,1:10),stdTmp(26,1:10)./sqrt(numTmp(26,1:10)))
plot(-3:6,meanTmp(27,1:10),'g:','Marker','.')
myErrorbar(-3:6,meanTmp(27,1:10),stdTmp(27,1:10)./sqrt(numTmp(27,1:10)))
legend({'pa>pr&<pr','pa>pr&<re','pa>pr&<un','pa>re&<pr','pa>re&<re',...
    'pa>re&<un','pa>un&<pr','pa>un&<re','pa>un&<un'})
% saveas(figHandle,'pause_f2fDisp_band1','fig');

meanTmp = sptPropInWindow(2).modeAnalysis.fracMode3.mean;
stdTmp = sptPropInWindow(2).modeAnalysis.fracMode3.std;
numTmp = sptPropInWindow(2).modeAnalysis.fracMode3.numPoints;
meanTmp(numTmp<25) = NaN;
stdTmp(numTmp<25) = NaN;

figHandle = figure; hold on
plot(-3:6,meanTmp(1,1:10),'b','Marker','.')
myErrorbar(-3:6,meanTmp(1,1:10),stdTmp(1,1:10)./sqrt(numTmp(1,1:10)))
plot(-3:6,meanTmp(2,1:10),'c','Marker','.')
myErrorbar(-3:6,meanTmp(2,1:10),stdTmp(2,1:10)./sqrt(numTmp(2,1:10)))
plot(-3:6,meanTmp(3,1:10),'r','Marker','.')
myErrorbar(-3:6,meanTmp(3,1:10),stdTmp(3,1:10)./sqrt(numTmp(3,1:10)))
plot(-3:6,meanTmp(4,1:10),'r--','Marker','.')
myErrorbar(-3:6,meanTmp(4,1:10),stdTmp(4,1:10)./sqrt(numTmp(4,1:10)))
plot(-3:6,meanTmp(5,1:10),'r:','Marker','.')
myErrorbar(-3:6,meanTmp(5,1:10),stdTmp(5,1:10)./sqrt(numTmp(5,1:10)))
plot(-3:6,meanTmp(6,1:10),'g','Marker','.')
myErrorbar(-3:6,meanTmp(6,1:10),stdTmp(6,1:10)./sqrt(numTmp(6,1:10)))
plot(-3:6,meanTmp(7,1:10),'g--','Marker','.')
myErrorbar(-3:6,meanTmp(7,1:10),stdTmp(7,1:10)./sqrt(numTmp(7,1:10)))
plot(-3:6,meanTmp(8,1:10),'g:','Marker','.')
myErrorbar(-3:6,meanTmp(8,1:10),stdTmp(8,1:10)./sqrt(numTmp(8,1:10)))
plot(-3:6,meanTmp(9,1:10),'k','Marker','.','LineWidth',2,'MarkerSize',10)
myErrorbar(-3:6,meanTmp(9,1:10),stdTmp(9,1:10)./sqrt(numTmp(9,1:10)))
legend({'pr>re','pr>un','pr>pa(s)>pr','pr>pa(s)>re','pr>pa(s)>un',...
    'pr>pa(l)>pr','pr>pa(l)>re','pr>pa(l)>un','all'})
% saveas(figHandle,'protrusion_f2fDisp_band2','fig');

figHandle = figure; hold on
plot(-3:6,meanTmp(10,1:10),'b','Marker','.')
myErrorbar(-3:6,meanTmp(10,1:10),stdTmp(10,1:10)./sqrt(numTmp(10,1:10)))
plot(-3:6,meanTmp(11,1:10),'c','Marker','.')
myErrorbar(-3:6,meanTmp(11,1:10),stdTmp(11,1:10)./sqrt(numTmp(11,1:10)))
plot(-3:6,meanTmp(12,1:10),'r','Marker','.')
myErrorbar(-3:6,meanTmp(12,1:10),stdTmp(12,1:10)./sqrt(numTmp(12,1:10)))
plot(-3:6,meanTmp(13,1:10),'r--','Marker','.')
myErrorbar(-3:6,meanTmp(13,1:10),stdTmp(13,1:10)./sqrt(numTmp(13,1:10)))
plot(-3:6,meanTmp(14,1:10),'r:','Marker','.')
myErrorbar(-3:6,meanTmp(14,1:10),stdTmp(14,1:10)./sqrt(numTmp(14,1:10)))
plot(-3:6,meanTmp(15,1:10),'g','Marker','.')
myErrorbar(-3:6,meanTmp(15,1:10),stdTmp(15,1:10)./sqrt(numTmp(15,1:10)))
plot(-3:6,meanTmp(16,1:10),'g--','Marker','.')
myErrorbar(-3:6,meanTmp(16,1:10),stdTmp(16,1:10)./sqrt(numTmp(16,1:10)))
plot(-3:6,meanTmp(17,1:10),'g:','Marker','.')
myErrorbar(-3:6,meanTmp(17,1:10),stdTmp(17,1:10)./sqrt(numTmp(17,1:10)))
plot(-3:6,meanTmp(18,1:10),'k','Marker','.','LineWidth',2,'MarkerSize',10)
myErrorbar(-3:6,meanTmp(18,1:10),stdTmp(18,1:10)./sqrt(numTmp(18,1:10)))
legend({'re>pr','re>un','re>pa(s)>pr','re>pa(s)>re','re>pa(s)>un',...
    're>pa(l)>pr','re>pa(l)>re','re>pa(l)>un','all'})
% saveas(figHandle,'retraction_f2fDisp_band2','fig');

figHandle = figure; hold on
plot(-3:6,meanTmp(19,1:10),'b','Marker','.')
myErrorbar(-3:6,meanTmp(19,1:10),stdTmp(19,1:10)./sqrt(numTmp(19,1:10)))
plot(-3:6,meanTmp(20,1:10),'b--','Marker','.')
myErrorbar(-3:6,meanTmp(20,1:10),stdTmp(20,1:10)./sqrt(numTmp(20,1:10)))
plot(-3:6,meanTmp(21,1:10),'b:','Marker','.')
myErrorbar(-3:6,meanTmp(21,1:10),stdTmp(21,1:10)./sqrt(numTmp(21,1:10)))
plot(-3:6,meanTmp(22,1:10),'r','Marker','.')
myErrorbar(-3:6,meanTmp(22,1:10),stdTmp(22,1:10)./sqrt(numTmp(22,1:10)))
plot(-3:6,meanTmp(23,1:10),'r--','Marker','.')
myErrorbar(-3:6,meanTmp(23,1:10),stdTmp(23,1:10)./sqrt(numTmp(23,1:10)))
plot(-3:6,meanTmp(24,1:10),'r:','Marker','.')
myErrorbar(-3:6,meanTmp(24,1:10),stdTmp(24,1:10)./sqrt(numTmp(24,1:10)))
plot(-3:6,meanTmp(25,1:10),'g','Marker','.')
myErrorbar(-3:6,meanTmp(25,1:10),stdTmp(25,1:10)./sqrt(numTmp(25,1:10)))
plot(-3:6,meanTmp(26,1:10),'g--','Marker','.')
myErrorbar(-3:6,meanTmp(26,1:10),stdTmp(26,1:10)./sqrt(numTmp(26,1:10)))
plot(-3:6,meanTmp(27,1:10),'g:','Marker','.')
myErrorbar(-3:6,meanTmp(27,1:10),stdTmp(27,1:10)./sqrt(numTmp(27,1:10)))
legend({'pa>pr&<pr','pa>pr&<re','pa>pr&<un','pa>re&<pr','pa>re&<re',...
    'pa>re&<un','pa>un&<pr','pa>un&<re','pa>un&<un'})
% saveas(figHandle,'pause_f2fDisp_band2','fig');

meanTmp = sptPropInWindow(3).modeAnalysis.fracMode3.mean;
stdTmp = sptPropInWindow(3).modeAnalysis.fracMode3.std;
numTmp = sptPropInWindow(3).modeAnalysis.fracMode3.numPoints;
meanTmp(numTmp<25) = NaN;
stdTmp(numTmp<25) = NaN;

figHandle = figure; hold on
plot(-3:6,meanTmp(1,1:10),'b','Marker','.')
myErrorbar(-3:6,meanTmp(1,1:10),stdTmp(1,1:10)./sqrt(numTmp(1,1:10)))
plot(-3:6,meanTmp(2,1:10),'c','Marker','.')
myErrorbar(-3:6,meanTmp(2,1:10),stdTmp(2,1:10)./sqrt(numTmp(2,1:10)))
plot(-3:6,meanTmp(3,1:10),'r','Marker','.')
myErrorbar(-3:6,meanTmp(3,1:10),stdTmp(3,1:10)./sqrt(numTmp(3,1:10)))
plot(-3:6,meanTmp(4,1:10),'r--','Marker','.')
myErrorbar(-3:6,meanTmp(4,1:10),stdTmp(4,1:10)./sqrt(numTmp(4,1:10)))
plot(-3:6,meanTmp(5,1:10),'r:','Marker','.')
myErrorbar(-3:6,meanTmp(5,1:10),stdTmp(5,1:10)./sqrt(numTmp(5,1:10)))
plot(-3:6,meanTmp(6,1:10),'g','Marker','.')
myErrorbar(-3:6,meanTmp(6,1:10),stdTmp(6,1:10)./sqrt(numTmp(6,1:10)))
plot(-3:6,meanTmp(7,1:10),'g--','Marker','.')
myErrorbar(-3:6,meanTmp(7,1:10),stdTmp(7,1:10)./sqrt(numTmp(7,1:10)))
plot(-3:6,meanTmp(8,1:10),'g:','Marker','.')
myErrorbar(-3:6,meanTmp(8,1:10),stdTmp(8,1:10)./sqrt(numTmp(8,1:10)))
plot(-3:6,meanTmp(9,1:10),'k','Marker','.','LineWidth',2,'MarkerSize',10)
myErrorbar(-3:6,meanTmp(9,1:10),stdTmp(9,1:10)./sqrt(numTmp(9,1:10)))
legend({'pr>re','pr>un','pr>pa(s)>pr','pr>pa(s)>re','pr>pa(s)>un',...
    'pr>pa(l)>pr','pr>pa(l)>re','pr>pa(l)>un','all'})
% saveas(figHandle,'protrusion_f2fDisp_band3','fig');

figHandle = figure; hold on
plot(-3:6,meanTmp(10,1:10),'b','Marker','.')
myErrorbar(-3:6,meanTmp(10,1:10),stdTmp(10,1:10)./sqrt(numTmp(10,1:10)))
plot(-3:6,meanTmp(11,1:10),'c','Marker','.')
myErrorbar(-3:6,meanTmp(11,1:10),stdTmp(11,1:10)./sqrt(numTmp(11,1:10)))
plot(-3:6,meanTmp(12,1:10),'r','Marker','.')
myErrorbar(-3:6,meanTmp(12,1:10),stdTmp(12,1:10)./sqrt(numTmp(12,1:10)))
plot(-3:6,meanTmp(13,1:10),'r--','Marker','.')
myErrorbar(-3:6,meanTmp(13,1:10),stdTmp(13,1:10)./sqrt(numTmp(13,1:10)))
plot(-3:6,meanTmp(14,1:10),'r:','Marker','.')
myErrorbar(-3:6,meanTmp(14,1:10),stdTmp(14,1:10)./sqrt(numTmp(14,1:10)))
plot(-3:6,meanTmp(15,1:10),'g','Marker','.')
myErrorbar(-3:6,meanTmp(15,1:10),stdTmp(15,1:10)./sqrt(numTmp(15,1:10)))
plot(-3:6,meanTmp(16,1:10),'g--','Marker','.')
myErrorbar(-3:6,meanTmp(16,1:10),stdTmp(16,1:10)./sqrt(numTmp(16,1:10)))
plot(-3:6,meanTmp(17,1:10),'g:','Marker','.')
myErrorbar(-3:6,meanTmp(17,1:10),stdTmp(17,1:10)./sqrt(numTmp(17,1:10)))
plot(-3:6,meanTmp(18,1:10),'k','Marker','.','LineWidth',2,'MarkerSize',10)
myErrorbar(-3:6,meanTmp(18,1:10),stdTmp(18,1:10)./sqrt(numTmp(18,1:10)))
legend({'re>pr','re>un','re>pa(s)>pr','re>pa(s)>re','re>pa(s)>un',...
    're>pa(l)>pr','re>pa(l)>re','re>pa(l)>un','all'})
% saveas(figHandle,'retraction_f2fDisp_band3','fig');

figHandle = figure; hold on
plot(-3:6,meanTmp(19,1:10),'b','Marker','.')
myErrorbar(-3:6,meanTmp(19,1:10),stdTmp(19,1:10)./sqrt(numTmp(19,1:10)))
plot(-3:6,meanTmp(20,1:10),'b--','Marker','.')
myErrorbar(-3:6,meanTmp(20,1:10),stdTmp(20,1:10)./sqrt(numTmp(20,1:10)))
plot(-3:6,meanTmp(21,1:10),'b:','Marker','.')
myErrorbar(-3:6,meanTmp(21,1:10),stdTmp(21,1:10)./sqrt(numTmp(21,1:10)))
plot(-3:6,meanTmp(22,1:10),'r','Marker','.')
myErrorbar(-3:6,meanTmp(22,1:10),stdTmp(22,1:10)./sqrt(numTmp(22,1:10)))
plot(-3:6,meanTmp(23,1:10),'r--','Marker','.')
myErrorbar(-3:6,meanTmp(23,1:10),stdTmp(23,1:10)./sqrt(numTmp(23,1:10)))
plot(-3:6,meanTmp(24,1:10),'r:','Marker','.')
myErrorbar(-3:6,meanTmp(24,1:10),stdTmp(24,1:10)./sqrt(numTmp(24,1:10)))
plot(-3:6,meanTmp(25,1:10),'g','Marker','.')
myErrorbar(-3:6,meanTmp(25,1:10),stdTmp(25,1:10)./sqrt(numTmp(25,1:10)))
plot(-3:6,meanTmp(26,1:10),'g--','Marker','.')
myErrorbar(-3:6,meanTmp(26,1:10),stdTmp(26,1:10)./sqrt(numTmp(26,1:10)))
plot(-3:6,meanTmp(27,1:10),'g:','Marker','.')
myErrorbar(-3:6,meanTmp(27,1:10),stdTmp(27,1:10)./sqrt(numTmp(27,1:10)))
legend({'pa>pr&<pr','pa>pr&<re','pa>pr&<un','pa>re&<pr','pa>re&<re',...
    'pa>re&<un','pa>un&<pr','pa>un&<re','pa>un&<un'})
% saveas(figHandle,'pause_f2fDisp_band3','fig');

