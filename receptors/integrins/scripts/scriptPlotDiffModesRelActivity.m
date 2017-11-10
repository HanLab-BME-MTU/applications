meanTmp = sptPropInWindow(1).modeAnalysis.fracMode1.mean;
stdTmp = sptPropInWindow(1).modeAnalysis.fracMode1.std;
numTmp = sptPropInWindow(1).modeAnalysis.fracMode1.numPoints;
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


meanTmp = sptPropInWindow(1).modeAnalysis.fracMode2.mean;
stdTmp = sptPropInWindow(1).modeAnalysis.fracMode2.std;
numTmp = sptPropInWindow(1).modeAnalysis.fracMode2.numPoints;
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


meanTmp = sptPropInWindow(1).modeAnalysis.fracMode4.mean;
stdTmp = sptPropInWindow(1).modeAnalysis.fracMode4.std;
numTmp = sptPropInWindow(1).modeAnalysis.fracMode4.numPoints;
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
