function conditionBundleStatistics(scoringProcessCell,names,outputDirPlot)
% load scoring data associated to cells. A Cell of process list describes the conditions

if(~iscell(scoringProcessCell))
    scoringProcessCell={scoringProcessCell};
end

mkdirRobust(outputDirPlot);

% load score set in cells of cells
scoreCell=cell(1,length(scoringProcessCell));
for cIdx=1:length(scoringProcessCell)
	condScoreCell=cell(1,length(scoringProcessCell{cIdx}));
	for pIdx=1:length(scoringProcessCell{cIdx})
		tmp=load(scoringProcessCell{cIdx}(pIdx).outFilePaths_{1});
		condScoreCell{pIdx}=tmp.zScores;
	end
	scoreCell{cIdx}=condScoreCell;
end

c={'r','b','g','y','k'};

% Plot zScores distribution
scoresBin=-10:0.2:10;
histCell=cell(1,length(scoringProcessCell));
hold on;
for cIdx=1:length(scoringProcessCell)
	for pIdx=1:length(scoringProcessCell{cIdx})
		counts=histcounts(scoreCell{cIdx}{pIdx},scoresBin);
		histCell{cIdx}=[histCell{cIdx}; counts];
	end
	
	if(size(histCell{cIdx},1)>1)
        y=histCell{cIdx};
        shadedErrorBar(scoresBin(1:end-1),mean(y),std(y),c{cIdx},1);       
    else
	    plot(scoresBin(1:end-1),histCell{cIdx},[c{cIdx} '-']);
    end
end
legend(names)
hold off;
print([outputDirPlot  'scoreDist.png'],'-dpng');
print([outputDirPlot  'scoreDist.eps'],'-depsc');


% Plot fibered KT count (zscore>2)
fiberCount=[];
fiberGroup=[];
thresh=2
figure
hold on;
for cIdx=1:length(scoringProcessCell)
	for pIdx=1:length(scoringProcessCell{cIdx})
        scores=scoreCell{cIdx}{pIdx};
		counts=sum(scores(~isnan(scores))>thresh);
		fiberCount=[fiberCount counts];
        fiberGroup=[fiberGroup cIdx];
    end	
end
boxplot(fiberCount,fiberGroup,'labels',names)
ylabel('Est. Bundled Kinetochore count (Z-score > 2)');
hold off;

print([outputDirPlot  'fiberCount.png'],'-dpng');
print([outputDirPlot  'fiberCount.eps'],'-depsc');

% Plot fibered KT count (zscore>2)
percentagePosScore=zeros(1,length(scoringProcessCell));
percentageNegScore=zeros(1,length(scoringProcessCell));
figure
hold on;
for cIdx=1:length(scoringProcessCell)
    countPos=zeros(1,length(scoringProcessCell{cIdx}));
    countNeg=zeros(1,length(scoringProcessCell{cIdx}));
	for pIdx=1:length(scoringProcessCell{cIdx})
        scores=scoreCell{cIdx}{pIdx};
		countPos(pIdx)=sum(scores(~isnan(scores))>0);
        countNeg(pIdx)=sum(scores(~isnan(scores))<0);
    end	
    percentagePosScore(cIdx)=mean(100*countPos./(countPos+countNeg));
    percentageNegScore(cIdx)=mean(100*countNeg./(countPos+countNeg));
    
end
plot([1 4 8 16],percentagePosScore)
plot([1 4 8 16],percentageNegScore)

legend({'Pos. Score','Neg. Score'})
hold off;

print([outputDirPlot  'perc.png'],'-dpng');
print([outputDirPlot  'perc.eps'],'-depsc');



