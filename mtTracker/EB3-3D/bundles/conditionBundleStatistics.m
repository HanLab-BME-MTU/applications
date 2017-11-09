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
	parfor pIdx=1:length(scoringProcessCell{cIdx})
		tmp=load(scoringProcessCell{cIdx}(pIdx).outFilePaths_{1});
		condScoreCell{pIdx}=tmp.zScores(:);
	end
	scoreCell{cIdx}=condScoreCell;
end

%% Maria's format
allCell=[scoreCell{:}];
dataMatPadded = reformatDataCell(allCell);
condGroups=arrayfun(@(cIdx) cIdx*ones(1,length(scoreCell{cIdx})),1:length(scoreCell),'unif',0);
condGroups = horzcat(condGroups{:}); 
[FCondMean,FCondPooled,FPerCell]=plotDifferentPerspectives(dataMatPadded,condGroups);
printPNGEPSFIG(FCondMean,outputDirPlot,'ZScoresCondMean');
printPNGEPSFIG(FCondPooled,outputDirPlot,'ZScoresCondPooled');
printPNGEPSFIG(FPerCell,outputDirPlot,'ZScoresPerCell');

%%
c={'r','b','g','y','k'};
% Cellular average zScores distribution
avgScore=[];
avgScoreGroup=[];
for cIdx=1:length(scoringProcessCell)
	for pIdx=1:length(scoringProcessCell{cIdx})
		avgScore=[avgScore nanmean(scoreCell{cIdx}{pIdx}(:))];		
        avgScoreGroup=[avgScoreGroup cIdx];      
	end
end
F=figure;
legend(names);
notBoxPlot(avgScore,avgScoreGroup);
ylabel('Average Bundle ZScore');
printPNGEPSFIG(F,outputDirPlot,'avgBundleScorePerCell');

% Pooled zScores distribution
scoresBin=-10:0.2:10;
zScore=[];
zScoreGroup=[];
histCell=cell(1,length(scoringProcessCell));
[Handles,~,F]=setupFigure(1,2,2);
hold on;
for cIdx=1:length(scoringProcessCell)
	for pIdx=1:length(scoringProcessCell{cIdx})
		counts=histcounts(scoreCell{cIdx}{pIdx},scoresBin);
		histCell{cIdx}=[histCell{cIdx}; counts];
		zScore=[zScore scoreCell{cIdx}{pIdx}(:)'];
		zScoreGroup=[zScoreGroup cIdx*ones(size(scoreCell{cIdx}{pIdx}(:)'))];
	end
end
legend(names)
hold off;
axes(Handles(1))
boxplot(zScore,zScoreGroup);
[meanZScore,~,stdZscore]=statPerIndx(zScore,zScoreGroup);
axes(Handles(2))
H=shadedErrorBar(1:length(meanZScore),meanZScore,stdZscore);
ylabel('Pooled Bundle ZScore');
printPNGEPSFIG(F,outputDirPlot,'pooledBundleScorePerCell');


% Plot fibered KT count (zscore>2)
thresholds=[0.1 0.25 0.5 1 2];
for thresh=thresholds
fiberCount=[];
fiberGroup=[];
F=figure
hold on;
for cIdx=1:length(scoringProcessCell)
	for pIdx=1:length(scoringProcessCell{cIdx})
        scores=scoreCell{cIdx}{pIdx};
		counts=sum(scores(~isnan(scores))>thresh);
		fiberCount=[fiberCount counts];
        fiberGroup=[fiberGroup cIdx];
    end	
end
% boxplot(fiberCount,fiberGroup,'labels',names)
notBoxPlot(fiberCount,fiberGroup);
ylabel(['Est. Bundled Kinetochore count (Z-score > ' num2str(thresh)]);
printPNGEPSFIG(F,outputDirPlot,['fiberCount-' num2str(thresh)]);
end

% Plot fibered KT percentage (zscore>2)
percentagePosScore=zeros(1,length(scoringProcessCell));
percentageNegScore=zeros(1,length(scoringProcessCell));
F=figure;
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
minutes=[1 4 8 16];
plot(minutes(1:length(scoringProcessCell)),percentagePosScore)
plot(minutes(1:length(scoringProcessCell)),percentageNegScore)
legend({'Pos. Score','Neg. Score'})
hold off;
printPNGEPSFIG(F,outputDirPlot,'perc');


