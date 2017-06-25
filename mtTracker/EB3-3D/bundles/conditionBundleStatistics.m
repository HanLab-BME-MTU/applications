function conditionBundleStatistics(scoringProcessCell)
% load scoring data associated to cells. A Cell of process list describes the conditions

if(~iscell(scoringProcessCell))
    scoringProcessCell={scoringProcessCell};
end

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
hold off;
