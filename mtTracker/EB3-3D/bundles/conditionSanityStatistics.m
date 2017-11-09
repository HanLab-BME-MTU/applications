function conditionSanityStatistics(processCell,names,outputDirPlot)
% load scoring data associated to cells. A Cell of process list describes the conditions

if(~iscell(processCell))
    processCell={processCell};
end

mkdirRobust(outputDirPlot);

% load score set in cells of cells
medianSpeedCell=cell(1,length(processCell));
for cIdx=1:length(processCell)
	condScoreCell=cell(1,length(processCell{cIdx}));
	for pIdx=1:length(processCell{cIdx})
		tmp=load(processCell{cIdx}(pIdx).outFilePaths_{1});
		condScoreCell{pIdx}=tmp.medspeed;
	end
	medianSpeedCell{cIdx}=condScoreCell;
end

c={'r','b','g','y','k'};


% Plot fibered KT count (zscore>2)
speed=[];
group=[];
figure
hold on;
for cIdx=1:length(processCell)
	for pIdx=1:length(processCell{cIdx})
        cellSpeed=medianSpeedCell{cIdx}{pIdx};
		speed=[speed median(cellSpeed)];
        group=[group cIdx];
    end	
end
boxplot(speed,group,'labels',names)
ylabel('Median growth rate (\mu/s)');
hold off;

print([outputDirPlot  'fiberCount.png'],'-dpng');
print([outputDirPlot  'fiberCount.eps'],'-depsc');



