function conditionLifetimeStatistics(processTrackCell,names,outputDirPlot)
% load scoring data associated to cells. A Cell of process list describes the conditions

if(~iscell(processTrackCell))
    processTrackCell={processTrackCell};
end

mkdirRobust(outputDirPlot);
maxLifetime=1;
maxIntensity=1;


% load lifetime data set in cells of cells
lifetimeCell=cell(1,length(processTrackCell));
uncutLifetimeCell=cell(1,length(processTrackCell));
maxIntensitiesCell=cell(1,length(processTrackCell));

tic;
for cIdx=1:length(processTrackCell)
	condLifetimeCell=cell(1,length(processTrackCell{cIdx}));
 	condUCLifetimeCell=cell(1,length(processTrackCell{cIdx}));
 	condMaxIntensitiesCell=cell(1,length(processTrackCell{cIdx}));

	parfor pIdx=1:length(processTrackCell{cIdx})
		tmp=load(processTrackCell{cIdx}(pIdx).outFilePaths_{1});
		tracks=TracksHandle(tmp.tracksFinal);
		condLifetimeCell{pIdx}=[tracks.lifetime];

        measure=[tracks.lifetime];
        endFrames=[tracks.endFrame];
        noncutTracks=([tracks.startFrame]>1)&(endFrames<max(endFrames(:)))
        condUCLifetimeCell{pIdx}=measure(noncutTracks);
        condMaxIntensitiesCell{pIdx}=arrayfun(@(t) max(t.A),tracks(noncutTracks));
        disp(['valid tracks ' num2str(sum(noncutTracks)/length(tracks))]);
    end
    maxLifetime=max(maxLifetime,max([condLifetimeCell{:}]));
    maxIntensity=max(maxIntensity,max(vertcat(condMaxIntensitiesCell{:})));
	lifetimeCell{cIdx}=condLifetimeCell;
    uncutLifetimeCell{cIdx}=condUCLifetimeCell;
    maxIntensitiesCell{cIdx}=condMaxIntensitiesCell;
end
%%
toc
[Handle,~,F]=setupFigure(3,9,27,'AxesWidth',10,'AxesHeight',8,'DisplayMode','print','XSpace',[3 3 3 3],'YSpace',[5 5 5 5]);

%percs=[5 10 25 50 75 90 95];
% percs=[50 55 60 65 70 75 80];
% percs=[50 75 80 85 90];
%percs=linspace(50,95,13);
percs=linspace(1,70,10);

thresholds=prctile(vertcat(maxIntensitiesCell{1}{:}),percs);

c={'r','b','g','y','k'};

histCell=cell(1,length(lifetimeCell));
scoresBin=linspace(1,maxIntensity,100);
for cIdx=1:length(maxIntensitiesCell)
	for pIdx=1:length(maxIntensitiesCell{cIdx})
		measure=(maxIntensitiesCell{cIdx}{pIdx});

		[counts,edges,binIdx]=histcounts(measure,scoresBin);
		counts=counts/sum(counts);
		%[means,meds,stds,orderedIndex,counts] = statPerIndx(lifetimes,binIdx+1);

		histCell{cIdx}=[histCell{cIdx}; counts];
    end
    if(size(histCell{cIdx},1)>1)
       y=histCell{cIdx};
       axes(Handle(1));
       H=shadedErrorBar(scoresBin(1:end-1),(mean(y)),(std(y)),c{cIdx},1);
    else
	    plot(scoresBin(1:end-1),histCell{cIdx},[c{cIdx} '-']);
    end
    vline(thresholds);
end

for tIdx=1:length(thresholds)
    thresholdUncutLifetimeCell=cell(1,length(lifetimeCell));
	for cIdx=1:length(maxIntensitiesCell)
		for pIdx=1:length(maxIntensitiesCell{cIdx})
		    thresholdUncutLifetimeCell{cIdx}{pIdx}=uncutLifetimeCell{cIdx}{pIdx}(maxIntensitiesCell{cIdx}{pIdx}>thresholds(tIdx));
		end
	end
   	displayCondLifetime(thresholdUncutLifetimeCell,names,outputDirPlot,Handle(2*tIdx:2*tIdx+1));
    title(Handle(2*tIdx),[' IMax Thresh: ' num2str(thresholds(tIdx)) ', perc: ' num2str(percs(tIdx))]);
end

function displayCondLifetime(lifetimeCell,names,outputDirPlot,Handles)
c={'r','b','g','y','k'};
scoresBin=3:130;
lifetimeHistCell=cell(1,length(lifetimeCell));
shadedHandles=[];
hold on;
for cIdx=1:length(lifetimeCell)
	for pIdx=1:length(lifetimeCell{cIdx})
		lifetimes=(lifetimeCell{cIdx}{pIdx});
        
		[counts,edges,binIdx]=histcounts(lifetimes,scoresBin);
		counts=counts/sum(counts);
		%[means,meds,stds,orderedIndex,counts] = statPerIndx(lifetimes,binIdx+1);

		lifetimeHistCell{cIdx}=[lifetimeHistCell{cIdx}; counts];
    end
	if(size(lifetimeHistCell{cIdx},1)>1)
        y=lifetimeHistCell{cIdx};
       %% 
       axes(Handles(1));
       H=shadedErrorBar(scoresBin(1:end-1),(mean(y)),(std(y)),c{cIdx},1);
       axes(Handles(2));
       H=shadedErrorBar(scoresBin(1:end-1),log(mean(y)),log(std(y)),c{cIdx},1);
       shadedHandles=[shadedHandles H];
    else
	    plot(scoresBin(1:end-1),lifetimeHistCell{cIdx},[c{cIdx} '-']);
    end
end
lineToLegend=arrayfun(@(h) h.mainLine,shadedHandles,'unif',0);
legend(Handles(1),[lineToLegend{:}],names);
%ylim([0,4])
%xlim([-pi/2+0.05,1.5]);
xlabel('lft (s)');
ylabel('Count')
hold off;
print([outputDirPlot  'lifetime.png'],'-dpng');
print([outputDirPlot  'lifetime.eps'],'-depsc');



