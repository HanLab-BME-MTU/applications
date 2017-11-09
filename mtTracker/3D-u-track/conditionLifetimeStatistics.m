function conditionLifetimeStatistics(processTrackCell,names,outputDirPlot,varargin)
% load scoring daip = inputParser;
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('shade',true);
ip.parse(varargin{:});
p=ip.Results;

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
%% PLot maxIntensity vs Lifetime
[Handle,~,F]=setupFigure(1,1,1,'AxesWidth',10,'AxesHeight',8,'DisplayMode','print','XSpace',[2 2 2 2],'YSpace',[2 2 2 2]);
plotCumulLifetime(Handle,maxIntensitiesCell,uncutLifetimeCell)
printPNGEPSFIG(F,outputDirPlot,'lifetimeVSMaxInt')


%% Plot lft vs intensity for all celll
[Handle,~,F]=setupFigure(length(maxIntensitiesCell),6,6*length(maxIntensitiesCell),'AxesWidth',10,'AxesHeight',8,'DisplayMode','print','XSpace',[3 3 3 3],'YSpace',[5 5 5 5]);
plotCumulLifetimePerCell(Handle,maxIntensitiesCell,uncutLifetimeCell)
printPNGEPSFIG(F,outputDirPlot,'lifetimeVSMaxiPerCell')

%%
toc

%percs=[5 10 25 50 75 90 95];
% percs=[50 55 60 65 70 75 80];
% percs=[50 75 80 85 90];
%percs=linspace(50,95,13);
percs=[0 10 45 50 55 60 70];
nPlot=2+2*length(percs);
[Handle,~,F]=setupFigure(ceil(nPlot/8),8,nPlot,'AxesWidth',30,'AxesHeight',16,'DisplayMode','print','XSpace',[10 10 10 10],'YSpace',[10 10 10 10]);
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

plotCumulLifetime(Handle(2),maxIntensitiesCell,uncutLifetimeCell)
axes(Handle(2));
hline(thresholds);

%%
for tIdx=1:length(thresholds)
    aboveThresholdUncutLifetimeCell=cell(1,length(lifetimeCell));
    belowThresholdUncutLifetimeCell=cell(1,length(lifetimeCell));
	for cIdx=1:length(maxIntensitiesCell)
		for pIdx=1:length(maxIntensitiesCell{cIdx})
		    aboveThresholdUncutLifetimeCell{cIdx}{pIdx}=uncutLifetimeCell{cIdx}{pIdx}(maxIntensitiesCell{cIdx}{pIdx}>thresholds(tIdx));
            belowThresholdUncutLifetimeCell{cIdx}{pIdx}=uncutLifetimeCell{cIdx}{pIdx}(maxIntensitiesCell{cIdx}{pIdx}<=thresholds(tIdx));
		end
	end
   	displayCondLifetime([aboveThresholdUncutLifetimeCell belowThresholdUncutLifetimeCell],[names cellfun(@(n) [n ' < th.'],names,'unif',0)],outputDirPlot,Handle(2*tIdx+1:2*tIdx+2),p.shade);
    title(Handle(2*tIdx+1),[' IMax Th.: ' num2str(thresholds(tIdx)) ', perc: ' num2str(percs(tIdx))]);
end
printPNGEPSFIG(F,outputDirPlot,'lifetime')

%% scaling intensity
scaling=cell(1,length(processTrackCell));
scaledMaxIntensitiesCell=maxIntensitiesCell;
for cIdx=1:length(scaling)
    [scaling{cIdx}, offset, refIdx] = scaleEDFs(maxIntensitiesCell{cIdx}, 'Display', true,...
        'FigureName', ['Ch. ' num2str(cIdx) ' max. intensity scaling']);
    scaledMaxIntensitiesCell{cIdx}=arrayfun(@(i,s) (i{:})*s, scaledMaxIntensitiesCell{cIdx},scaling{cIdx},'unif',0);
end

for tIdx=1:length(thresholds)
    aboveThresholdUncutLifetimeCell=cell(1,length(lifetimeCell));
    belowThresholdUncutLifetimeCell=cell(1,length(lifetimeCell));
    for cIdx=1:length(maxIntensitiesCell)
        for pIdx=1:length(maxIntensitiesCell{cIdx})
            aboveThresholdUncutLifetimeCell{cIdx}{pIdx}=uncutLifetimeCell{cIdx}{pIdx}(scaledMaxIntensitiesCell{cIdx}{pIdx}>thresholds(tIdx));
            belowThresholdUncutLifetimeCell{cIdx}{pIdx}=uncutLifetimeCell{cIdx}{pIdx}(scaledMaxIntensitiesCell{cIdx}{pIdx}<=thresholds(tIdx));
        end
    end
    displayCondLifetime([aboveThresholdUncutLifetimeCell belowThresholdUncutLifetimeCell],[names cellfun(@(n) [n ' < th.'],names,'unif',0)],outputDirPlot,Handle(2*tIdx+1:2*tIdx+2),p.shade);
    title(Handle(2*tIdx+1),[' IMax Th.: ' num2str(thresholds(tIdx)) ', perc: ' num2str(percs(tIdx))]);
end
printPNGEPSFIG(F,outputDirPlot,'scaledLifetime')

function displayCondLifetime(lifetimeCell,names,outputDirPlot,Handles,shade)
c={'r','b','g','y','k'};
scoresBin=3:1:140;
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
	if((size(lifetimeHistCell{cIdx},1)>1))
        y=lifetimeHistCell{cIdx};
        if(shade)
            %%
            axes(Handles(1));
            H=shadedErrorBar(scoresBin(1:end-1),(mean(y)),(std(y)),c{cIdx},1);
            axes(Handles(2));
            H=shadedErrorBar(scoresBin(1:end-1),log(mean(y)),log(std(y)),c{cIdx},1);
            shadedHandles=[shadedHandles H];
        else
            plot(Handles(1),scoresBin(1:end-1),(mean(y)),[c{cIdx} '-']);
            plot(Handles(2),scoresBin(1:end-1),(log(mean(y))),[c{cIdx} '-']);
            shadedHandles=[shadedHandles Handles(2)];
        end
    else
	    plot(scoresBin(1:end-1),lifetimeHistCell{cIdx},[c{cIdx} '-']);
    end
end
if(shade)
    lineToLegend=arrayfun(@(h) h.mainLine,shadedHandles,'unif',0);
    legend(Handles(1),[lineToLegend{:}],names);
else
    legend(Handles(1),names);
end
ylim(Handles(1),[0,0.05]);
xlim(Handles(1),[3,140]);
xlim(Handles(2),[3,140]);

xlabel('lft (s)');
ylabel('Count')
hold off;



function plotCumulLifetimePerCell(Handle,maxIntensitiesCell,uncutLifetimeCell)
hold on
for cIdx=1:length(maxIntensitiesCell)
	for pIdx=1:length(maxIntensitiesCell{cIdx})
		measure=(maxIntensitiesCell{cIdx}{pIdx});
		lifetime=(uncutLifetimeCell{cIdx}{pIdx});
	    scatter(Handle((cIdx-1)*6+pIdx),lifetime,measure);
    end
end
xlabel('lifetime (f)')
ylabel('Max Intensity (A.U.)')
hold off

function plotCumulLifetime(Handle,maxIntensitiesCell,uncutLifetimeCell)
hold on
for cIdx=1:length(maxIntensitiesCell)
	for pIdx=1:length(maxIntensitiesCell{cIdx})
		measure=(maxIntensitiesCell{cIdx}{pIdx});
		lifetime=(uncutLifetimeCell{cIdx}{pIdx});
	    scatter(Handle,lifetime,measure);
    end
end
xlabel(Handle,'lifetime (f)')
ylabel(Handle,'Max Intensity (A.U.)')
hold off
