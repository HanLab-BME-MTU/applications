function lifetimeThresholdDebugging(processTrackCell,names,outputDirPlot,varargin)
% load scoring daip = inputParser;
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('shade',true);
ip.addParameter('render',true);
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
uncutTrackIDCell=cell(1,length(processTrackCell));
uncutTrackCell=cell(1,length(processTrackCell));

%% detecting cut tracks.
tic;
for cIdx=1:length(processTrackCell)
	condLifetimeCell=cell(1,length(processTrackCell{cIdx}));
 	condUCLifetimeCell=cell(1,length(processTrackCell{cIdx}));
 	condMaxIntensitiesCell=cell(1,length(processTrackCell{cIdx}));
 	condTrackIDCell=cell(1,length(processTrackCell{cIdx}));
 	condTrackCell=cell(1,length(processTrackCell{cIdx}));

	parfor pIdx=1:length(processTrackCell{cIdx})
		tmp=load(processTrackCell{cIdx}(pIdx).outFilePaths_{1});
		tracks=TracksHandle(tmp.tracksFinal);
		condLifetimeCell{pIdx}=[tracks.lifetime];
        condTrackIDCell{pIdx}=1:length(tracks);
        
        measure=[tracks.lifetime];
        endFrames=[tracks.endFrame];
        noncutTracks=([tracks.startFrame]>1)&(endFrames<max(endFrames(:)));
        condUCLifetimeCell{pIdx}=measure(noncutTracks);
        
        condTrackCell{pIdx}=tracks(noncutTracks);
        condTrackIDCell{pIdx}= condTrackIDCell{pIdx}(noncutTracks);
        condMaxIntensitiesCell{pIdx}=arrayfun(@(t) max(t.A),tracks(noncutTracks));
        disp(['valid tracks ' num2str(sum(noncutTracks)/length(tracks))]);
    end
    maxLifetime=max(maxLifetime,max([condLifetimeCell{:}]));
    maxIntensity=max(maxIntensity,max(vertcat(condMaxIntensitiesCell{:})));
	lifetimeCell{cIdx}=condLifetimeCell;
    uncutLifetimeCell{cIdx}=condUCLifetimeCell;
    maxIntensitiesCell{cIdx}=condMaxIntensitiesCell;
    uncutTrackIDCell{cIdx}=condTrackIDCell;
    uncutTrackCell{cIdx}=condTrackCell;
end

%% PLot maxIntensity vs Lifetime and show selected tracks
minLft=0;
maxLft=30;
minInt=50;
maxInt=maxIntensity;

[Handle,~,F]=setupFigure(1,1,1,'AxesWidth',10,'AxesHeight',8,'DisplayMode','print','XSpace',[2 2 2 2],'YSpace',[2 2 2 2]);
plotCumulLifetime(Handle,maxIntensitiesCell,uncutLifetimeCell)
hold on;
plot(Handle,[minLft maxLft],[minInt maxInt],'-r','linewidth',2);
hold off; 
printPNGEPSFIG(F,outputDirPlot,'lifetimeSelected')

%% select short and bright tracks
cutFunction=(@(i,l) (i-minInt) > (l*(maxInt-minInt)/(maxLft-minLft)));
brightShortTracks=uncutTrackCell;
brightShortTracksIDX=uncutTrackIDCell;

for cIdx=1:length(maxIntensitiesCell)
	for pIdx=1:length(maxIntensitiesCell{cIdx})
		iMax=(maxIntensitiesCell{cIdx}{pIdx});
        lft=(uncutLifetimeCell{cIdx}{pIdx});
        brightShortTracks{cIdx}{pIdx}=brightShortTracks{cIdx}{pIdx}(cutFunction(iMax',lft));
        brightShortTracksIDX{cIdx}{pIdx}=(cutFunction(iMax',lft));        
    end
end

%%
toc
percs=linspace(0,50,8);
%percs=[0 10 15 20 45 50];
nPlot=2+1*length(percs);

%% Display bona fide CCP vs other struct using normal intensity on problematic tracks
% [Handle,~,F]=setupFigure(ceil(nPlot/8),8,nPlot,'AxesWidth',30,'AxesHeight',16,'DisplayMode','print','XSpace',[10 10 10 10],'YSpace',[10 10 10 10],'Name','lifetime-suppressBrightAndShort');
% thresholds=prctile(vertcat(maxIntensitiesCell{1}{:}),percs);
% c={'r','b','g','y','k'};
% for tIdx=1:length(thresholds)
%     aboveThresholdUncutLifetimeCell=cell(1,length(lifetimeCell));
%     belowThresholdUncutLifetimeCell=cell(1,length(lifetimeCell));
% 
% 	for cIdx=1:length(maxIntensitiesCell)
% 		for pIdx=1:length(maxIntensitiesCell{cIdx})
%             nonSelectedTrack=~brightShortTracksIDX{cIdx}{pIdx}';
% 		    aboveThresholdUncutLifetimeCell{cIdx}{pIdx}=uncutLifetimeCell{cIdx}{pIdx}((maxIntensitiesCell{cIdx}{pIdx}>thresholds(tIdx))&nonSelectedTrack);
%             belowThresholdUncutLifetimeCell{cIdx}{pIdx}=uncutLifetimeCell{cIdx}{pIdx}((maxIntensitiesCell{cIdx}{pIdx}<=thresholds(tIdx))&nonSelectedTrack);
% 		end
% 	end
%    	displayCondLifetime([aboveThresholdUncutLifetimeCell belowThresholdUncutLifetimeCell],[names cellfun(@(n) [n ' < th.'],names,'unif',0)],outputDirPlot,Handle(tIdx+2),p.shade,tIdx==1);
%     title(Handle(tIdx+2),[' IMax Th.: ' num2str(thresholds(tIdx)) ', perc: ' num2str(percs(tIdx))]);
% end
% printPNGEPSFIG(F,outputDirPlot,'lifetime-suppressBrightAndShort')

%% scaling intensity
scaling=cell(1,length(processTrackCell));
scaledMaxIntensitiesCell=maxIntensitiesCell;
for cIdx=1:length(scaling)
    [scaling{cIdx}, offset, refIdx] = scaleEDFs(maxIntensitiesCell{cIdx}, 'Display', true,...
        'FigureName', ['Ch. ' num2str(cIdx) ' max. intensity scaling']);
    scaledMaxIntensitiesCell{cIdx}=arrayfun(@(i,s) (i{:})*s, scaledMaxIntensitiesCell{cIdx},scaling{cIdx},'unif',0);
end

P=28;
T=prctile(vertcat(scaledMaxIntensitiesCell{1}{:}),P);
aboveThresholdUncutLifetimeCell=cell(1,length(lifetimeCell));
belowThresholdUncutLifetimeCell=cell(1,length(lifetimeCell));
for cIdx=1:length(scaledMaxIntensitiesCell)
    for pIdx=1:length(scaledMaxIntensitiesCell{cIdx})
        aboveThreshold=(scaledMaxIntensitiesCell{cIdx}{pIdx}>T);
        aboveThresholdUncutLifetimeCell{cIdx}{pIdx}=uncutLifetimeCell{cIdx}{pIdx}(aboveThreshold);
        belowThresholdUncutLifetimeCell{cIdx}{pIdx}=uncutLifetimeCell{cIdx}{pIdx}(~aboveThreshold);
    end
end
[Handle,~,F]=setupFigure(1,2,2,'AspectRatio',1,'AxesWidth',10,'DisplayMode','print','Name','lifetime-scaled-IMax (A.U.)/%: ','DisplayMode','screen');
plotCumulLifetime(Handle(1),maxIntensitiesCell,uncutLifetimeCell)
axes(Handle(1));
hline(T,'-');
displayCondLifetime([aboveThresholdUncutLifetimeCell belowThresholdUncutLifetimeCell],[names cellfun(@(n) [n ' < th.'],names,'unif',0)],outputDirPlot,Handle(2),p.shade,true);
printPNGEPSFIG(F,outputDirPlot,['lifetimeScaled-P' num2str(P)])


%% Display bona fide CCP vs other struct using scaled intensity.
[Handle,~,F]=setupFigure(ceil(nPlot/10),10,nPlot,'AspectRatio',1,'AxesWidth',50,'AxesHeight',20,'DisplayMode','print','XSpace',[40 40 40 40],'YSpace',[20 10 20 10],'Name',['lifetimeScaled-P' num2str(P)],'DisplayMode','screen');
%[Handle,~,F]=setupFigure(ceil(nPlot/10),10,nPlot,'AspectRatio',1, 'AxesWidth',30,'DisplayMode','multiple','Name','lifetime-scaled','DisplayMode','screen');
thresholds=prctile(vertcat(scaledMaxIntensitiesCell{1}{:}),percs);

plotCumulLifetime(Handle(2),maxIntensitiesCell,uncutLifetimeCell)
axes(Handle(2));
hline(thresholds,'-');

c={'r','b','g','y','k'};
for tIdx=1:length(thresholds)
    aboveThresholdUncutLifetimeCell=cell(1,length(lifetimeCell));
    belowThresholdUncutLifetimeCell=cell(1,length(lifetimeCell));

	for cIdx=1:length(scaledMaxIntensitiesCell)
		for pIdx=1:length(scaledMaxIntensitiesCell{cIdx})
            aboveThreshold=(scaledMaxIntensitiesCell{cIdx}{pIdx}>thresholds(tIdx));
		    aboveThresholdUncutLifetimeCell{cIdx}{pIdx}=uncutLifetimeCell{cIdx}{pIdx}(aboveThreshold);
            belowThresholdUncutLifetimeCell{cIdx}{pIdx}=uncutLifetimeCell{cIdx}{pIdx}(~aboveThreshold);
		end
	end
   	displayCondLifetime([aboveThresholdUncutLifetimeCell belowThresholdUncutLifetimeCell],[names cellfun(@(n) [n ' < th.'],names,'unif',0)],outputDirPlot,Handle(2+tIdx),p.shade,tIdx==1);
    title(Handle(tIdx+2),[num2str( ceil(thresholds(tIdx))) '/' num2str(ceil(percs(tIdx)))]);
end
printPNGEPSFIG(F,outputDirPlot,'lifetimeScaled')


%% render below and above threshold
if(p.render)
    tic;
    bonafideCCP=uncutTrackCell;
    otherStructure=uncutTrackCell;
    for cIdx=1:length(maxIntensitiesCell)
        processOverlayCell=cell(1,length(maxIntensitiesCell{cIdx}));
        for pIdx=1:length(maxIntensitiesCell{cIdx})
            aboveThreshold=(scaledMaxIntensitiesCell{cIdx}{pIdx}>prctile(vertcat(scaledMaxIntensitiesCell{1}{:}),50));
            bonafideCCP{cIdx}{pIdx}=bonafideCCP{cIdx}{pIdx}(aboveThreshold);
            otherStructure{cIdx}{pIdx}=otherStructure{cIdx}{pIdx}(~aboveThreshold);
            MDCrop=processTrackCell{cIdx}(pIdx).getOwner();
            processOverlay=ExternalProcess(MDCrop,'overlayProjTracksMovie');
            overlayProjTracksMovie(MDCrop.getPackage(100).getProcess(4),'tracksOrProcess',bonafideCCP{cIdx}{pIdx},'process',processOverlay,'name','bonafideCCP','colormap',[255 0 0]);
            overlayProjTracksMovie(processOverlay,'tracksOrProcess',otherStructure{cIdx}{pIdx},'process',processOverlay,'name','bonafideCCP-otherStructure','colormap',[0 0 255]);
            processOverlayCell{pIdx}=processOverlay;
        end
        printProcMIPArray(processOverlayCell,[outputDirPlot filesep 'MIPArray']) 

    end
    disp('Rendering time') 
    toc;
end


function displayCondLifetime(lifetimeCell,names,outputDirPlot,Handles,shade,plotLegend)
c={'r','b','g','y','k'};
scoresBin=3:1:140;
lifetimeHistCell=cell(1,length(lifetimeCell));
shadedHandles=[];
hold(Handles(1),'on');
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
            %axes(Handles(2));
            %H=shadedErrorBar(scoresBin(1:end-1),log(mean(y)),log(std(y)),c{cIdx},1);
            shadedHandles=[shadedHandles H];
        else
            plot(Handles(1),scoresBin(1:end-1),(mean(y)),[c{cIdx} '-']);
            %plot(Handles(2),scoresBin(1:end-1),(log(mean(y))),[c{cIdx} '-']);
            shadedHandles=[shadedHandles Handles(1)];
        end
    else
	    plot(scoresBin(1:end-1),lifetimeHistCell{cIdx},[c{cIdx} '-']);
    end
end
hold(Handles(1),'off');
if(plotLegend)
if(shade)
    lineToLegend=arrayfun(@(h) h.mainLine,shadedHandles,'unif',0);
    legend(Handles(1),[lineToLegend{:}],names);
else
    legend(Handles(1),names,'Location','northeast');
end
end
ylim(Handles(1),[0,0.05]);
xlim(Handles(1),[3,140]);
%xlim(Handles(2),[3,140]);

xlabel(Handles(1),'lft (s)');
ylabel(Handles(1),'Count');



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
ylabel('Max Int. (A.U.)')
hold off

function printPNGEPSFIG(fhandle,outputDirPlot,filename)
set(0,'CurrentFigure',fhandle);
print([outputDirPlot filesep filename '.png'],'-dpng');
print([outputDirPlot filesep filename '.eps'],'-depsc');
saveas(gcf,[outputDirPlot filesep filename '.fig'])

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
ylabel(Handle,'Max Int. (A.U.)')
hold off
