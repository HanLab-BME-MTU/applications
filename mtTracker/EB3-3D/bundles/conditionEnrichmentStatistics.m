function conditionEnrichmentStatistics(packageList,names,outputDirPlot,varargin)
% packageList describe a condition (either experimental or algorithmical).

ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('shade',true);
ip.addParameter('printCellName',true);
ip.addParameter('projProcesses',[]);
ip.parse(varargin{:});
p=ip.Results;


projProcessesLabFrame=(arrayfun(@(pack) pack.getProcess(4),packageList,'unif',0))';

enrichmentProcesses=(arrayfun(@(pack) pack.getProcess(6),packageList,'unif',0))';
projProcessesAstral=(arrayfun(@(pack) pack.getProcess(7),packageList,'unif',0))';
projProcessesPolar=(arrayfun(@(pack) pack.getProcess(8),packageList,'unif',0))';

projProcessesAstralMovie=(arrayfun(@(pack) pack.getProcess(9),packageList,'unif',0))';
projProcessesPolarMovie=(arrayfun(@(pack) pack.getProcess(10),packageList,'unif',0))';

%% Sumarizing movies
% nProcess=length(enrichmentProcesses);
% processRenderArray=printProcMIPArray([ num2cell([[projProcessesPolarMovie{:}]])], ... 
%         [outputDirPlot filesep 'PolarArrayMovie'],'MIPIndex',[3*ones(1,2*nProcess)],'maxHeight',1200, ... 
%         'MIPSize','auto');

%% A clusterfuck of data collection below...
if(~iscell(enrichmentProcesses))
    enrichmentProcesses={enrichmentProcesses};
end

mkClrDir(outputDirPlot);

% load score set in cells of cells
elevationCell=cell(1,length(enrichmentProcesses));
densityCell=cell(1,length(enrichmentProcesses));
distanceCell=cell(1,length(enrichmentProcesses));
ampCell=cell(1,length(enrichmentProcesses));
timeCell=cell(1,length(enrichmentProcesses));

% Define astral and polar region in a naive way
astralThresh=0;
polarThresh=0;
distBinning=[];

%% Cumulative plots rendering movie per movie
% printProjProcessArray(projProcesses,'OrthoProj',outputDirPlot,{'All'}, ... 
%             'splitMovie',true,'forceSize',[],'outFile',4);
% printProcessArray({{[projProcesses{:}]} {[projProcessesPolar{:}]}},{[enrichmentProcesses{:}]},3,outputDirPlot,names(1));
% printProcessArray({{[projProcesses{:}]} {[projProcessesPolar{:}]}},{[enrichmentProcesses{:}]},4,outputDirPlot,{[names{:}]});
% printProcessArray({projProcesses projProcessesPolar},enrichmentProcesses,5,outputDirPlot,names);
% printProcessArray({projProcesses projProcessesPolar},enrichmentProcesses,6,outputDirPlot,names);
% printProcessArray({projProcesses projProcessesPolar},enrichmentProcesses,7,outputDirPlot,names);

nonEmptyProcess=find(cellfun(@(p) ~isempty(p),enrichmentProcesses))';


%% Condition structure for easier interpretation
for cIdx=1:length(enrichmentProcesses)
    condition.name=names{cIdx};
    condElevationCell=cell(1,length(enrichmentProcesses{cIdx}));
    condDistanceCell=cell(1,length(enrichmentProcesses{cIdx}));
	condDensityCell=cell(1,length(enrichmentProcesses{cIdx}));
    condAmpCell=cell(1,length(enrichmentProcesses{cIdx}));
    condTimeCell=cell(1,length(enrichmentProcesses{cIdx}));
	for pIdx=1:length(enrichmentProcesses{cIdx})
		tmp=load(enrichmentProcesses{cIdx}(pIdx).outFilePaths_{2});
		condElevationCell{pIdx}=tmp.elevations;
		condDensityCell{pIdx}=tmp.densities;
        condDistanceCell{pIdx}=tmp.poleDistances;
        polarDistBinning=tmp.polarDistBinning;
        astralDistBinning=tmp.astralDistBinning;
        condAmpCell{pIdx}=tmp.amps;
        condTimeCell{pIdx}=tmp.timesSteps;
        astralThresh=tmp.astralThresh;
        polarThresh=tmp.polarThresh;
	end
	elevationCell{cIdx}=condElevationCell;
	densityCell{cIdx}=condDensityCell;
    distanceCell{cIdx}=condDistanceCell;
    ampCell{cIdx}=condAmpCell;
    timeCell{cIdx}=condTimeCell;
end



% Count and density function of distance to poles
% Overlay decay 
AstralK=1/2;
PolarK=0;

%distBinning=linspace(30,70,20);
coneVolumes=2*pi*(1-cos(pi/2-polarThresh))*polarDistBinning.^3/3;
polarVolumes=(coneVolumes(2:end)-coneVolumes(1:(end-1)));
coneVolumes=2*pi*(1-cos(pi/2+astralThresh))*astralDistBinning.^3/3;
astralVolumes=(coneVolumes(2:end)-coneVolumes(1:(end-1)));

astralXTicks=astralDistBinning(1:end-1)*0.1;
polarXTicks=polarDistBinning(1:end-1)*0.1;

relativeVolumeAstral=astralVolumes/astralVolumes(1);
relativeVolumePolar=polarVolumes/polarVolumes(1);

meansPolarDensitiesCell=cell(1,length(enrichmentProcesses));
meansAstralDensitiesCell=cell(1,length(enrichmentProcesses));
meansPolarIntCell=cell(1,length(enrichmentProcesses));
meansAstralIntCell=cell(1,length(enrichmentProcesses));
coneNormCountPolarCell=cell(1,length(enrichmentProcesses));
coneNormCountAstralCell=cell(1,length(enrichmentProcesses));
countPolarCell=cell(1,length(enrichmentProcesses));
countAstralCell=cell(1,length(enrichmentProcesses));
normAmpPolarCell=cell(1,length(enrichmentProcesses));
normAmpAstralCell=cell(1,length(enrichmentProcesses));
cmap=winter(2*length(enrichmentProcesses));
cmapAstral=winter(length(enrichmentProcesses));
cmapPolar=autumn(length(enrichmentProcesses));

for cIdx=nonEmptyProcess
    conditionEnrichmenData(cIdx).name=names{cIdx};
    moviesStruct=[];
    for pIdx=1:length(enrichmentProcesses{cIdx})
        allDist=vertcat(distanceCell{cIdx}{pIdx}{:});
        allDens=vertcat(densityCell{cIdx}{pIdx}{:});
        allElevs=vertcat(elevationCell{cIdx}{pIdx}{:});
        allAmp=vertcat(ampCell{cIdx}{pIdx}{:});

        tmp=load(enrichmentProcesses{cIdx}(pIdx).outFilePaths_{2});
        allAzim=vertcat(tmp.azimuths{:});
        allKTAzim=vertcat(tmp.KTazimuths{:});
        allKTElev=vertcat(tmp.KTelevations{:});
        allKTDist=vertcat(tmp.KTpoleDistances{:});
        allKTTime=vertcat(tmp.KTTimesSteps{:});
        polarDetectionsMask=tmp.polarDetectionsMask;
        astralDetectionsMask=tmp.astralDetectionsMask;

        % separate astral and polar and sort function of distances
        [polarCounts,polarEdges,polarBinIdx]=histcounts(allDist(polarDetectionsMask),polarDistBinning);
        [meanPolarDensities,medPolarDensities,stdPolarDensities]=statPerIndx(allDens(polarDetectionsMask),polarBinIdx,1:length(polarDistBinning)-1);
        [meanPolarAmp,medPolarAmp,stdPolarAmp,~,polarCount,sumPolarAmp]=statPerIndx(allAmp(polarDetectionsMask),polarBinIdx,1:length(polarDistBinning)-1);
        coneVolumes=2*pi*(1-cos(pi/2-polarThresh))*polarDistBinning.^3/3;
        polarVolumes=(coneVolumes(2:end)-coneVolumes(1:(end-1)));
        normalizedPolarAmp=sumPolarAmp./polarVolumes;
        normalizedPolarCount=polarCounts./polarVolumes;

        [astralCounts,astralEdges,astralBinIdx]=histcounts(allDist(astralDetectionsMask),astralDistBinning);
        [meanAstralDensities,medAstralDensities,stdAstralDensities]=statPerIndx(allDens(astralDetectionsMask),astralBinIdx,1:length(astralDistBinning)-1);
        [meanAstralAmp,medAstralAmp,stdAstralAmp,~,astralCount,sumAstralAmp]=statPerIndx(allAmp(astralDetectionsMask),astralBinIdx,1:length(astralDistBinning)-1);
        coneVolumes=2*pi*(1-cos(pi/2+astralThresh))*astralDistBinning.^3/3;
        astralVolumes=(coneVolumes(2:end)-coneVolumes(1:(end-1)));
        normalizedAstralAmp=sumAstralAmp./astralVolumes;
        normalizedAstralCount=astralCounts./astralVolumes;

        MD=enrichmentProcesses{cIdx}(pIdx).getOwner();


        % no scaling at first
        % allDens=allDens/mean(allDens(allElevs<0));  

        meansPolarDensitiesCell{cIdx}=[meansPolarDensitiesCell{cIdx}; meanPolarDensities];
        meansAstralDensitiesCell{cIdx}=[meansAstralDensitiesCell{cIdx}; meanAstralDensities];
        coneNormCountPolarCell{cIdx}=[coneNormCountPolarCell{cIdx}; normalizedPolarCount];
        coneNormCountAstralCell{cIdx}=[coneNormCountAstralCell{cIdx}; normalizedAstralCount];
        countPolarCell{cIdx}=[countPolarCell{cIdx}; polarCounts];
        countAstralCell{cIdx}=[countAstralCell{cIdx}; astralCounts];
        normAmpAstralCell{cIdx}=[normAmpAstralCell{cIdx}; normalizedAstralAmp];
        normAmpPolarCell{cIdx}=[normAmpPolarCell{cIdx}; normalizedPolarAmp];
        meansAstralIntCell{cIdx}=[meansAstralIntCell{cIdx}; meanAstralAmp];
        meansPolarIntCell{cIdx}=[meansPolarIntCell{cIdx}; meanPolarAmp];

        movie.movieData=enrichmentProcesses{cIdx}(pIdx).getOwner().movieDataPath_;
        movie.poleDist=vertcat(distanceCell{cIdx}{pIdx}{:}); 
        movie.binningDistance=polarDistBinning; 
        movie.binnedCountAstral=astralCounts; 
        movie.binnedCountPolar=polarCounts; 
        movie.binnedCountAstralPerFrame=astralCounts/MD.nFrames_; 
        movie.binnedCountPolarPerFrame=polarCounts/MD.nFrames_;

        KTPerFrame=mean(cellfun(@(k) length(k),tmp.KTTimesSteps));

        movie.binnedCountAstralPerFramePerKT=(astralCounts/KTPerFrame)/MD.nFrames_; 
        movie.binnedCountPolarPerFramePerKT=(polarCounts/KTPerFrame)/MD.nFrames_;  
        
        [KTDistCounts,edges,index]=histcounts(allKTDist,[polarDistBinning]);

        KTDistCounts=KTDistCounts./MD.nFrames_;
        KTROICounts=fliplr(cumsum(fliplr(KTDistCounts)));
        movie.binnedCountAstralPerFramePerKTAtDist=(astralCounts./KTROICounts)/MD.nFrames_; 
        movie.binnedCountPolarPerFramePerKTAtDist=(polarCounts./KTROICounts)/MD.nFrames_;  
        


        astalAt2UM=mean(astralCounts((astralDistBinning>20)&(astralDistBinning<26)));
        polarAt2UM=mean(polarCounts((polarDistBinning>20)&(polarDistBinning<26)));

        % astalAt2UM=mean(astralCounts((astralDistBinning==20)));
        % polarAt2UM=mean(polarCounts((polarDistBinning==20)));

        movie.binnedCountAstralNorm2um=astralCounts/astalAt2UM; 
        movie.binnedCountPolarNorm2um=polarCounts/polarAt2UM;  

        KTCountWeigth=KTROICounts./KTROICounts(1);
        movie.binnedCountAstralNorm2umWeighted=astralCounts./(astalAt2UM); 
        movie.binnedCountPolarNorm2umWeighted=polarCounts./(polarAt2UM*KTCountWeigth);  

        movie.elevation=vertcat(elevationCell{cIdx}{pIdx}{:});
        movie.intensity=vertcat(ampCell{cIdx}{pIdx}{:});
        movie.density=vertcat(densityCell{cIdx}{pIdx}{:});;
        movie.azimuths=allAzim;
        movie.KTazimuths=allKTAzim;
        movie.KTelevations=allKTElev;
        movie.KTdist=allKTDist;
        movie.KTROICounts=KTROICounts;
        movie.KTDistCounts=KTDistCounts;
        movie.KTtime=allKTTime;
        movie.time=vertcat(timeCell{cIdx}{pIdx}{:});
        moviesStruct=[moviesStruct movie];
    end
    conditionEnrichmenData(cIdx).movies=moviesStruct;    
end

%% Retrieving projection process
projProcessesLabFrame=(arrayfun(@(pack) pack.getProcess(4),packageList,'unif',0))';
projProcessesPolarMovie=(arrayfun(@(pack) pack.getProcess(10),packageList,'unif',0))';
%% Plotting analysis data
processGraphCell=cell(length(conditionEnrichmenData),4);
for cIdx=1:length(conditionEnrichmenData)
    processGraphCell{cIdx,1}=plotCountVsDistance(conditionEnrichmenData,astralXTicks,polarXTicks,names,outputDirPlot,'type','per frame','lightCond',cIdx);
    processGraphCell{cIdx,2}=plotCountVsDistance(conditionEnrichmenData,astralXTicks,polarXTicks,names,outputDirPlot,'type','normalized 2um','lightCond',cIdx);
    processGraphCell{cIdx,3}=plotCountVsDistance(conditionEnrichmenData,astralXTicks,polarXTicks,names,outputDirPlot,'type','KT ROI Count','lightCond',cIdx);
    processGraphCell{cIdx,4}=plotCountVsDistance(conditionEnrichmenData,astralXTicks,polarXTicks,names,outputDirPlot,'type','weighted','lightCond',cIdx);
end
nProcess=length(enrichmentProcesses);
processRenderArray=printProcMIPArray( [projProcessesLabFrame projProcessesAstral projProcessesPolar processGraphCell(:,1)  processGraphCell(:,2) projProcessesPolarMovie processGraphCell(:,3) processGraphCell(:,4)], ... 
        [outputDirPlot filesep 'GraphFullPolarArrayMovie'],'MIPIndex',[1*ones(1,nProcess) 3*ones(1,nProcess)  3*ones(1,nProcess) ones(1,nProcess)  ones(1,nProcess) 3*ones(1,nProcess) ones(1,nProcess)  ones(1,nProcess)],'maxHeight',1200, ... 
        'MIPSize','auto','printFrame',1:50);
figure;imshow(imread(sprintfPath(processRenderArray.outFilePaths_{1},1)));
return;

plotCountVsDistance(conditionEnrichmenData,astralXTicks,polarXTicks,names,outputDirPlot,'type','per frame');
plotCountVsDistance(conditionEnrichmenData,astralXTicks,polarXTicks,names,outputDirPlot,'type','normalized 2um');
plotCountVsDistance(conditionEnrichmenData,astralXTicks,polarXTicks,names,outputDirPlot,'type','weighted');

%% For each condition
processGraphCell=cell(length(conditionEnrichmenData),4);
for cIdx=1:length(conditionEnrichmenData)
    processGraphCell{cIdx,1}=plotCountVsDistance(conditionEnrichmenData,astralXTicks,polarXTicks,names,outputDirPlot,'type','per frame','lightCond',cIdx);
    processGraphCell{cIdx,2}=plotCountVsDistance(conditionEnrichmenData,astralXTicks,polarXTicks,names,outputDirPlot,'type','normalized 2um','lightCond',cIdx);
    processGraphCell{cIdx,3}=plotCountVsDistance(conditionEnrichmenData,astralXTicks,polarXTicks,names,outputDirPlot,'type','KT ROI Count','lightCond',cIdx);
    % processGraphCell{cIdx,4}=plotCountVsDistance(conditionEnrichmenData,astralXTicks,polarXTicks,names,outputDirPlot,'type','perFramePerKT','lightCond',cIdx);
    processGraphCell{cIdx,4}=plotCountVsDistance(conditionEnrichmenData,astralXTicks,polarXTicks,names,outputDirPlot,'type','weighted','lightCond',cIdx);
end

%% drawing each condition in a split fashion
for cIdx=1:length(conditionEnrichmenData)
    processRenderArray=printProcMIPArray( {processGraphCell{cIdx,:} projProcessesAstral{cIdx} projProcessesPolar{cIdx}}, ... 
            [outputDirPlot filesep 'GraphPolarArray-Split'],'MIPIndex',[ones(1,size(processGraphCell,2)) 3*ones(1,2)],'maxHeight',1200, ... 
            'MIPSize','auto','MIParrangement','auto');
    movefile(sprintfPath(processRenderArray.outFilePaths_{1},1),[outputDirPlot filesep 'GraphPolarArray-Split' filesep 'GraphPolarArray-' num2str(cIdx) '.png']);
end

%% Show cell growth compared with cell
nProcess=length(enrichmentProcesses);
processRenderArray=printProcMIPArray( [processGraphCell num2cell([[projProcessesAstral{:}]' [projProcessesPolar{:}]'])], ... 
        [outputDirPlot filesep 'GraphAstralPolarArray'],'MIPIndex',[ones(1,size(processGraphCell,2)*nProcess) 3*ones(1,2*nProcess)],'maxHeight',1200, ... 
        'MIPSize','auto','MIParrangement','auto');

nProcess=length(enrichmentProcesses);
processRenderArray=printProcMIPArray( [processGraphCell num2cell([[projProcessesAstralMovie{:}]' [projProcessesPolarMovie{:}]'])], ... 
        [outputDirPlot filesep 'GraphAstralPolarArrayMovie'],'MIPIndex',[ones(1,size(processGraphCell,2)*nProcess) 3*ones(1,2*nProcess)],'maxHeight',1200, ... 
        'MIPSize','auto');





return;

%% Display elevation vs density data for MT only 
cmap=cool(length(enrichmentProcesses));
c={'r','b','g','y','k','o'};
scoresBin=-pi/2:0.1:pi/2;
meansDensitiesCell=cell(1,length(enrichmentProcesses));
Handles=[];
[Handle,~,F]=setupFigure(1,1,1,'AxesWidth',10,'AxesHeight',6);
hold on;
for cIdx=nonEmptyProcess
    for pIdx=1:length(enrichmentProcesses{cIdx})
        allElevs=vertcat(elevationCell{cIdx}{pIdx}{:});
        allDens=vertcat(densityCell{cIdx}{pIdx}{:});
        
        % scale astral density to one
        allDens=allDens/mean(allDens(allElevs<0));
        
        [counts,edges,binIdx]=histcounts(allElevs,scoresBin);
        [means,meds,stds,orderedIndex,counts] = statPerIndx(allDens,binIdx+1);

        meansDensitiesCell{cIdx}=[meansDensitiesCell{cIdx}; means];
    end
    if(~isempty((meansDensitiesCell{cIdx})))
        if(size(meansDensitiesCell{cIdx},1)>1)
            y=meansDensitiesCell{cIdx};
            if(p.shade)
                H=shadedErrorBar(scoresBin(1:end),mean(y),std(y),{'Color',cmap(cIdx,:)});
                Handles=[Handles H];
            else
                plot(scoresBin(1:end),mean(y),[c{1} '-'],'Color',cmap(cIdx,:));
            end
        else
            plot(scoresBin(1:end),meansDensitiesCell{cIdx},[c{1} '-'],'Color',cmap(cIdx,:));
        end
    end
end
lineToLegend=arrayfun(@(h) h.mainLine,Handles,'unif',0)
legend([lineToLegend{:}],names(nonEmptyProcess),'Location','northwest')
ylim([0,5])
xlim([-pi/2+0.05,1.5]);
xlabel('Elevation (rad)');
ylabel('Comet density (1/\mu m^3)')
hold off;
printPNGEPSFIG(F,outputDirPlot,'enrichment');
close(F);

%% Mean intensity vs pole distance.
[Handle,~,F]=setupFigure(1,1,1,'AxesWidth',10,'AxesHeight',6);
hold on;
shadeHandles=[];
predHandles=[];
minDist=2;  
fitXTicks=astralXTicks(astralXTicks>=minDist);
astralK=0; 
for cIdx=nonEmptyProcess

    %X = linsolve([exp(-XTicks*AstralK)'],mean(meansAstralIntCell{cIdx})');
    [f,B]=catastropheRateFit(mean(meansAstralIntCell{cIdx}(:,astralXTicks>=minDist)),fitXTicks);
    predictedAstral= f(B,astralXTicks);
    astralK=B(2); 

    % A=sum(mean(meansAstralIntCell{cIdx}))./sum(exp(-astralXTicks*AstralK));
    % predictedAstral= A*exp(-astralXTicks*AstralK);

    %X = linsolve([exp(-XTicks*AstralK)'],mean(meansPolarIntCell{cIdx})');
    A=sum(mean(meansPolarIntCell{cIdx}))./sum(exp(-polarXTicks*AstralK));
    predictedPolar= A*exp(-polarXTicks*AstralK);    
    
    if(~isempty((meansAstralIntCell{cIdx})))
        shadeHandle=plotOrShade(Handle(1),astralXTicks,meansAstralIntCell{cIdx},p.shade,cmapAstral(cIdx,:));
        shadeHandles=[shadeHandles shadeHandle];
        % pH=plot(Handle(1),astralXTicks,predictedAstral);
        % predHandles=[predHandles pH];
    end    
    
    if(~isempty((meansPolarIntCell{cIdx})))
        shadeHandle=plotOrShade(Handle(1),polarXTicks,meansPolarIntCell{cIdx},p.shade,cmapPolar(cIdx,:));
        shadeHandles=[shadeHandles shadeHandle];
        % pH=plot(Handle(1),polarXTicks,predictedPolar);
        % predHandles=[predHandles pH];
    end
end
lineToLegend=arrayfun(@(h) h.mainLine,shadeHandles,'unif',0);
legend([lineToLegend{:} predHandles],{'Astral','Polar',['predAstral K: ' num2str(astralK)]},'Location','northeast');
% ylim([0,4])
xlim([0 7]);
xlabel('Pole distance (\mu m)');
ylabel(Handle,'Mean Intensities');
hold off;
printPNGEPSFIG(F,outputDirPlot,'meanIntensities');
%close(F);


%% Proximal and distal counts
[Handle,~,F]=setupFigure(1,1,1,'AxesWidth',10,'AxesHeight',6);
set(F,'Visible', 'off'); 
shadeHandles=[];
predHandles=[];
for cIdx=nonEmptyProcess
        lcountAstralCell=countAstralCell{cIdx};
        lcountPolarCell=countPolarCell{cIdx};        

        decayAmps=sum((lcountAstralCell),2)./sum(exp(-astralXTicks*AstralK));
        predictedAstral= decayAmps*exp(-astralXTicks*AstralK);
        lowDistances=(astralXTicks<=median(astralXTicks));
        highDistances=(astralXTicks>=median(astralXTicks));
        
        if(~isempty((lcountAstralCell)))
            notBoxPlot([sum(lcountAstralCell(:,lowDistances),2) sum(lcountAstralCell(:,highDistances),2) sum(predictedAstral(:,highDistances),2)]);
        end
end
ylabel('Astral Comet counts');
Handle.XTickLabel={'Proximal','Distal','Predicted'};
printPNGEPSFIG(F,outputDirPlot,'cometCountBoxAstral');

%% Proximal and distal counts in Polar area
[Handle,~,F]=setupFigure(1,1,1,'AxesWidth',10,'AxesHeight',6,'XSpace',[3 0.75 0.5]);
set(F,'Visible', 'off'); 
shadeHandles=[];
predHandles=[];
validDistances=polarXTicks(polarXTicks>=minDist);
for cIdx=nonEmptyProcess
        lcountPolarCell=countPolarCell{cIdx};
        lowDistances=(polarXTicks<=median(validDistances))&(polarXTicks>=minDist);
        highDistances=(polarXTicks>=median(validDistances))&(polarXTicks>=minDist);
        
        if(~isempty((lcountPolarCell)))
            notBoxPlot([sum(lcountPolarCell(:,lowDistances),2) sum(lcountPolarCell(:,highDistances),2)]);
        end
end
ylabel(Handle,'Polar Comet counts');
Handle.XTickLabel={'Proximal','Distal'};
printPNGEPSFIG(F,outputDirPlot,'cometCountBoxPolar');

% Polar distance vs cone normalized comet count
[Handle,~,F]=setupFigure(1,1,1,'AxesWidth',10,'AxesHeight',6);
hold on;
shadeHandles=[];
predHandles=[];
for cIdx=nonEmptyProcess  
%     coneNormCountAstralCell{cIdx}=coneNormCountAstralCell{cIdx}./repmat(coneNormCountAstralCell{cIdx}(:,1),1,size(coneNormCountAstralCell{cIdx},2));
%     coneNormCountPolarCell{cIdx}=coneNormCountPolarCell{cIdx}./repmat(coneNormCountPolarCell{cIdx}(:,1),1,size(coneNormCountPolarCell{cIdx},2));
%     
    A=sum(mean(coneNormCountAstralCell{cIdx}))./sum(exp(-astralXTicks*AstralK));
    predictedAstralNormCount= A*exp(-astralXTicks*AstralK);
    
    if(~isempty((coneNormCountAstralCell{cIdx})))
        shadeHandle=plotOrShade(Handle(1),astralXTicks,coneNormCountAstralCell{cIdx},p.shade,cmapAstral(cIdx,:));
        shadeHandles=[shadeHandles shadeHandle];
        pH=plot(Handle(1),astralXTicks,predictedAstralNormCount);
        predHandles=[predHandles pH];
    end    

    if(~isempty((coneNormCountPolarCell{cIdx})))
        shadeHandle=plotOrShade(Handle(1),polarXTicks,coneNormCountPolarCell{cIdx},p.shade,cmapPolar(cIdx,:));
        shadeHandles=[shadeHandles shadeHandle];
    end
end
lineToLegend=arrayfun(@(h) h.mainLine,shadeHandles,'unif',0)
legend([lineToLegend{:} predHandles],{'Astral','Polar','predAstral'},'Location','northeast')
% ylim([0,4])
xlim([0 7]);
xlabel('Pole distance (\mu m)');
ylabel(Handle,'Normalized comet counts per vol (1/\mu m^3)');
hold off;
printPNGEPSFIG(F,outputDirPlot,'normCometCount');
%close(F);

%% Pole distance vs comet density
[H,~,F]=setupFigure(1,1,1,'AxesWidth',10,'AxesHeight',6);
hold on;
shadeHandles=[];
predHandles=[];
for cIdx=nonEmptyProcess
         A=sum(mean(meansAstralDensitiesCell{cIdx}))./sum(exp(-astralXTicks*AstralK)./relativeVolumeAstral);
        predictedAstralAmp= A*exp(-astralXTicks*AstralK)./relativeVolumeAstral;
        
        if(~isempty((meansAstralDensitiesCell{cIdx})))
            shadeHandles=plotOrShade(H(1),astralXTicks,meansAstralDensitiesCell{cIdx},p.shade,cmapAstral(cIdx,:));
            shadeHandles=[shadeHandles shadeHandle];
            pH=plot(H(1),astralXTicks,predictedAstralAmp);
            predHandles=[predHandles pH];
        end
        
        if(~isempty((meansPolarDensitiesCell{cIdx})))
            shadeHandle=plotOrShade(H(1),polarXTicks,meansPolarDensitiesCell{cIdx},p.shade,cmapPolar(cIdx,:));
            shadeHandles=[shadeHandles shadeHandle];
        end
end
lineToLegend=arrayfun(@(h) h.mainLine,shadeHandles,'unif',0)
legend([lineToLegend{:} predHandles],{'Astral','Polar','predAstral'},'Location','northeast')
% ylim([0,4])
xlim([0 7]);
xlabel('Pole distance (\mu m)');
ylabel(H,'Comet density (1/\mu m^3)');
hold off;
printPNGEPSFIG(F,outputDirPlot,'SphereDensityVsDistance');
%close(F)

%  Pole distance vs normalized intensity 
[HAmp,~,FAmp]=setupFigure(1,1,1,'AxesWidth',10,'AxesHeight',6);
hold on;
shadeHandles=[];
predHandles=[];
astralK=0;
for cIdx=nonEmptyProcess
        if(~isempty((normAmpAstralCell{cIdx})))
            shadeHandles=plotOrShade(HAmp(1),astralXTicks,normAmpAstralCell{cIdx},p.shade,cmapAstral(cIdx,:));
            shadeHandles=[shadeHandles shadeHandle];
            % pH=plot(HAmp(1),astralXTicks,predictedAstralAmp);
            % predHandles=[predHandles pH];
        end
        
        if(~isempty((normAmpPolarCell{cIdx})))
            shadeHandle=plotOrShade(HAmp(1),polarXTicks,normAmpPolarCell{cIdx},p.shade,cmapPolar(cIdx,:));
            shadeHandles=[shadeHandles shadeHandle];
        end
end
lineToLegend=arrayfun(@(h) h.mainLine,shadeHandles,'unif',0)
legend([lineToLegend{:} predHandles],{'Astral','Polar',['predAstral K:' num2str(astralK)]},'Location','northeast')
% ylim([0,4])
xlim([0 7]);
xlabel('Pole distance (\mu m)');
ylabel(HAmp,'Normalized intensity');
hold off;
printPNGEPSFIG(F,outputDirPlot,'normAmpVsDistance');

%% Count function of elevation
radianBin=-pi/2:0.1:pi/2;
for cIdx=nonEmptyProcess
    for pIdx=1:length(conditionEnrichmenData(cIdx).movies)
        [counts,edges,binIdx]=histcounts(conditionEnrichmenData(cIdx).movies(pIdx).elevation,radianBin);
        conditionEnrichmenData(cIdx).movies(pIdx).binnedCountvsElevation=counts;
        conditionEnrichmenData(cIdx).movies(pIdx).binningElevation=radianBin;
    end
end
save(fullfile(outputDirPlot,'conditionEnrichmenData.mat'),'conditionEnrichmenData');

%% Projection display

%% For each conditon, display elevation vs density
for cIdx=nonEmptyProcess
    rsize=[300 400];
    cellPlate=cell(ceil(length(enrichmentProcesses{cIdx})/2),2);
	for pIdx=1:length(enrichmentProcesses{cIdx})
        img=imread(sprintfPath(projProcessesAstral{cIdx}(pIdx).outFilePaths_{2},1));
        plotImg=imread(enrichmentProcesses{cIdx}(pIdx).outFilePaths_{1});
        img=imresize(img,rsize);
        plotImg=imresize(plotImg,rsize);
        if(p.printCellName)
            movieName=enrichmentProcesses{cIdx}(pIdx).getOwner().outputDirectory_;
            movieName=strrep(movieName,'/analysis','');
            movieName=strrep(movieName,'/deskew','');
            movieName=strsplit(movieName,'/');
            movieName=[movieName{end-1} ' ' movieName{end}];
            %plotImg=AddTextToImage(plotImg,movieName,size(plotImg));
            plotImg=AddTextToImage(plotImg,movieName,[1 50],[0 0 0], 'Arial',16);
        end
        cellPlate{ceil(pIdx/2),2-mod(pIdx,2)}=[plotImg img];
    end
    emptyMovie=cellfun(@(c) isempty(c),cellPlate);
    if(any(emptyMovie(:)))
        cellPlate{emptyMovie}=uint8(zeros(rsize(1),2*rsize(2),3));
    end
    imwrite(horzcat(vertcat(cellPlate{:,1}),vertcat(cellPlate{:,2})),[outputDirPlot 'plate_' names{cIdx} '.png']);
end

% %% For each conditon, display pole dist vs density
% for cIdx=nonEmptyProcess
%     rsize=[200 400];
%     cellPlate=cell(ceil(length(enrichmentProcesses{cIdx})/2),2);
% 	for pIdx=1:length(enrichmentProcesses{cIdx})
%         img=imread(sprintfPath(projProcessesAstral{cIdx}(pIdx).outFilePaths_{2},1));
%         plotImg=imread(enrichmentProcesses{cIdx}(pIdx).outFilePaths_{4});
%         img=imresize(img,rsize);
%         plotImg=imresize(plotImg,rsize);
%         cellPlate{ceil(pIdx/2),2-mod(pIdx,2)}=[plotImg img];
%         % imwrite([plotImg img],[outputDirPlot 'plate_density_' names{cIdx} '-' num2str(pIdx) '.png']);
%     end
%     emptyMovie=cellfun(@(c) isempty(c),cellPlate);
%     if(any(emptyMovie(:)))    
%         cellPlate{emptyMovie}=uint8(zeros(rsize(1),2*rsize(2),3));
%     end
%     imwrite(horzcat(vertcat(cellPlate{:,1}),vertcat(cellPlate{:,2})),[outputDirPlot 'plate_density_' names{cIdx} '.png']);
% end


printProcessArray({{[projProcessesAstral{:}]} {[projProcessesPolar{:}]}},{[enrichmentProcesses{:}]},4,outputDirPlot,{[names{:}]});
printProcessArray({{[projProcessesAstral{:}]} {[projProcessesPolar{:}]}},{[enrichmentProcesses{:}]},3,outputDirPlot,{[names{:}]});
printProcessArray({projProcessesAstral projProcessesPolar},enrichmentProcesses,5,outputDirPlot,names);
printProcessArray({projProcessesAstral projProcessesPolar},enrichmentProcesses,6,outputDirPlot,names);
printProcessArray({projProcessesAstral projProcessesPolar},enrichmentProcesses,7,outputDirPlot,names);


function render=printProcessArray(processProjCellArray,enrichmentProcesses,outputIndex,outputDirPlot,names)
%% For each conditon, display intensity vs pole dist 
for cIdx=1:length(enrichmentProcesses)
    rsize=[200 400];
    cellPlate=cell(ceil(length(enrichmentProcesses{cIdx})/2),2);
	for pIdx=1:length(enrichmentProcesses{cIdx})
        plotImg=imread(enrichmentProcesses{cIdx}(pIdx).outFilePaths_{outputIndex});
        plotImg=imresize(plotImg,rsize);
        img=plotImg;
        for pcIdx=1:length(processProjCellArray)
            pr=processProjCellArray{pcIdx}{cIdx}(pIdx);
            tmpimg=imread(sprintfPath(pr.outFilePaths_{2},1));
            tmpimg=imresize(tmpimg,rsize);
            img=[img tmpimg];
        end
        cellPlate{ceil(pIdx/2),2-mod(pIdx,2)}=img;
        [~,graphName,~]=fileparts(enrichmentProcesses{cIdx}(pIdx).outFilePaths_{outputIndex});
        % imwrite(img,[outputDirPlot 'plate_' graphName '_' names{cIdx} '-' num2str(pIdx) '.png']);
    end
    emptyMovie=cellfun(@(c) isempty(c),cellPlate);
    if(any(emptyMovie(:)))    
        cellPlate{emptyMovie}=uint8(zeros(rsize(1),3*rsize(2),3));
    end
    [~,graphName,~]=fileparts(enrichmentProcesses{cIdx}(1).outFilePaths_{outputIndex});
    render=horzcat(vertcat(cellPlate{:,1}),vertcat(cellPlate{:,2}));
    imwrite(render,[outputDirPlot 'plate_' graphName '_' names{cIdx} '.png']);
end

function processRenderer=plotCountVsDistance(conditionEnrichmenData,astralXTicks,polarXTicks,names,outputDirPlot,varargin)
    ip = inputParser;
    ip.CaseSensitive = false;
    ip.KeepUnmatched = true;
    ip.addParameter('shade',false);
    ip.addParameter('type','per frame');
    ip.addParameter('lightCond',[]);
    ip.parse(varargin{:});
    p=ip.Results;

    cmap=winter(2*length(conditionEnrichmenData));
    cmapAstral=winter(2*length(conditionEnrichmenData)); % X2 to avoid light colors
    cmapPolar=autumn(2*length(conditionEnrichmenData));

    [Handle,~,F]=setupFigure(1,1,1,'AxesWidth',5,'AxesHeight',5,'AspectRatio',1);
    set(F,'Visible', 'off'); 
    hold on;
    shadeHandles=[];
    AstralShadeHandles=[];
    PolarShadeHandles=[];
    predHandles=[];
    astralK=0;
    minDist=2;  
    fitXTicks=astralXTicks(astralXTicks>=minDist);
    for cIdx=1:length(conditionEnrichmenData)
        switch p.type
        case 'per frame'
            lcountAstralCell=vertcat(conditionEnrichmenData(cIdx).movies.binnedCountAstralPerFrame);
            lcountPolarCell=vertcat(conditionEnrichmenData(cIdx).movies.binnedCountPolarPerFrame); 
        case 'perFramePerKT'
            % For each frame, KT is normalize according to the number of
            lcountAstralCell=vertcat(conditionEnrichmenData(cIdx).movies.binnedCountAstralPerFramePerKT);
            lcountPolarCell=vertcat(conditionEnrichmenData(cIdx).movies.binnedCountPolarPerFramePerKT);
        case 'perFramePerKTAtDist'
            % For each frame, KT is normalize according to the number of
            lcountAstralCell=vertcat(conditionEnrichmenData(cIdx).movies.binnedCountAstralPerFramePerKTAtDist);
            lcountPolarCell=vertcat(conditionEnrichmenData(cIdx).movies.binnedCountPolarPerFramePerKTAtDist);    
        case 'normalized 2um'
            % For each frame, KT is normalize according to the number of
            lcountAstralCell=vertcat(conditionEnrichmenData(cIdx).movies.binnedCountAstralNorm2um);
            lcountPolarCell=vertcat(conditionEnrichmenData(cIdx).movies.binnedCountPolarNorm2um);  
        case 'weighted'
            % For each frame, KT is normalize according to the number of
            lcountAstralCell=vertcat(conditionEnrichmenData(cIdx).movies.binnedCountAstralNorm2umWeighted);
            lcountPolarCell=vertcat(conditionEnrichmenData(cIdx).movies.binnedCountPolarNorm2umWeighted);      
        case 'KT ROI Count'
            % For each frame, KT is normalize according to the number of
            lcountAstralCell=[];
            lcountPolarCell=vertcat(conditionEnrichmenData(cIdx).movies.KTROICounts);       
        end


        % A=sum(mean(lcountPolarCell))./sum(exp(-polarXTicks*PolarK));
        % predictedPolar= A*exp(-polarXTicks*PolarK);
        if(isempty(p.lightCond)||(ismember(cIdx,p.lightCond)))
            transparency=1;
            thickness=1;
        else
            transparency=0.7;
            thickness=0.5;
        end

        if(~isempty((lcountAstralCell)))
            % [f,B]=catastropheRateFitCount(mean(lcountAstralCell(:,astralXTicks>=minDist)),fitXTicks);
            % predictedAstral= f(B,astralXTicks);
            % astralK=B(2);
            % % store fitting results
            % conditionEnrichmenData(cIdx).predictedAstralCount=predictedAstral;
            % conditionEnrichmenData(cIdx).catastropheModel=f;
            % conditionEnrichmenData(cIdx).catastropheModelParam=B;
            % pH=plot(Handle(1),astralXTicks,predictedAstral,'Color',cmapAstral(cIdx,:));
            % predHandles=[predHandles pH];

            shadeHandle=plotOrShade(Handle(1),astralXTicks,lcountAstralCell,p.shade,cmapAstral(cIdx,:),transparency,thickness);
            shadeHandles=[shadeHandles shadeHandle];

        end    

        if(~isempty((lcountPolarCell)))
            shadeHandle=plotOrShade(Handle(1),polarXTicks,lcountPolarCell,p.shade,cmapPolar(cIdx,:),transparency,thickness);
            shadeHandles=[shadeHandles shadeHandle];
        end
    end
    hold off;
    shadeHandles=arrayfun(@(h) h.mainLine,shadeHandles,'unif',0);
    polarLineToLegend=arrayfun(@(h) h.mainLine,PolarShadeHandles,'unif',0);

    c2m=(@(x) [x{:}]);
    if(isempty(p.lightCond))
    % legend([shadeHandles{:}], ... 
    %         c2m(cellfun(@(n) {['Astral-' n],['Polar-' n]},names,'unif',0)), ...
    %         'Location','northwest');
    end
    %legend([lineToLegend{:} predHandles],{'Astral','Polar',['predAstral K:' num2str(astralK)]},'Location','northwest')
    % ylim([0,4])
    xlim(Handle,[0 7]);
    xlabel('Pole distance (\mu m)');
    ylabel(Handle,['+Tips count ' p.type]);
    printPNGEPSFIG(F,outputDirPlot,['cometCount' p.type '_' num2str(p.lightCond)]);
    %% create projection process saving independant projection location
    processRenderer=ExternalProcess(MovieList(),'Renderer');
    processRenderer.setOutFilePaths({[outputDirPlot filesep  'cometCount' p.type '_' num2str(p.lightCond) '.png']});


function shadeHandles=plotOrShade(h,XTicks,measure,shade,cmap,transparency,thickness)
c={'r','b','g','y','k','o'};
shadeHandles=[];
if(size(measure,1)>1)
    y=measure;
    if(shade)
        axes(h);
        shadeHandles=shadedErrorBar(XTicks,mean(y),std(y),{'Color',cmap(1,:)},1);
    else
        p=plot(h,XTicks,mean(y),[c{1} '-'],'Color',cmap(1,:));
        p.Color(4)=transparency;
        p.LineWidth=thickness;
    end
else
    p=plot(h,XTicks,measure,[c{1} '-'],'Color',cmap(1,:));
    p.Color(4)=transparency;
    p.LineWidth=thickness;
end

function [f,B]=catastropheRateFitCount( dataPoint,fitXTicks)

    f = @(b,x) b(1).*exp(-b(2).*x);
    nrmrsd = @(b) norm(dataPoint - f(b,fitXTicks));     % Residual Norm Cost Function
    B0 = [max(dataPoint) 0.5 ];                                            % Choose Appropriate Initial Estimates
    [B,rnrm] = fminsearch(nrmrsd, B0); 
    disp(['Estimated Catastrophe Rate: ' num2str(B(2))])

function [f,B]=catastropheRateFit( dataPoint,fitXTicks)

    f = @(b,x) b(1).*exp(-b(2).*x) + b(3);
    nrmrsd = @(b) norm(dataPoint - f(b,fitXTicks));     % Residual Norm Cost Function
    B0 = [max(dataPoint) 0.5 min(dataPoint)];                                            % Choose Appropriate Initial Estimates
    [B,rnrm] = fminsearch(nrmrsd, B0); 
    disp(['Estimated Catastrophe Rate: ' num2str(B(2))])