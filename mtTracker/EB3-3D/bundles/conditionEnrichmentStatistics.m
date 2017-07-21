function conditionEnrichmentStatistics(enrichmentProcesses,projProcesses,names,outputDirPlot)
% load scoring data associated to cells. A Cell of process list describes the conditions

if(~iscell(enrichmentProcesses))
    enrichmentProcesses={enrichmentProcesses};
end

mkdirRobust(outputDirPlot);

% load score set in cells of cells
elevationCell=cell(1,length(enrichmentProcesses));
densityCell=cell(1,length(enrichmentProcesses));

for cIdx=1:length(enrichmentProcesses)
	condElevationCell=cell(1,length(enrichmentProcesses{cIdx}));
	condDensityCell=cell(1,length(enrichmentProcesses{cIdx}));
	for pIdx=1:length(enrichmentProcesses{cIdx})
		tmp=load(enrichmentProcesses{cIdx}(pIdx).outFilePaths_{2});
		condElevationCell{pIdx}=tmp.elevations;
		condDensityCell{pIdx}=tmp.densities;
	end
	elevationCell{cIdx}=condElevationCell;
	densityCell{cIdx}=condDensityCell;
end

c={'r','b','g','y','k'};
scoresBin=-pi/2:0.05:pi/2;


meansDensitiesCell=cell(1,length(enrichmentProcesses));
Handles=[];
[Handle,~,F]=setupFigure(1,1,1,'AxesWidth',10,'AxesHeight',6);
hold on;
for cIdx=1:length(enrichmentProcesses)
	for pIdx=1:length(enrichmentProcesses{cIdx})
		allElevs=vertcat(elevationCell{cIdx}{pIdx}{:});
	    allDens=vertcat(densityCell{cIdx}{pIdx}{:});

		[counts,edges,binIdx]=histcounts(allElevs,scoresBin);
		[means,meds,stds,orderedIndex,counts] = statPerIndx(allDens,binIdx+1);

		meansDensitiesCell{cIdx}=[meansDensitiesCell{cIdx}; means];
	end
	
	if(size(meansDensitiesCell{cIdx},1)>1)
        y=meansDensitiesCell{cIdx};
        H=shadedErrorBar(scoresBin(1:end),mean(y),std(y),c{cIdx},1);
        Handles=[Handles H];
    else
	    plot(scoresBin(1:end),meansDensitiesCell{cIdx},[c{cIdx} '-']);
    end
end
lineToLegend=arrayfun(@(h) h.mainLine,Handles,'unif',0)
legend([lineToLegend{:}],names)
ylim([0,4])
xlim([-pi/2+0.05,1.5]);
xlabel('Elevation (rad)');
ylabel('Comet density (1/\mu m^3)')
hold off;
print([outputDirPlot  'enrichment.png'],'-dpng');
print([outputDirPlot  'enrichment.eps'],'-depsc');

%% For each conditon, display a density and cell image
for cIdx=1:length(enrichmentProcesses)
    rsize=[300 400];
    cellPlate=cell(ceil(length(enrichmentProcesses{cIdx})/2),2);
	for pIdx=1:length(enrichmentProcesses{cIdx})
        img=imread(sprintfPath(projProcesses{cIdx}(pIdx).outFilePaths_{1},projProcesses{cIdx}(pIdx).getOwner().nFrames_));
        plotImg=imread(enrichmentProcesses{cIdx}(pIdx).outFilePaths_{1});
        img=imresize(img,rsize);
        plotImg=imresize(plotImg,rsize);
        cellPlate{ceil(pIdx/2),2-mod(pIdx,2)}=[plotImg img];
    end
    emptyMovie=cellfun(@(c) isempty(c),cellPlate);
    if(any(emptyMovie(:)))
        cellPlate{emptyMovie}=uint8(zeros(rsize(1),2*rsize(2),3));
    end
    imwrite(horzcat(vertcat(cellPlate{:,1}),vertcat(cellPlate{:,2})),[outputDirPlot 'plate_' num2str(cIdx) '.png']);
end


%% For each conditon, display a total count and cell image
for cIdx=1:length(enrichmentProcesses)
    rsize=[200 400];
    cellPlate=cell(ceil(length(enrichmentProcesses{cIdx})/2),2);
	for pIdx=1:length(enrichmentProcesses{cIdx})
        img=imread(sprintfPath(projProcesses{cIdx}(pIdx).outFilePaths_{1},projProcesses{cIdx}(pIdx).getOwner().nFrames_));
        plotImg=imread(enrichmentProcesses{cIdx}(pIdx).outFilePaths_{3});
        img=imresize(img,rsize);
        plotImg=imresize(plotImg,rsize);
        cellPlate{ceil(pIdx/2),2-mod(pIdx,2)}=[plotImg img];
    end
    emptyMovie=cellfun(@(c) isempty(c),cellPlate);
    if(any(emptyMovie(:)))    
        cellPlate{emptyMovie}=uint8(zeros(rsize(1),2*rsize(2),3));
    end
    imwrite(horzcat(vertcat(cellPlate{:,1}),vertcat(cellPlate{:,2})),[outputDirPlot 'plate_totalCount_' num2str(cIdx) '.png']);
end

