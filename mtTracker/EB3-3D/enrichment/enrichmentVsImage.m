function enrichmentVsImage(enrichmentProcess,names,outputDirPlot)
% load scoring data associated to cells. A Cell of process list describes the conditions

if(~iscell(enrichmentProcess))
    enrichmentProcess={enrichmentProcess};
end

mkdirRobust(outputDirPlot);

% load score set in cells of cells
elevationCell=cell(1,length(enrichmentProcess));
densityCell=cell(1,length(enrichmentProcess));

for cIdx=1:length(enrichmentProcess)
	condElevationCell=cell(1,length(enrichmentProcess{cIdx}));
	condDensityCell=cell(1,length(enrichmentProcess{cIdx}));
	for pIdx=1:length(enrichmentProcess{cIdx})
		tmp=load(enrichmentProcess{cIdx}(pIdx).outFilePaths_{2});
		condElevationCell{pIdx}=tmp.elevations;
		condDensityCell{pIdx}=tmp.densities;
	end
	elevationCell{cIdx}=condElevationCell;
	densityCell{cIdx}=condDensityCell;
end

c={'r','b','g','y','k'};
scoresBin=-pi/2:0.05:pi/2;


meansDensitiesCell=cell(1,length(enrichmentProcess));
Handles=[];
[Handle,~,F]=setupFigure(1,1,1,'AxesWidth',10,'AxesHeight',6);
hold on;
for cIdx=1:length(enrichmentProcess)
	for pIdx=1:length(enrichmentProcess{cIdx})
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




