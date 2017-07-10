function conditionEnrichmentStatistics(enrichmentProcess,names,outputDirPlot)
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
hold on;
for cIdx=1:length(enrichmentProcess)
	for pIdx=1:length(enrichmentProcess{cIdx})
		allElevs=[elevationCell{cIdx}{pIdx}{:}];
	    allDens=[densityCell{cIdx}{pIdx}{:}];

		[counts,edges,binIdx]=histcounts(allElevs,scoresBin);
		[means,meds,stds,orderedIndex,counts] = statPerIndx(allDens,binIdx+1);

		meansDensitiesCell{cIdx}=[meansDensitiesCell{cIdx}; means];
	end
	
	if(size(histCell{cIdx},1)>1)
        y=meansDensitiesCell{cIdx};
        shadedErrorBar(scoresBin(1:end-1),mean(y),std(y),c{cIdx},1);       
    else
	    plot(scoresBin(1:end-1),meansDensitiesCell{cIdx},[c{cIdx} '-']);
    end
end
legend(names)
hold off;
print([outputDirPlot  'enrichment.png'],'-dpng');
print([outputDirPlot  'enrichment.eps'],'-depsc');




