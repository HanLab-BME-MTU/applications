function conditionAstralTracksStatistics(processCell,names,outputDirPlot,varargin)
% load scoring data associated to cells. A Cell of process list describes the conditions
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('shade',true);
ip.parse(varargin{:});
p=ip.Results;

if(~iscell(processCell))
    processCell={processCell};
end

mkdirRobust(outputDirPlot);

% load score set in cells of cells
measureCell=cell(1,length(processCell));
accelerationCell=cell(1,length(processCell));

for cIdx=1:length(processCell)
	condScoreCell=cell(1,length(processCell{cIdx}));
	for pIdx=1:length(processCell{cIdx})
		tmp=load(processCell{cIdx}(pIdx).outFilePaths_{1});
		condScoreCell{pIdx}=tmp.distance;
	end
	measureCell{cIdx}=condScoreCell;
	condScoreCell=cell(1,length(processCell{cIdx}));
	for pIdx=1:length(processCell{cIdx})
		tmp=load(processCell{cIdx}(pIdx).outFilePaths_{1});
        try
            condScoreCell{pIdx}=tmp.speedStd;
        catch
        end;
	end
	accelerationCell{cIdx}=condScoreCell;    
end

c={'r','b','g','y','k'};


% Plot mean astral length
allMeasure=[];
group=[];
figure
hold on;
for cIdx=1:length(processCell)
	for pIdx=1:length(processCell{cIdx})
        measures=measureCell{cIdx}{pIdx};
		allMeasure=[allMeasure mean(measures)];
        group=[group cIdx];
    end	
end
boxplot(allMeasure,group,'labels',names)
ylabel('Avg length (\mu m)');
hold off;

print([outputDirPlot  'astralLength.png'],'-dpng');
print([outputDirPlot  'astralLength.eps'],'-depsc');

% Plot astral length histograms
allMeasure=[];
scoresBin=0:0.1:7;
group=[];
astralHist=cell(1,length(processCell));
for cIdx=1:length(processCell)
    conditionAstralCataData.name=names{cIdx};
    moviesStruct=[];
	for pIdx=1:length(processCell{cIdx})
        measures=measureCell{cIdx}{pIdx};
		allMeasure=[allMeasure measures'];
		[counts,edges,binIdx]=histcounts(allMeasure,scoresBin);
		astralHist{cIdx}=[astralHist{cIdx}; counts];
        movie.movieData=processCell{cIdx}(pIdx).getOwner().movieDataPath_;
        movie.trackLengthBin=scoresBin; 
        movie.trackLengthCount=counts; 
        moviesStruct=[moviesStruct movie];
    end
    conditionAstralCataData.movies=moviesStruct;  
end
save(fullfile(outputDirPlot,'conditionAstralCataData.mat'),'conditionAstralCataData');

[Handles,~,F]=setupFigure(length(processCell),1,length(processCell));
for cIdx=1:length(astralHist)
    if(~isempty((astralHist{cIdx})))
        if(size(astralHist{cIdx},1)>1)
            y=astralHist{cIdx};
            if(p.shade)
            	axes(Handles(cIdx))
                H=shadedErrorBar(scoresBin(1:end-1),mean(y),std(y));
      
            else
                plot(scoresBin(1:end-1),mean(y),[c{1} '-']);
            end
        else
            plot(scoresBin(1:end),astralHist{cIdx},[c{1} '-']);
        end
    end
end
arrayfun(@(h) xlabel(h,'Avg length (\mu m)'),Handles);
arrayfun(@(h) ylabel(h,'MT count'),Handles);

hold off;

print([outputDirPlot  'astralLengthHist.png'],'-dpng');
print([outputDirPlot  'astralLengthHist.eps'],'-depsc');


%% Maria's format
allCell=[accelerationCell{:}];
dataMatPadded = reformatDataCell(allCell);
condGroups=arrayfun(@(cIdx) cIdx*ones(1,length(accelerationCell{cIdx})),1:length(accelerationCell),'unif',0);
condGroups = horzcat(condGroups{:}); 
[FCondMean,FCondPooled,FPerCell]=plotDifferentPerspectives(dataMatPadded,condGroups);
printPNGEPSFIG(FCondMean,outputDirPlot,'AccelerationCondMean');
printPNGEPSFIG(FCondPooled,outputDirPlot,'AccelerationCondPooled');
printPNGEPSFIG(FPerCell,outputDirPlot,'AccelerationPerCell');



