function [ h ] = generateDistanceDensityPlot( MD, laminChannel, histoneChannel )
%generateDistanceDensityPlot 

MD.sanityCheck;

reader = CellReader(MD.getReader());

procID = MD.getProcessIndex('lamins.classes.LaminsAnalysisProcess');
proc = MD.processes_{procID};
path = proc.outFilePaths_{laminChannel};
params = proc.getParameters();
if(~isfield(params,'analysisDate'))
    params = proc.getDefaultParams();
end
data = load([path filesep 'skeletons_' params.analysisDate '.mat']);

I.lamins = reader{laminChannel,data.tz(1)};
I.histones = reader{histoneChannel,data.tz(1)};

outpath = [path filesep 'DistanceDensityPlots'];
mkdir(outpath);

name = strrep(MD.getFilename,'_','\_');

h(1) = figure;
configurePaper(h(1));
title(name);
imshowpair(imadjust(I.lamins),imadjust(I.histones));
data.S3{data.tz(1)}.drawEdgesAsLines([],'c');
savefig(h(1),[outpath filesep 'merge.fig']);
print(h(1),[outpath filesep 'merge.png'],'-dpng');

h(2) = data.S3{data.tz(1)}.distanceDensityPlot(I.histones);
configurePaper(h(2));
title(name);
savefig(h(2),[outpath filesep 'histone_euclidean.fig']);
print(h(2),[outpath filesep 'histone_euclidean.png'],'-dpng');

h(3) = data.S3{data.tz(1)}.distanceDensityPlot(I.histones,'metric','cityblock');
configurePaper(h(3));
title(name);
savefig(h(3),[outpath filesep 'histone_cityblock.fig']);
print(h(3),[outpath filesep 'histone_cityblock.png'],'-dpng');

h(4) = data.S3{data.tz(1)}.distanceDensityPlot(I.lamins);
configurePaper(h(4));
title(name);
savefig(h(4),[outpath filesep 'lamin_euclidean.fig']);
print(h(4),[outpath filesep 'lamin_euclidean.png'],'-dpng');

h(5) = data.S3{data.tz(1)}.distanceDensityPlot(I.lamins,'metric','cityblock');
configurePaper(h(5));
title(name);
savefig(h(5),[outpath filesep 'lamin_cityblock.fig']);
print(h(5),[outpath filesep 'lamin_cityblock.png'],'-dpng');

h(6) = data.S3{data.tz(1)}.distanceDensityPlot(I.histones,'metric','roundEuclidean');
configurePaper(h(2));
title(name);
savefig(h(6),[outpath filesep 'histone_round_euclidean.fig']);
print(h(6),[outpath filesep 'histone_round_euclidean.png'],'-dpng');

h(7) = data.S3{data.tz(1)}.distanceDensityPlot(I.lamins,'metric','roundEuclidean');
configurePaper(h(7));
title(name);
savefig(h(7),[outpath filesep 'lamin_round_euclidean.fig']);
print(h(7),[outpath filesep 'lamin_round_euclidean.png'],'-dpng');

h(8) = figure;
configurePaper(h(8));
title(name);
imshow(imadjust(I.lamins),[]);
data.S3{data.tz(1)}.drawEdgesAsLines([],'r');
savefig(h(8),[outpath filesep 'lamin_merge.fig']);
print(h(8),[outpath filesep 'lamin_merge.png'],'-dpng');



end
function configurePaper(h)
    set(h,'PaperSize',[ 3 3]);
end

