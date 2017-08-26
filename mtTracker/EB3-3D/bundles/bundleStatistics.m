function [handles,hFigMapped,zscore]= bundleStatistics(MD,varargin)
%Plot and compare building 
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addParameter('kinBundle',[]);
ip.addParameter('kinBundleName',[]);
ip.addParameter('bundleMTRange',[]);
ip.addParameter('mappedMTField','capturedMT');
ip.addParameter('bundledMTField','fiber');
ip.addParameter('plotName',[]);
ip.parse(MD,varargin{:});
p=ip.Results;

%%
outputDirBundle=[MD.outputDirectory_ filesep 'Kin' filesep 'bundles'];
if(isempty(p.kinBundle))
    tmp=load([outputDirBundle filesep 'kin-MT-bundle.mat'],'kinTracks');
    kinTracksCell={tmp.kinTracks};
else
    kinTracksCell=p.kinBundle;
end

[handlesMapped,~,hFigMapped]=setupFigure(1,2,2,'AspectRatio',1,'AxesWidth',5,'XSPace',[2 2.5 1.5]);
[handlesNull,~,hFigNull]=setupFigure(3,5,15,'AspectRatio',1,'AxesWidth',5,'XSPace',[2 2.5 1.5]);

handlesNull(7:8)=handlesMapped;
outputDirPlot=[outputDirBundle filesep 'plot' filesep];
mkdirRobust(outputDirPlot);
[handles,hFig,zscore]= displayBundleStatistics('kinBundle',kinTracksCell,'plotHandleArray',handlesNull,varargin{:});
print([outputDirPlot p.plotName '_bundleStat.png'],'-dpng');
print([outputDirPlot p.plotName '_bundleStat.eps'],'-depsc');
close(hFigNull);


