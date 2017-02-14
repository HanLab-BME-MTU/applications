function [handles,hFig]= bundleStatistics(MD,varargin)
%Plot and compare building 
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addParameter('kinBundle',[]);
ip.addParameter('kinBundleName',[]);
ip.addParameter('bundleMTRange',[10 35]);
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

[handles,~,hFig]=setupFigure(1,2,2,'AspectRatio',1,'AxesWidth',4);

outputDirPlot=[outputDirBundle filesep 'plot' filesep];
system(['mkdir ' outputDirPlot]);

[handles,hFig]= displayBundleStatistics('kinBundle',kinTracksCell,'kinBundleName',p.kinBundleName,'bundleMTRange',p.bundleMTRange);

print([outputDirPlot 'avgMTPerKin-kinCount.png'],'-dpng');
print([outputDirPlot 'avgMTPerKin-kinCount.eps'],'-depsc');



