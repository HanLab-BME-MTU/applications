function bundleStatistics(MD,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addParameter('printAll',false, @islogical);
ip.addParameter('testKinIdx',[19 46 156],@isnumeric);
ip.parse(MD,varargin{:});

%%
outputDirBundle=[MD.outputDirectory_ filesep 'Kin' filesep 'bundles'];
tmp=load([outputDirBundle filesep 'kin-MT-bundle.mat'],'kinTracks');
kinTracks=tmp.kinTracks;

%testKinIdx=p.testKinIdx;


outputDirPlot=[outputDirBundle filesep 'plot' filesep];
system(['mkdir -p ' outputDirPlot]);

%%  Kinetochore fiber evolution with time
% Function if Kin average temporal location
kinFiberMask=arrayfun(@(k) isempty(k.fiber),kinTracks);
kinTrajTime=arrayfun(@(k) (k.f(ceil(end/2))),kinTracks);
fiberTime=kinTrajTime(kinFiberMask);
figure()
hist(fiberTime,150);


%% Number of fiber at each given frame. 
fiberCount=zeros(1,kinTracks.numTimePoints);
kinetochoreCount=zeros(1,kinTracks.numTimePoints);
for k=1:length(kinTracks)
    fiberCount(kinTracks(k).f)=fiberCount(kinTracks(k).f)+double(length(kinTracks(k).fiber));
    kinetochoreCount(kinTracks(k).f)=kinetochoreCount(kinTracks(k).f)+1;
end
figure();
plot(linspace(0,kinTracks.numTimePoints*MD.timeInterval_,kinTracks.numTimePoints), fiberCount./kinetochoreCount);
xlabel('Frame count');
ylabel('avg MT per bundle');
print([outputDirPlot 'avgMTPerPlot.png'],'-dpng');
print([outputDirPlot 'avgMTPerPlot.eps'],'-depsc');

%%
figure();
plot(linspace(0,kinTracks.numTimePoints*MD.timeInterval_,kinTracks.numTimePoints), kinetochoreCount);
xlabel('Frame count');
ylabel('kinetochore Count');
print([outputDirPlot 'kinCount.png'],'-dpng');
print([outputDirPlot 'kinCount.eps'],'-depsc');

%% Timing of each microtubule
numTimePoint=zeros(1,2*kinTracks.numTimePoints);
diffTimingCell=cell(1,length(kinTracks));
endTimingCell=cell(1,length(kinTracks));

for k=1:length(kinTracks)
    bundledMTs=kinTracks(k).catchingMT(kinTracks(k).fiber>0);
    diffTiming=zeros(1,length(bundledMTs)-2);
    endTiming=zeros(1,length(bundledMTs)-2);
    for mtIdx=2:(length(bundledMTs)-1)
        diffTiming(mtIdx)= bundledMTs(mtIdx+1).t(end) - bundledMTs(mtIdx).t(end);
        endTiming(mtIdx)= bundledMTs(mtIdx).t(end) - bundledMTs(1).t(end);      
    end
    numTimePoint(diffTiming+kinTracks.numTimePoints)=numTimePoint(diffTiming+kinTracks.numTimePoints)+1;
    diffTimingCell{k}=diffTiming;
    endTimingCell{k}=endTiming;
end
%%
h=setupFigure(1,2,2);
diffTiming=cell2mat(diffTimingCell);
endTiming=cell2mat(endTimingCell);
scatter(h(1),endTiming,diffTiming);
H=h(1);
xlim(H,[-5 20])
xlabel(H,'Frame count after first bundled MT');
ylabel(H,'Frame count until next bundled MT');

H=h(2);
plot(H,linspace(-kinTracks.numTimePoints*MD.timeInterval_,kinTracks.numTimePoints*MD.timeInterval_,2*kinTracks.numTimePoints), numTimePoint);
xlim(H,[-5 20])
xlabel(H,'Relative Frame count');
ylabel(H,'frequency');
print([outputDirPlot 'timing.png'],'-dpng');
print([outputDirPlot 'timing.eps'],'-depsc');

