function [handles,hFig]= displayBundleStatistics(varargin)
%Plot and compare building
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('kinBundle',[]);
ip.addParameter('kinBundleName',[]);
ip.addParameter('plotHandleArray',[]);
ip.addParameter('bundleMTRange',[]);
ip.parse(varargin{:});
p=ip.Results;

%%
kinTracksCell=p.kinBundle;

if(isempty(p.plotHandleArray))
    [handles,~,hFig]=setupFigure(3,5,15,'AspectRatio',1,'AxesWidth',5,'XSPace',[2 2.5 1.5]);
else
    handles=p.plotHandleArray;
end

%% Counting events at each frame
%% Here we refer to pseudo-capture as capture.

% Stats on living Kin, over the kin lifetime.
bundledMTCountOverKinLft=zeros(length(kinTracksCell),kinTracksCell{1}.numTimePoints());
bundledPerCaptureRatioOverKinLft=zeros(length(kinTracksCell),kinTracksCell{1}.numTimePoints());
capturedMTCountOverKinLft=zeros(length(kinTracksCell),kinTracksCell{1}.numTimePoints());
capturedKinOverKinLft=zeros(length(kinTracksCell),kinTracksCell{1}.numTimePoints());
bundledKinOverKinLft=zeros(length(kinTracksCell),kinTracksCell{1}.numTimePoints());

% Stats on living Kin, measured cumulatively.
%cumulBundledPerCaptureRatioPerKLivingKin=zeros(length(kinTracksCell),kinTracksCell{1}.numTimePoints());
cumulCapturedMTCountOverKinLft=zeros(length(kinTracksCell),kinTracksCell{1}.numTimePoints());
cumulBundledMTCountPerLivingKin=zeros(length(kinTracksCell),kinTracksCell{1}.numTimePoints());

% Here we assume that after detection of a pseudo-capturing MT or a bundling MT microtubule
% the Kin is considered captured/bunlded
capturedKin=zeros(length(kinTracksCell),kinTracksCell{1}.numTimePoints());
bundledKin=zeros(length(kinTracksCell),kinTracksCell{1}.numTimePoints());

% Stats on Living MT
livingCapturedMTCountPerKin=zeros(length(kinTracksCell),kinTracksCell{1}.numTimePoints());

livingBundledMTCount=zeros(length(kinTracksCell),kinTracksCell{1}.numTimePoints());
%livingBundledMTlivingBundledMTRatioCountPerKin=zeros(length(kinTracksCell),kinTracksCell{1}.numTimePoints());

kinetochoreCount=zeros(length(kinTracksCell),kinTracksCell{1}.numTimePoints());

for i=1:length(kinTracksCell)
    kinTracks=kinTracksCell{i};
    for k=1:length(kinTracks)
        kin=kinTracks(k);
        capturedMTCountOverKinLft(i,kin.f)=capturedMTCountOverKinLft(i,kin.f)+double(length(kin.fiber));
        if(~isempty(kin.fiber))
            capturedKinOverKinLft(i,kin.f)=capturedKinOverKinLft(i,kin.f)+1;
            bundledPerCaptureRatioOverKinLft(i,kin.f)=bundledPerCaptureRatioOverKinLft(i,kin.f)+double(length(find(kin.fiber)))/double(length(kin.fiber));
        end
        if(any(kin.fiber>0))
            bundledKinOverKinLft(i,kin.f)=bundledKinOverKinLft(i,kin.f)+1;
        end
        bundledMTCountOverKinLft(i,kin.f)=bundledMTCountOverKinLft(i,kin.f)+double(length(find(kin.fiber)));
        kinetochoreCount(i,kin.f)=kinetochoreCount(i,kin.f)+1;
        bundleStartFrame=kin.endFrame;
        captureStartFrame=kin.endFrame;
        for mIdx=1:length(kin.fiber)
            mt=kin.catchingMT(mIdx);
            if(kin.fiber(mIdx)>0)
              livingBundledMTCount(i,mt.f)=livingBundledMTCount(i,mt.f)+1;
              cumulBundledMTCountPerLivingKin(i,mt.startFrame:kin.endFrame)=cumulBundledMTCountPerLivingKin(i,mt.startFrame:kin.endFrame)+1;
              if(mt.startFrame< bundleStartFrame)
                bundleStartFrame=mt.startFrame;
              end
            end
            if(mt.startFrame< captureStartFrame)
              captureStartFrame=mt.startFrame;
            end
            livingCapturedMTCountPerKin(i,mt.f)=livingCapturedMTCountPerKin(i,mt.f)+1;
            cumulCapturedMTCountOverKinLft(i,mt.startFrame:kin.endFrame)=cumulCapturedMTCountOverKinLft(i,mt.startFrame:kin.endFrame)+1;
        end
        capturedKin(i,captureStartFrame:kin.endFrame)=capturedKin(i,captureStartFrame:kin.endFrame)+1;
        bundledKin(i,bundleStartFrame:kin.endFrame)=bundledKin(i,bundleStartFrame:kin.endFrame)+1;
     % livingCapturedMTCount(i,kin.catchingMT(mIdx).f)=
    %  livingBundledMTlivingBundledMTRatioCountPerKin(i,kin.catchingMT(mIdx).f)=livingBundledMTlivingBundledMTRatioCountPerKin(i,kin.catchingMT(mIdx).f)+;

    end
end

% Kinetochore counts
H=handles(1);
plot((H),linspace(0,kinTracks.numTimePoints(),kinTracks.numTimePoints()), kinetochoreCount);
xlabel((H),'Frame count');
ylabel((H),'kinetochore Count');
legend(H,p.kinBundleName);

H=handles(2);
plot(H,linspace(0,kinTracksCell{1}.numTimePoints(),kinTracksCell{1}.numTimePoints()), capturedKinOverKinLft);
xlabel(H,'Frame count');
ylabel(H,{'Living kin with captured MT'});

H=handles(3);
plot(H,linspace(0,kinTracksCell{1}.numTimePoints(),kinTracksCell{1}.numTimePoints()), bundledKinOverKinLft);
xlabel(H,'Frame count');
ylabel(H,{'Living Kin with','bundled MT'});

H=handles(4);
plot(H,linspace(0,kinTracksCell{1}.numTimePoints(),kinTracksCell{1}.numTimePoints()), capturedKin);
xlabel(H,'Frame count');
ylabel(H,{'Captured kin count'});

H=handles(5);
plot(H,linspace(0,kinTracksCell{1}.numTimePoints(),kinTracksCell{1}.numTimePoints()), bundledKin);
xlabel(H,'Frame count');
ylabel(H,{'Bundled kin count'});

% pseudo-capture MT
offset=5;
H=handles(offset+1)
plot(H,linspace(0,kinTracksCell{1}.numTimePoints(),kinTracksCell{1}.numTimePoints()), capturedMTCountOverKinLft./kinetochoreCount);
xlabel(H,'Frame count');
ylabel(H,{'Captured MT count', 'over kin lifetime.'});

H=handles(offset+2);
plot(H,linspace(0,kinTracksCell{1}.numTimePoints(),kinTracksCell{1}.numTimePoints()), cumulCapturedMTCountOverKinLft./kinetochoreCount);
xlabel(H,'Frame count');
ylabel(H,{'Cumulated captured MT', 'count over kin lifetime.'});


H=handles(offset+3)
plot(H,linspace(0,kinTracksCell{1}.numTimePoints(),kinTracksCell{1}.numTimePoints()), livingCapturedMTCountPerKin./kinetochoreCount);
xlabel(H,'Frame count');
ylabel(H,{'Living captured MT', 'per kin'});


% Bundled MT counts
offset=10;
H=handles(offset+1);
plot(H,linspace(0,kinTracksCell{1}.numTimePoints(),kinTracksCell{1}.numTimePoints()), bundledMTCountOverKinLft./kinetochoreCount);
xlabel(H,'Frame count');
ylabel(H,{'Bundle MT count', 'over kin lifetime.'});

H=handles(offset+2);
plot(H,linspace(0,kinTracksCell{1}.numTimePoints(),kinTracksCell{1}.numTimePoints()), cumulBundledMTCountPerLivingKin./kinetochoreCount);
xlabel(H,'Frame count');
ylabel(H,{'cumulative Bundle MT', 'count over kin lifetime.'});

H=handles(offset+3);
plot(H,linspace(0,kinTracksCell{1}.numTimePoints(),kinTracksCell{1}.numTimePoints()), livingBundledMTCount./kinetochoreCount);
xlabel(H,'Frame count');
ylabel(H,{'Living bundled MT count', 'per kin'});

% Bundle to MT Ratio (TODO)

% plot(handles(8),linspace(0,kinTracksCell{1}.numTimePoints(),kinTracksCell{1}.numTimePoints()), livingBundledMTlivingBundledMTRatioCountPerKin./kinetochoreCount);
% xlabel(handles(8),'Frame count');
% ylabel(handles(8),{'Living bundled to captured MT ratio', 'per kin'});

% plot(handles(6),linspace(0,kinTracksCell{1}.numTimePoints(),kinTracksCell{1}.numTimePoints()), bundledPerCaptureRatioOverKinLft./kinetochoreCount);
% xlabel(handles(6),'Frame count');
% ylabel(handles(6),{'Bundled to captured ratio', 'over kin lifetime'});



%%


%%  Kinetochore fiber evolution with time
% Function if Kin average temporal location
% fiberTime=[];
% for i=1:length(kinTracksCell)
%     kinTracks=kinTracksCell{i};
%     kinFiberMask=arrayfun(@(k) isempty(k.fiber),kinTracks);
%     kinTrajTime=arrayfun(@(k) (k.f(ceil(end/2))),kinTracks);
%     fiberTime=[fiberTime; kinTrajTime(kinFiberMask)];
% end
% subplot(1,2,1);
% hist(fiberTime,150);
% legend(p.kinBundleName);


%% Timing of each microtubule
% kinTracks=kinTracksCell{1};
% numTimePoint=zeros(1,2*kinTracks.numTimePoints());
% diffTimingCell=cell(1,length(kinTracks));
% endTimingCell=cell(1,length(kinTracks));
%
% for k=1:length(kinTracks)
%     bundledMTs=kin.catchingMT(kin.fiber>0);
%     diffTiming=zeros(1,length(bundledMTs)-2);
%     endTiming=zeros(1,length(bundledMTs)-2);
%     for mtIdx=2:(length(bundledMTs)-1)
%         diffTiming(mtIdx)= bundledMTs(mtIdx+1).t(end) - bundledMTs(mtIdx).t(end);
%         endTiming(mtIdx)= bundledMTs(mtIdx).t(end) - bundledMTs(1).t(end);
%     end
%     numTimePoint(diffTiming+kinTracks.numTimePoints())=numTimePoint(diffTiming+kinTracks.numTimePoints())+1;
%     diffTimingCell{k}=diffTiming;
%     endTimingCell{k}=endTiming;
% end

% h=setupFigure(1,2,2);
% diffTiming=cell2mat(diffTimingCell);
% endTiming=cell2mat(endTimingCell);
% scatter(h(1),endTiming,diffTiming);
% H=h(1);
% xlim(H,[-5 20])
% xlabel(H,'Frame count after first bundled MT');
% ylabel(H,'Frame count until next bundled MT');
%
% H=h(2);
% plot(H,linspace(-kinTracks.numTimePoints()*MD.timeInterval_,kinTracks.numTimePoints()*MD.timeInterval_,2*kinTracks.numTimePoints()), numTimePoint);
% xlim(H,[-5 20])
% xlabel(H,'Relative Frame count');
% ylabel(H,'frequency');
% print([outputDirPlot 'timing.png'],'-dpng');
% print([outputDirPlot 'timing.eps'],'-depsc');
