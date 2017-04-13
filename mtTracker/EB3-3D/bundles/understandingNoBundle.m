load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/movieData.mat')



%% load spindle ref data 
processSpindleRef=MD.getProcess(1);
processRandKin=MD.getProcess(2);

kinTracks=load(processSpindleRef.outFilePaths_{2}); kinTracks=kinTracks.kinTracks;
EB3TracksInliers=load(processSpindleRef.outFilePaths_{3}); EB3TracksInliers=EB3TracksInliers.EB3TracksInliers;
kinTracksInliers=load(processSpindleRef.outFilePaths_{4}); kinTracksInliers=kinTracksInliers.kinTracksInliers;
save('~/kinTracks.mat','kinTracks');
save('~/kinTracksInliers.mat','kinTracksInliers');


%%
% create random kinetochore and measure diff
randomDist=20;
randKinTracks=randomizeKinSpindleRef(MD,randomDist,'processDetectPoles',MD.getPackage(1).getProcess(1));
randKinTracksInliers=randKinTracks(logical(arrayfun(@(t) t.inliers(1),kinTracks)));
MTBundling(MD,kinTracksInliers,EB3TracksInliers,randKinTracksInliers,'kinRange',1:100);
save('~/randKinTracksInliers.mat','randKinTracksInliers');
save('~/randKinTracks.mat','randKinTracks');


%% load random No bundle
randKinTracksNoBundleData=load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/Kin/randomized/randKinTracksNoBundle.mat')
randKinTracksNB=randKinTracksNoBundleData.randKinTracks;
save('~/randKinTracksNoBundle.mat','randKinTracksNB');
randKinTracksNBInliers=randKinTracksNB(logical(arrayfun(@(t) t.inliers(1),kinTracks)))
save('~/randKinTracksNoBundleInliers.mat','randKinTracksNBInliers');
%%

MTBundling(MD,kinTracksInliers,EB3TracksInliers,randKinTracksNBInliers,'kinRange',1:100);

%% Amira shows that XYZ is not the issue but that the NoBundle azimuth are not the same. 
outputDirAmira=[MD.outputDirectory_ filesep 'Kin' filesep 'randomizedSurprise'];

%%
EB3TermAzi=arrayfun(@(t) t.azimuth(1,end),EB3TracksInliers);
EB3TermRho=arrayfun(@(t) t.rho(1,end),EB3TracksInliers);
EB3TermElev=arrayfun(@(t) t.elevation(1,end),EB3TracksInliers);
amiraWriteTracks([outputDirAmira filesep 'EB3.am'],EB3TracksInliers,'cumulativeOnly',true,'edgeProp',{{'Azi',EB3TermAzi},{'rho',EB3TermRho},{'elev',EB3TermElev}})

%%
KinTermAzi=arrayfun(@(t) t.azimuth(1,end),kinTracksInliers);
KinTermRho=arrayfun(@(t) t.rho(1,end),kinTracksInliers);
kinTermElev=arrayfun(@(t) t.elevation(1,end),kinTracksInliers);

rKinTermAzi=arrayfun(@(t) t.azimuth(1,end),randKinTracksInliers);
rKinTermRho=arrayfun(@(t) t.rho(1,end),randKinTracksInliers);
rkinTermElev=arrayfun(@(t) t.elevation(1,end),randKinTracksInliers);

rNBKinTermAzi=arrayfun(@(t) t.azimuth(1,end),randKinTracksNBInliers);
rNBKinTermRho=arrayfun(@(t) t.rho(1,end),randKinTracksNBInliers);
rNBkinTermElev=arrayfun(@(t) t.elevation(1,end),randKinTracksNBInliers);

amiraWriteTracks([outputDirAmira filesep 'kin.am'],kinTracksInliers,'cumulativeOnly',true,'edgeProp',{{'Azi',KinTermAzi},{'rho',KinTermRho},{'elev',kinTermElev}})
amiraWriteTracks([outputDirAmira filesep 'random.am'],randKinTracksInliers,'cumulativeOnly',true,'edgeProp',{{'Azi',rKinTermAzi},{'rho',rKinTermRho},{'elev',rkinTermElev}})
amiraWriteTracks([outputDirAmira filesep 'randomNB.am'],randKinTracksNBInliers,'cumulativeOnly',true,'edgeProp',{{'Azi',rNBKinTermAzi},{'rho',rNBKinTermRho},{'elev',rNBkinTermElev}})

%%



kinTracksInliers(1).azimuth(:,1:5)
randKinTracksInliers(1).azimuth(:,1:5)
randKinTracksNBInliers(1).azimuth(:,1:5)

kinTracksInliers(1).elevation(:,1:5)
randKinTracksInliers(1).elevation(:,1:5)
randKinTracksNBInliers(1).elevation(:,1:5)

%%
kinTracksInliers(1).rho(:,1:5)
randKinTracksInliers(1).rho(:,1:5)
randKinTracksNBInliers(1).rho(:,1:5)




%% Compute or load randomization (sanity check)
randKinTracks=randomizeKinSpindleRef(MD,randomDist,'processDetectPoles',MD.getPackage(1).getProcess(1));

amiraWriteTracks([outputDirAmira filesep 'kin_' num2str(kIdx) '.am'],trackSet,'cumulativeOnly',true,'edgeProp',{{'kinEB3',trackType},{'bundleSize',bundleInfo},{'bundlingMTKin',bundlingMTKinInfo}})


%%
randKinTracks=load(processRandKin.outFilePaths_{1});  randKinTracks=randKinTracks.randKinTracks;

%%
randKinTracksInlier=randKinTracks(logical(arrayfun(@(t) t.inliers(1),kinTracks)));

%%
MTBundling(MD,kinTracksInliers,EB3TracksInliers,randKinTracksInlier);

%%
% It turns out that the response to bundling study is very dissimetric, the
% first few kin show strong differences but the not the overall responses.
% Let s test various part of the kinetochore.
kinIndxBoundariesLow=1:100:length(kinTracksInliers);
kinIndxBoundariesHigh=100:100:length(kinTracksInliers);

parfor b=1:length(kinIndxBoundariesHigh) 
    MTBundling(MD,kinTracksInliers,EB3TracksInliers,randKinTracksInlier, 'kinRange',kinIndxBoundariesLow(b):kinIndxBoundariesHigh(b),'name',['kin-' num2str(kinIndxBoundariesLow(b)) '-' num2str(kinIndxBoundariesHigh(b))]);
end