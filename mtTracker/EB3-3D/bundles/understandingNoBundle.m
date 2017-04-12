load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/movieData.mat')



%% load spindle ref data 
processSpindleRef=MD.getProcess(1);
processRandKin=MD.getProcess(2);

kinTracks=load(processSpindleRef.outFilePaths_{2}); kinTracks=kinTracks.kinTracks;
EB3TracksInliers=load(processSpindleRef.outFilePaths_{3}); EB3TracksInliers=EB3TracksInliers.EB3TracksInliers;
kinTracksInliers=load(processSpindleRef.outFilePaths_{4}); kinTracksInliers=kinTracksInliers.kinTracksInliers;

%%
% create random kinetochore and measure diff
randomDist=20;
randKinTracks=randomizeKinSpindleRef(MD,randomDist,'processDetectPoles',MD.getPackage(1).getProcess(1));
randKinTracksInliers=randKinTracks(logical(arrayfun(@(t) t.inliers(1),kinTracks)));
MTBundling(MD,kinTracksInliers,EB3TracksInliers,randKinTracksInliers,'kinRange',1:100);


%% load random No bundle
randKinTracksNoBundleData=load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/Kin/randomized/randKinTracksNoBundle.mat')
randKinTracksNB=randKinTracksNoBundleData.randKinTracks;
save('~/randKinTracksNoBundle.mat','randKinTracksNB');
randKinTracksNBInliers=randKinTracksNB(logical(arrayfun(@(t) t.inliers(1),kinTracks)))
save('~/randKinTracksNoBundleInliers.mat','randKinTracksNBInliers');
%%

MTBundling(MD,kinTracksInliers,EB3TracksInliers,randKinTracksNBInliers,'kinRange',1:100);

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