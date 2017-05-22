% Testing our new maapping and randomization tested visually on full cells
% that have an expected behavior


%%
%%%%%%% First toy Cell
MD=MovieData.loadMatFile('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/movieData.mat');

%% load spindle ref data and uniform Randomization 
processDetectPoles=MD.getPackage(1).getProcess(1);
processAddSpindleRef=MD.getPackage(1).getProcess(2);
processRandKinAndISO=MD.getPackage(1).getProcess(3);

poleMovieInfo=load(processDetectPoles.outFilePaths_{1}); poleMovieInfo=poleMovieInfo.poleMovieInfo;
[poleRefsISO,P1,P2]=buildSpindleRef(poleMovieInfo,1);

EB3Tracks=load(processAddSpindleRef.outFilePaths_{1}); EB3Tracks=EB3Tracks.EB3Tracks;
kinTracks=load(processAddSpindleRef.outFilePaths_{2}); kinTracks=kinTracks.kinTracks;
EB3TracksInliers=load(processAddSpindleRef.outFilePaths_{3}); EB3TracksInliers=EB3TracksInliers.EB3TracksInliers;
kinTracksInliers=load(processAddSpindleRef.outFilePaths_{4}); kinTracksInliers=kinTracksInliers.kinTracksInliers;

%% load process indenpendant data
outputDirProj=[MD.outputDirectory_ filesep 'Kin' filesep 'track' filesep  ];
kinTrackData=load([outputDirProj  filesep 'tracksLabRef.mat']);
kinTracksISO=kinTrackData.tracksLabRef;
kinTracksISOInliers=kinTracksISO(logical(arrayfun(@(t) t.inliers(1),kinTracks)));

outputDirProj=[MD.outputDirectory_ filesep 'EB3' filesep 'track' filesep  ];
tmp=load([outputDirProj filesep 'tracksLabRef.mat']);
EB3TracksISO=tmp.tracksLabRef;
EB3TracksISOInliers=EB3TracksISO(logical(arrayfun(@(t) t.inliers(1),EB3Tracks)));

%%
%%%%%% Uniform randomization (select good kin bad kin)
%% Compute uniform random data (OPT)
maxRandomDist=20;
processUniformRandom=ExternalProcess(MD,'randomizeTracks');
randomUniform=randomizeTracks(MD,maxRandomDist,'randomType','uniform','tracks',kinTracksISO,'process',processUniformRandom);
MD.getPackage(1).setProcess(3,processUniformRandom);

%% load randomization data
tmp=load(processUniformRandom.outFilePaths_{1});
randKinTracksISO=tmp.randKinTracks;
randKinTracksISOInliers=randKinTracksISO(logical(arrayfun(@(t) t.inliers(1),kinTracks)));

%% Mapping Kin based (+ angle)
tic;
mapMTToKin(kinTracksInliers(1:100),EB3TracksInliers,1000,'distType','normalDistAndAngle','kinDistCutoff',[-2000,100],'position','end');
mapMTToKin(randKinTracksInliers(1:100),EB3TracksInliers,1000,'distType','normalDistAndAngle','kinDistCutoff',[-2000,100],'position','end');
toc;

bundleStatistics(MD,'kinBundle',{kinTracksInliers(1:100),randKinTracksInliers(1:100)}, ... 
                'kinBundleName',{'Inlier','RandomInlier'}, ... 
                'plotName','mappedEnd','mappedMTField','associatedMTP1');
            
%% Mapping manifold based (tube). 
buildFiberManifoldAndMapMT(P1,kinTracksISOInliers(1:100),EB3TracksISOInliers,7,'kinDistCutoff',[-20,20]);
buildFiberManifoldAndMapMT(P1,randKinTracksISOInliers(1:100),EB3TracksISOInliers,7,'kinDistCutoff',[-20,20]);
bundleStatistics(MD,'kinBundle',{kinTracksISOInliers(1:100),randKinTracksISOInliers(1:100)}, ... 
                'kinBundleName',{'Inlier','RandomInlier'}, ... 
                'plotName','mappedEnd','mappedMTField','associatedMT');


%% Select the "N Best" and "N worst" kin
mappedNBDiff=zeros(1,length(kinTracksInliers(1:100)));
for k=1:length(mappedNBDiff)
    mappedNBDiff(k)=length(kinTracksISOInliers(k).associatedMT)-length(randKinTracksISOInliers(k).associatedMT);
end
[diffVal,sortedKinIdx]=sort(mappedNBDiff);

N=10;
capturingKinIdx=sortedKinIdx((end-N):end);
nonCapturingKinIdx=sortedKinIdx(1:N);

bundleStatistics(MD,'kinBundle',{kinTracksISOInliers(capturingKinIdx),randKinTracksISOInliers(capturingKinIdx)}, ...
    'kinBundleName',{'Inlier','RandomInlier'},'plotName','goodKin','mappedMTField','associatedMT');
bundleStatistics(MD,'kinBundle',{kinTracksISOInliers(nonCapturingKinIdx),randKinTracksISOInliers(nonCapturingKinIdx)}, ...
    'kinBundleName',{'Inlier','RandomInlier'},'plotName','nonCapturingUnifRandom','mappedMTField','associatedMT');

%% Render results in Amira 
tic;
outputDirBundle=[MD.outputDirectory_ filesep 'Kin' filesep 'bundles'];
outputDirAmira=[outputDirBundle filesep 'kinMappingCount' filesep];
amiraWriteTracks([outputDirAmira filesep 'kin_.am'],kinTracksInliers(1:100), ...
    'cumulativeOnly',true,'edgeProp',{{'mappedNumber',mappedNBDiff}})
toc

%% project kin and random ROI on random
processProjCapturingKin=projectKinAndRandom(MD,processAddSpindleRef, ...
    processDetectPoles,processUniformRandom,capturingKinIdx,'name','capturingKin','showRand',true);
processProjNonCapturingKin=projectKinAndRandom(MD,processAddSpindleRef, ...
    processDetectPoles,processUniformRandom,nonCapturingKinIdx,'name','nonCapturingKin','showRand',true);

%% overlaying
tic;
overlayProjTracksList(MD,processProjCapturingKin,kinTracksISOInliers(capturingKinIdx),kinTracksISOInliers(capturingKinIdx),randKinTracksISOInliers(capturingKinIdx),randKinTracksISOInliers(capturingKinIdx),P1,'mappedTrackField','associatedMT','name','mapped-randunif')
overlayProjTracksList(MD,processProjNonCapturingKin,kinTracksISOInliers(nonCapturingKinIdx),kinTracksISOInliers(nonCapturingKinIdx),randKinTracksISOInliers(nonCapturingKinIdx),randKinTracksISOInliers(nonCapturingKinIdx),P1,'mappedTrackField','associatedMT','name','mapped-randunif')
toc;


%% View selected kin
kinTest=[2];
tic;
procProjectSelectKin=projectKinAndRandom(MD,processAddSpindleRef, ...
    processDetectPoles,processUniformRandom,kinTest,'name','capturingKin','showRand',false);
procProjectSelectKinRand=projectKinAndRandom(MD,processAddSpindleRef, ...
    processDetectPoles,processUniformRandom,kinTest,'name','capturingKin','showRand',true);
overlayProjTracksList(MD,procProjectSelectKin,kinTracksISOInliers(kinTest),kinTracksISOInliers(kinTest),randKinTracksISOInliers(kinTest),randKinTracksISOInliers(kinTest),P1,'mappedTrackField','associatedMT','name','mapped-randunif')
overlayProjTracksList(MD,procProjectSelectKinRand,kinTracksISOInliers(kinTest),kinTracksISOInliers(kinTest),randKinTracksISOInliers(kinTest),randKinTracksISOInliers(kinTest),P1,'mappedTrackField','associatedMT','name','mapped-randunif-non-rand')
toc;

%% Test new randomization on targeted kin and then full cell
%% Manifold set building 
[manifoldsAllKinP1,subManifoldsAllKinP1]=buildFiberManifold(P1,kinTracksISOInliers,EB3TracksISOInliers,10,'kinDistCutoff',[-20,20]);
[manifoldsAllKinP2,subManifoldsAllKinP2]=buildFiberManifold(P2,kinTracksISOInliers,EB3TracksISOInliers,10,'kinDistCutoff',[-20,20]);

%% captured 
randCapturedKinTracks=randomizeTracks(MD,maxRandomDist,'tracks',kinTracksISOInliers(capturingKinIdx),'mappingDist',mappingDist,'dynManifoldsCell',[subManifoldsAllKinP1]);
buildFiberManifoldAndMapMT(P1,randCapturedKinTracks,EB3TracksISOInliers,5,'kinDistCutoff',[-20,20]);
bundleStatistics(MD,'kinBundle',{kinTracksISOInliers(capturingKinIdx),randCapturedKinTracks},'kinBundleName',{'Inlier','RandomInlier'},'plotName','mappedEnd','mappedMTField','associatedMT');


%% Non-captured
mappingDist=10;
randNonCapturedKinTracks=randomizeTracks(MD,maxRandomDist,'tracks',kinTracksISOInliers(nonCapturingKinIdx),'mappingDist',mappingDist,'dynManifoldsCell',[subManifoldsAllKinP1]);
buildFiberManifoldAndMapMT(P1,randNonCapturedKinTracks,EB3TracksISOInliers,5,'kinDistCutoff',[-20,20]);
bundleStatistics(MD,'kinBundle',{kinTracksISOInliers(nonCapturingKinIdx),randNonCapturedKinTracks},'kinBundleName',{'Inlier','RandomInlier'},'plotName','randomAntispace','mappedMTField','associatedMT');

processProjPyramid=projectKinAndRandom(MD,kinTracksISOInliers(nonCapturingKinIdx),processDetectPoles,randNonCapturedKinTracks,...
                'name','nonCapturing-randAntiSpace','showRand',true);

overlayProjTracksList(MD,processProjPyramid,kinTracksISOInliers(nonCapturingKinIdx),kinTracksISOInliers(nonCapturingKinIdx), ... 
                        randNonCapturedKinTracks,randNonCapturedKinTracks,P1,'mappedTrackField','associatedMT')

%% Randomized non-space all kin
procSupervisedRandom=ExternalProcess(MD,'randomizeTracks');
mapppingDist=10;
randAntiSpaceKinTracks=randomizeTracks(MD,maxRandomDist,'tracks',kinTracksISOInliers(1:100),'mappingDist',mappingDist,'dynManifoldsCell',[subManifoldsAllKinP1],'process',procSupervisedRandom);
buildFiberManifoldAndMapMT(P1,randAntiSpaceKinTracks(1:100),EB3TracksISOInliers,5,'kinDistCutoff',[-20,20]);
bundleStatistics(MD,'kinBundle',{kinTracksISOInliers(1:100),randAntiSpaceKinTracks(1:100)},'kinBundleName',{'Inlier','RandomInlier'},'plotName','randomAntispaceHundredFirst','mappedMTField','associatedMT');


%% MC
maxRandomDist=20;
processUniformRandomMC=ExternalProcess(MD,'randomizeTracks');
[randTracksCell]=randomizeTracksMC(MD,maxRandomDist,'randomType','uniform','tracks',kinTracksISO,'process',processUniformRandomMC,'simuNumber',1000);
%%
load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/Kin/randomized/MC-1000-unif-randKinTracks-mapped.mat')

%%
tic;
for sIdx=1:length(randTracksCell(1:2))
    buildFiberManifoldAndMapMT(P1,randTracksCell{sIdx}(1:100),EB3TracksISOInliers,5,'kinDistCutoff',[-20,20]);
 end
toc;
% tic;
%     procFolder=[processUniformRandomMC.getOwner().outputDirectory_  filesep 'Kin' filesep 'randomized' filesep];
%     mkdirRobust(procFolder);
%     save([procFolder 'MC-' '-' 'unif' '-2-randKinTracks-mapped.mat'],'randTracksCell', '-v7.3');
% toc;

%%
tic;
randTracksCellTruncate=cell(1,length(randTracksCell));
for sIdx=1:length(randTracksCell(1:2))
    randTracksCellTruncate{sIdx}=randTracksCell{sIdx}(1:100);
end

bundleStatistics(MD,'kinBundle',randTracksCellTruncate,'plotName','unifMC','mappedMTField','associatedMT');
toc;