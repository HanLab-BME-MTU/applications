    % Understanding why our basic MT Bundling does not show differences in the
% between random and kin using local visualization.
% Possible: better mapping or randomization.
% 1 - recompute the best/worst kin
% 2 - single kin computation for debug purposes
% 3 - selected kin computation for generalization

%%
%%%%%%% loading data
%load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/movieData.mat')
MD=MovieData.loadMatFile('C:\Users\Philippe\project-local\externBetzig\analysis\adavid\smallSample\prometaphase\earlyCell1_12\movieData.mat')
MDDecon=MovieData.loadMatFile('C:\Users\Philippe\project-local\externBetzig\analysis\adavid\smallSample\prometaphase\earlyCell1_12\deconv\analysis\movieData.mat');
MDZ=MovieData.loadMatFile('C:\Users\Philippe\project-local\externBetzig\analysis\adavid\smallSample\prometaphase\earlyCell1_12\movieDataPixelSize.mat');

%% load spindle ref data
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
%%%%%% Summarizing uniform randomization issues (select good kin bad kin)
%% Compute random data (OPT)
randomDist=20;
processRandKinAndISO=ExternalProcess(MD,'randomizeTracks');
randomizeKinSpindleRef(MD,randomDist,'processDetectPoles',MD.getPackage(1).getProcess(1),'process',processRandKinAndISO);
MD.getPackage(1).setProcess(3,processRandKinAndISO);

%% load randomization data
tmp=load(processRandKinAndISO.outFilePaths_{1});
randKinTracks=tmp.randKinTracks;
randKinTracksISO=tmp.randKinTracksISO;
randKinTracksInliers=randKinTracks(logical(arrayfun(@(t) t.inliers(1),kinTracks)));
randKinTracksISOInliers=randKinTracksISO(logical(arrayfun(@(t) t.inliers(1),kinTracks)));

%% Mapping Kin based (+ angle)
tic;
mapMTToKin(kinTracksInliers(1:100),EB3TracksInliers,1000,'distType','normalDistAndAngle','kinDistCutoff',[-2000,100],'position','end');
mapMTToKin(randKinTracksInliers(1:100),EB3TracksInliers,1000,'distType','normalDistAndAngle','kinDistCutoff',[-2000,100],'position','end');
toc;

bundleStatistics(MD,'kinBundle',{kinTracksInliers(1:100),randKinTracksInliers(1:100)}, ... 
                'kinBundleName',{'Inlier','RandomInlier'}, ... 
                'plotName','mappedEnd','mappedMTField','associatedMTP1');
            
%% Mapping manifold based (just tube). 
buildFiberManifoldAndMapMT(P1,kinTracksISOInliers(1:100),EB3TracksISOInliers,10,'kinDistCutoff',[-20,20]);
buildFiberManifoldAndMapMT(P1,randKinTracksISOInliers(1:100),EB3TracksISOInliers,10,'kinDistCutoff',[-20,20]);


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

%% Render results in Amira and then 
tic;
outputDirBundle=[MD.outputDirectory_ filesep 'Kin' filesep 'bundles'];
outputDirAmira=[outputDirBundle filesep 'kinMappingCount' filesep];
amiraWriteTracks([outputDirAmira filesep 'kin_.am'],kinTracksInliers(1:100), ...
    'cumulativeOnly',true,'edgeProp',{{'mappedNumber',mappedNBDiff}})
toc

%% project kin and random ROI on random
processProjCapturingKin=projectKinAndRandom(MD,processAddSpindleRef, ...
    processDetectPoles,processRandKinAndISO,capturingKinIdx,'name','capturingKin','showRand',true);
processProjNonCapturingKin=projectKinAndRandom(MD,processAddSpindleRef, ...
    processDetectPoles,processRandKinAndISO,nonCapturingKinIdx,'name','nonCapturingKin','showRand',true);

%% overlaying
tic;
overlayProjTracksList(MD,processProjCapturingKin,kinTracksISOInliers(capturingKinIdx),kinTracksISOInliers(capturingKinIdx),randKinTracksISOInliers(capturingKinIdx),randKinTracksISOInliers(capturingKinIdx),P1,'mappedTrackField','associatedMT','name','mapped-randunif')
overlayProjTracksList(MD,processProjNonCapturingKin,kinTracksISOInliers(nonCapturingKinIdx),kinTracksISOInliers(nonCapturingKinIdx),randKinTracksISOInliers(nonCapturingKinIdx),randKinTracksISOInliers(nonCapturingKinIdx),P1,'mappedTrackField','associatedMT','name','mapped-randunif')
toc;

%%%%%% Projection/debug on a single kin
%% example kinIdx and associated references
kinIdx=11;
zTrack=kinTracksISOInliers(kinIdx);
randZTrack=randKinTracksISOInliers(kinIdx);
refs=buildRefsFromTracks(P1,zTrack);
refsRand=buildRefsFromTracks(P1,randZTrack);

tic;
processProjKinAndRandom=ExternalProcess(MD,'rawProj');
project1D(MD,[P1,randZTrack],'dynPoligonREF',[refs(1).applyBase(P1,[]) refs(1).applyBase(zTrack,[]) refs(1).applyBase(randZTrack,[])],'FoF',refs(1), ...
    'name',['bundle-P1-kin-' num2str(kinIdx) '-R'],'channelRender','grayRed','saveSingleProj',true,'processSingleProj',processProjKinAndRandom,'intMinPrctil',[1 50]);%, ...
toc;
%%
tic;
processProjPyramid=projectKinAndRandom(MD,processAddSpindleRef,processDetectPoles,processRandKinAndISO,kinIdx,...
                'name','testPyramid','showRand',true);
toc;
%% projection on deconv (FAILED because cropped)
tic;
processProjKinAndRandom=ExternalProcess(MD,'rawProj');
project1D(MDDecon,[P1,randZTrack],'dynPoligonREF',[refs(1).applyBase(P1,[]) refs(1).applyBase(zTrack,[]) refs(1).applyBase(randZTrack,[])],'FoF',refs(1), ...
    'name',['deconv-P1-kin-' num2str(kinIdx) '-R'],'channelRender','grayRed','saveSingleProj',true,'processSingleProj',processProjKinAndRandom,'intMinPrctil',[1 50]);%, ...
toc;

%%
tic
myColormap=uint8( ...
    [[0 0 255]; ... % blue "all" tracks
    [0 255 0]; ... % green mapped tracks
    [255 0 0]; ... % kinetochore tracks
    ]);
overlayProjTracksMovie(processProjKinAndRandom,'tracks',[refs(1).applyBase(zTrack,[]); refs(1).applyBase(randZTrack,[])] ,'colorIndx',[3 2],'colormap',myColormap,'name','random');
toc;

% Mapping
tic;
mapMTToKin(kinTracksInliers(kinIdx),EB3TracksInliers,1000,'distType','normalDistAndAngle','kinDistCutoff',[-2000,100],'position','end');
mapMTToKin(randKinTracksInliers(kinIdx),EB3TracksInliers,1000,'distType','normalDistAndAngle','kinDistCutoff',[-2000,100],'position','end');
toc;
bundleStatistics(MD,'kinBundle',{kinTracksInliers(kinIdx),randKinTracksInliers(kinIdx)},'kinBundleName',{'Inlier','RandomInlier'},'plotName','mappedEnd','mappedMTField','associatedMTP1');

% Overlay mapping result locally
tic
mappedMT=kinTracksInliers(kinIdx).associatedMTP1;
mappedRef=refs(1).applyBase(EB3TracksISO([mappedMT.index]),[]);
mappedSpindleRef=poleRefsISO(1).applyBase(EB3TracksISO([mappedMT.index]),[]);

mappedRandMT=randKinTracksInliers(kinIdx).associatedMTP1;
mappedRandRef=refs(1).applyBase(EB3TracksISO([mappedRandMT.index]),[]);
mappedRandSpindleRef=poleRefsISO(1).applyBase(EB3TracksISO([mappedRandMT.index]),[]);

colorIndx=[ones(1,length(mappedMT)) 2*ones(1,length(mappedRandMT))];
overlayProjTracksMovie(processProjKinAndRandom,'tracks',[mappedRef;mappedRandRef;refs(1).applyBase(zTrack,[]); refs(1).applyBase(randZTrack,[])] ,'colorIndx',[colorIndx 3 3],'colormap',myColormap,'name','mapped');
toc;

%% Overlay mapping results on the full Cell
%% Full MIP and random targeted
tic;
origTracks=TracksHandle();
origTracks.x=ones(1,P1.numTimePoints());
origTracks.y=ones(1,P1.numTimePoints());
origTracks.z=ones(1,P1.numTimePoints());
origTracks.startFrame=1;
origTracks.endFrame=P1.numTimePoints();
maxTracks=TracksHandle();
maxTracks.x=MD.imSize_(2)*ones(1,P1.numTimePoints());
maxTracks.y=MD.imSize_(1)*ones(1,P1.numTimePoints());
maxTracks.z=ones(1,P1.numTimePoints())*MD.zSize_*MD.pixelSizeZ_/MD.pixelSize_;
maxTracks.startFrame=1;
maxTracks.endFrame=P1.numTimePoints();
toc;

tic;
dynFullPolygon=[poleRefsISO(1).applyBaseToTrack(origTracks,[]) poleRefsISO(1).applyBaseToTrack(maxTracks,[])];
processProjSpindleRef=ExternalProcess(MD,'rawProj');
project1D(  MD,[P1,zTrack],'FoF',poleRefsISO(1),'dynPoligonREF',dynFullPolygon, ...
            'name','fullSpindleRef','channelRender','grayRed','saveSingleProj',true, ...
            'processSingleProj',processProjSpindleRef,'intMinPrctil',[1 10],'intMaxPrctil',[100 100]);
toc;
%%
tic;
overlayProjTracksMovie(processProjSpindleRef,'tracks',[mappedSpindleRef ;mappedRandSpindleRef;poleRefsISO(1).applyBase(zTrack,[]); poleRefsISO(1).applyBase(randZTrack,[])] ,'colorIndx',[colorIndx 3 3],'colormap',myColormap,'name','mapped');
toc;

%% Amira rendering
tic;
outputDirBundle=[MD.outputDirectory_ filesep 'Kin' filesep 'bundles'];
outputDirAmira=[outputDirBundle filesep 'mappingRand' filesep];
amiraWriteTracks([outputDirAmira filesep 'kin_' num2str(kinIdx) '.am'],mappedRandMT,'cumulativeOnly',true)
toc

%% testing submanifold-based mapping 
tic;
processSingleProjPyramid=projectKinAndRandom(MD,processAddSpindleRef,processDetectPoles,processRandKinAndISO,kinIdx,...
                'name','testPyramid','showRand',true);
toc;

%%
tic;
buildFiberManifoldAndMapMT(P1,kinTracksISOInliers(kinIdx),EB3TracksISOInliers,10,'kinDistCutoff',[-20,20]);
buildFiberManifoldAndMapMT(P1,randKinTracksISOInliers(kinIdx),EB3TracksISOInliers,10,'kinDistCutoff',[-20,20]);
toc;

overlayProjTracksList(MD,processSingleProjPyramid,kinTracksISOInliers(kinIdx),kinTracksISOInliers(kinIdx),randKinTracksISOInliers(kinIdx),randKinTracksISOInliers(kinIdx),P1,'mappedTrackField','associatedMT','name','mapped-randunif')
%%
tic;
buildFiberManifoldAndMapMT(P1,kinTracksISOInliers(kinIdx),EB3TracksISOInliers,10,'kinDistCutoff',[-20,20],'manifoldEntry',true);
buildFiberManifoldAndMapMT(P1,randKinTracksISOInliers(kinIdx),EB3TracksISOInliers,10,'kinDistCutoff',[-20,20],'manifoldEntry',true);
toc; 

overlayProjTracksList(MD,processSingleProjPyramid,kinTracksISOInliers(kinIdx),kinTracksISOInliers(kinIdx),randKinTracksISOInliers(kinIdx),randKinTracksISOInliers(kinIdx),P1,'mappedTrackField','associatedMT','name','mappedDir-randunif')

%%
%%%%%%%%%% New randomization based trained on a few example
%%
kinTest=[32 11 19 21 66];
tic;
processProjPyramid=projectKinAndRandom(MD,processAddSpindleRef,processDetectPoles,processRandKinAndISO,kinTest,...
                'name','testPyramid','showRand',true);
toc;

%%
buildFiberManifoldAndMapMT(P1,kinTracksISOInliers(kinTest),EB3TracksISOInliers,10,'kinDistCutoff',[-20,20]);
buildFiberManifoldAndMapMT(P1,randKinTracksISOInliers(kinTest),EB3TracksISOInliers,10,'kinDistCutoff',[-20,20]);
overlayProjTracksList(MD,processProjPyramid,kinTracksISOInliers(kinTest),kinTracksISOInliers(kinTest),randKinTracksISOInliers(kinTest),randKinTracksISOInliers(kinTest),P1,'mappedTrackField','associatedMT');

buildFiberManifoldAndMapMT(P1,kinTracksISOInliers(kinTest),EB3TracksISOInliers,10,'kinDistCutoff',[-20,20],'manifoldEntry',true);
buildFiberManifoldAndMapMT(P1,randKinTracksISOInliers(kinTest),EB3TracksISOInliers,10,'kinDistCutoff',[-20,20],'manifoldEntry',true);
overlayProjTracksList(MD,processProjPyramid,kinTracksISOInliers(kinTest),kinTracksISOInliers(kinTest),randKinTracksISOInliers(kinTest),randKinTracksISOInliers(kinTest),P1,'mappedTrackField','associatedMT','name','mappedDir');


%% new randomization, single tracks
[manifoldsAllKinP1,subManifoldsAllKinP1]=buildFiberManifold(P1,kinTracksISOInliers,EB3TracksISOInliers,10,'kinDistCutoff',[-20,20]);
[manifoldsAllKinP2,subManifoldsAllKinP2]=buildFiberManifold(P2,kinTracksISOInliers,EB3TracksISOInliers,10,'kinDistCutoff',[-20,20]);

%%
kinIdx=66;
kinTrack=kinTracksISOInliers(kinIdx);
maxRandomDist=30;
mappingDist=10;
[randKinTrack,searchTracks]=randomizeDenseTrackBruteForce(kinTrack,[subManifoldsAllKinP1,subManifoldsAllKinP2],mappingDist,maxRandomDist);

% Describe the Dynamical ROI ( A pyramid that descibe the maximum distance for randomization)
% build pyramid
refKP1=buildRefsFromTracks(P1,kinTrack);
baseX=refKP1.getTracksFromBaseVector('X').getMultCoord(20);
baseY=refKP1.getTracksFromBaseVector('Y').getMultCoord(20);
dynROI=[P1  ...
    kinTrack.getAddCoord(baseX) kinTrack.getAddCoord(baseY) ...
    kinTrack.getAddCoord(baseX.getMultCoord(-1))  kinTrack.getAddCoord(baseY.getMultCoord(-1))];
arrayfun(@(t) refKP1.applyBase(t,'P1K'),dynROI,'unif',0);
insetROI=[P1 randKinTrack];
% Project around the pyramid and show inset.
processSinglePyramid=ExternalProcess(MD,'randProj');
project1D(MD,insetROI,'dynPoligonREF',[dynROI.P1K],'FoF',refKP1, ...
    'name',['testPyramid-Rand-P1-kin-' num2str(kinIdx) '-R'],'channelRender','grayRed','processSingleProj',processSinglePyramid,'intMinPrctil',[1 50]);

overlayProjTracksMovie(processSinglePyramid,'tracks',[refKP1(1).applyBase(searchTracks,[]); refKP1(1).applyBase(kinTrack,[])] ,'colorIndx',[1 3],'colormap',myColormap,'name','BFRandsearchTracks');
overlayProjTracksMovie(processSinglePyramid,'tracks',[refKP1(1).applyBase(randKinTrack,[]); refKP1(1).applyBase(kinTrack,[])] ,'colorIndx',[1 3],'colormap',myColormap,'name','BFRand');

buildFiberManifoldAndMapMT(P1,randKinTrack,EB3TracksISOInliers,10,'kinDistCutoff',[-20,20]);
buildFiberManifoldAndMapMT(P1,kinTrack,EB3TracksISOInliers,10,'kinDistCutoff',[-20,20]);
overlayProjTracksList(MD,processSinglePyramid,kinTrack,kinTrack,randKinTrack,randKinTrack,P1,'mappedTrackField','associatedMT')

%% Randomize captured and non-capture and check diff 
%% Captured 
randCapturedKinTracks=randomizeTracks(MD,maxRandomDist,'tracks',kinTracksISOInliers(capturingKinIdx),'mappingDist',mappingDist,'dynManifoldsCell',[subManifoldsAllKinP1]);
buildFiberManifoldAndMapMT(P1,randCapturedKinTracks,EB3TracksISOInliers,10,'kinDistCutoff',[-20,20]);
bundleStatistics(MD,'kinBundle',{kinTracksISOInliers(capturingKinIdx),randCapturedKinTracks},'kinBundleName',{'Inlier','RandomInlier'},'plotName','mappedEnd','mappedMTField','associatedMT');

%% Non-captured
randNonCapturedKinTracks=randomizeTracks(MD,maxRandomDist,'tracks',kinTracksISOInliers(nonCapturingKinIdx),'mappingDist',mappingDist,'dynManifoldsCell',[subManifoldsAllKinP1]);
buildFiberManifoldAndMapMT(P1,randNonCapturedKinTracks,EB3TracksISOInliers,10,'kinDistCutoff',[-20,20]);
bundleStatistics(MD,'kinBundle',{kinTracksISOInliers(nonCapturingKinIdx),randNonCapturedKinTracks},'kinBundleName',{'Inlier','RandomInlier'},'plotName','randomAntispace','mappedMTField','associatedMT');

processProjPyramid=projectKinAndRandom(MD,kinTracksISOInliers(nonCapturingKinIdx),processDetectPoles,randNonCapturedKinTracks,...
                'name','nonCapturing-randAntiSpace','showRand',true);

overlayProjTracksList(MD,processProjPyramid,kinTracksISOInliers(nonCapturingKinIdx),kinTracksISOInliers(nonCapturingKinIdx), ... 
                        randNonCapturedKinTracks,randNonCapturedKinTracks,P1,'mappedTrackField','associatedMT')

%% Randomized non-space all kin
procSupervisedRandom=ExternalProcess(MD,'randomizeTracks');
randAntiSpaceKinTracks=randomizeTracks(MD,maxRandomDist,'tracks',kinTracksISOInliers(1:100),'mappingDist',mappingDist,'dynManifoldsCell',[subManifoldsAllKinP1],'process',procSupervisedRandom);
buildFiberManifoldAndMapMT(P1,randAntiSpaceKinTracks(1:100),EB3TracksISOInliers,10,'kinDistCutoff',[-20,20]);
bundleStatistics(MD,'kinBundle',{kinTracksISOInliers(1:100),randAntiSpaceKinTracks(1:100)},'kinBundleName',{'Inlier','RandomInlier'},'plotName','randomAntispaceHundredFirst','mappedMTField','associatedMT');
