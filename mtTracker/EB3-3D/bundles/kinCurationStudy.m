% Understanding why our basic MT Bundling does not show differences in the 
% between random and kin using local visualization. 
% Possible: better mapping or randomization.

%load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/movieData.mat')
load('C:\Users\Philippe\project-local\externBetzig\analysis\adavid\smallSample\prometaphase\earlyCell1_12\movieData.mat')



%% load spindle ref data
processDetectPoles=MD.getPackage(1).getProcess(1);
processAddSpindleRef=MD.getPackage(1).getProcess(2);
processRandKin=MD.getPackage(1).getProcess(3);

kinTracks=load(processAddSpindleRef.outFilePaths_{2}); kinTracks=kinTracks.kinTracks;
EB3TracksInliers=load(processAddSpindleRef.outFilePaths_{3}); EB3TracksInliers=EB3TracksInliers.EB3TracksInliers;
kinTracksInliers=load(processAddSpindleRef.outFilePaths_{4}); kinTracksInliers=kinTracksInliers.kinTracksInliers;

%% load randomization data
randKinTracks=load(processRandKin.outFilePaths_{1}); randKinTracks=randKinTracks.randKinTracks;
randKinTracksInliers=randKinTracks(logical(arrayfun(@(t) t.inliers(1),kinTracks)));




%% load process indenpendant data 
outputDirProj=[MD.outputDirectory_ filesep 'Kin' filesep 'track' filesep  ];
kinTrackData=load([outputDirProj  filesep 'tracksLabRef.mat']);
kinTracksISO=kinTrackData.tracksLabRef;
kinTracksISOInliers=kinTracksISO(logical(arrayfun(@(t) t.inliers(1),kinTracks)));

outputDirProj=[MD.outputDirectory_ filesep 'EB3' filesep 'track' filesep  ];
tmp=load([outputDirProj filesep 'tracksLabRef.mat']);
EB3TracksISO=tmp.tracksLabRef;

%%
processSpindleRef=ExternalProcess(MD,'buildSpindleRef',@(p) buildSpindleRef('processDetectPoles',processDetectPoles,'process',p));
processSpindleRef.run();

%% Compute or load randomization (sanity check)
randomDist=20;
processRandKinAndISO=ExternalProcess(MD,'rawProj');
[randKinTracks,randKinTracksISO]=randomizeKinSpindleRef(MD,randomDist,'processDetectPoles',MD.getPackage(1).getProcess(1),'process',processRandKinAndISO);
randKinTracksInliers=randKinTracks(logical(arrayfun(@(t) t.inliers(1),kinTracks)));
randKinTracksISOInliers=randKinTracksISO(logical(arrayfun(@(t) t.inliers(1),kinTracks)));

%% Start with a single kinetochore, obserbe random kin and capture  behavior compared with normal kinetochore 
kinIdx=32;
MTBundling(MD,kinTracksInliers,EB3TracksInliers,randKinTracksInliers,'kinRange',kinIdx);

%% Collect capture and bundling results and project in ad-hoc ref.
zTrack=kinTracksISOInliers(kinIdx);
randZTrack=randKinTracksISOInliers(kinIdx);

poleMovieInfo=load(processDetectPoles.outFilePaths_{1}); poleMovieInfo=poleMovieInfo.poleMovieInfo;
[poleRefsISO,P1,P2]=buildSpindleRef(poleMovieInfo,1);
refs=buildRefsFromTracks(P1,zTrack);
refsRand=buildRefsFromTracks(P1,randZTrack);

%EB3TracksISOREF=refs(1).applyBase(EB3TracksISO,'P1K')
catchingMTNMREF=kinTracksInliers(kinIdx).catchingMTKinRef;
capturedIndx=[catchingMTNMREF.index];
capturedTrackISOREF=refs(1).applyBase(EB3TracksISO(capturedIndx),'P1K');

% this is going to be represented in kinPole ref, hence muse be projected in adequate environment.
catchingRandMTNMREF=randKinTracksInliers(kinIdx).catchingMTKinRef;
capturedRandIndx=[catchingRandMTNMREF.index];
capturedRandTrackISOREF=refs(1).applyBase(EB3TracksISO(capturedRandIndx),'P1K');


%% CreateFullMIP in spindle ref
% still observing kin-random but also localizing the random positions. 
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
%%
tic;
dynFullPolygon=[poleRefsISO(1).applyBaseToTrack(origTracks,[]) poleRefsISO(1).applyBaseToTrack(maxTracks,[])];
processProjSpindleRef=ExternalProcess(MD,'rawProj');
project1D(  MD,[P1,zTrack],'FoF',poleRefsISO(1),'dynPoligonREF',dynFullPolygon, ... 
            'name','fullSpindleRef','channelRender','grayRed','saveSingleProj',true,'processSingleProj',processProjSpindleRef);
toc;


%% Produce projection for the kinetochore only
tic;
processSingleProj=ExternalProcess(MD,'rawProj');
project1D(MD,[P1,zTrack],'dynPoligonREF',[refs(1).applyBase(P1,[]) refs(1).applyBase(zTrack,[])],'FoF',refs(1), ... 
    'name',['bundle-P1-kin-' num2str(kinIdx)],'channelRender','grayRed','saveSingleProj',true,'processSingleProj',processSingleProj);%, ... 
    %'tracks',EB3TracksISOREF,'colorIndx',colorIndx,'colormap',myColormap);
toc;
%%
tic
myColormap=uint8( ... 
    [[0 0 255]; ... % blue "all" tracks
    [0 255 0]; ... % green mapped tracks
    [255 0 0]; ... % kinetochore tracks   
    ]);
colorIndx=ones(1,length(capturedTrackISOREF));
overlayProjTracksMovie(processSingleProj,'tracks',[capturedTrackISOREF; refs(1).applyBase(zTrack,[])],'colorIndx',[colorIndx 3],'colormap',myColormap);
toc

%% Observe kin and random at the same time
% still observing kin-random but also localizing the random positions. 
tic;
processProjKinAndRandom=ExternalProcess(MD,'rawProj');
project1D(MD,[P1,zTrack],'dynPoligonREF',[refs(1).applyBase(P1,[]) refs(1).applyBase(zTrack,[]) refs(1).applyBase(randZTrack,[])],'FoF',refs(1), ... 
    'name',['bundle-P1-kin-' num2str(kinIdx) '-R'],'channelRender','grayRed','saveSingleProj',true,'processSingleProj',processProjKinAndRandom);%, ...
toc;
%%
tic
colorIndx=[ones(1,length(capturedTrackISOREF)) 2*ones(1,length(capturedRandTrackISOREF))];
overlayProjTracksMovie(processProjKinAndRandom,'tracks',[capturedTrackISOREF; capturedRandTrackISOREF;refs(1).applyBase(zTrack,[]); refs(1).applyBase(randZTrack,[])] ,'colorIndx',[colorIndx 3 3],'colormap',myColormap,'name','captured');
toc

%%
tic;
mapMTToKin(kinTracksInliers(kinIdx),EB3TracksInliers,1000,'distType','normalDistAndAngle','kinDistCutoff',[-2000,100],'position','end');
mapMTToKin(randKinTracksInliers(kinIdx),EB3TracksInliers,1000,'distType','normalDistAndAngle','kinDistCutoff',[-2000,100],'position','end');
toc;
%%
bundleStatistics(MD,'kinBundle',{kinTracksInliers(kinIdx),randKinTracksInliers(kinIdx)},'kinBundleName',{'Inlier','RandomInlier'},'plotName','mappedEnd','mappedMTField','associatedMTP1');

%%
tic
mappedMT=kinTracksInliers(kinIdx).associatedMTP1;
mappedRef=refs(1).applyBase(EB3TracksISO([mappedMT.index]),[]); 
mappedSpindleRef=poleRefsISO(1).applyBase(EB3TracksISO([mappedMT.index]),[]); 

mappedRandMT=randKinTracksInliers(kinIdx).associatedMTP1;
mappedRandRef=refs(1).applyBase(EB3TracksISO([mappedRandMT.index]),[]); 
mappedRandSpindleRef=poleRefsISO(1).applyBase(EB3TracksISO([mappedRandMT.index]),[]); 

colorIndx=[ones(1,length(mappedMT)) 2*ones(1,length(mappedRandMT))];
overlayProjTracksMovie(processProjKinAndRandom,'tracks',[mappedRef;mappedRandRef;refs(1).applyBase(zTrack,[]); refs(1).applyBase(randZTrack,[])] ,'colorIndx',[colorIndx 3 3],'colormap',myColormap,'name','mapped');
%%
overlayProjTracksMovie(processProjSpindleRef,'tracks',[mappedSpindleRef ;mappedRandSpindleRef;poleRefsISO(1).applyBase(zTrack,[]); poleRefsISO(1).applyBase(randZTrack,[])] ,'colorIndx',[colorIndx 3 3],'colormap',myColormap,'name','mapped');
toc;

%%
tic;
outputDirBundle=[MD.outputDirectory_ filesep 'Kin' filesep 'bundles'];
outputDirAmira=[outputDirBundle filesep 'mappingRand' filesep];
amiraWriteTracks([outputDirAmira filesep 'kin_' num2str(kinIdx) '.am'],mappedRandMT,'cumulativeOnly',true)
toc


%%%% test on 100 first kin 
%% Mapping
tic;
mapMTToKin(kinTracksInliers(1:100),EB3TracksInliers,1000,'distType','normalDistAndAngle','kinDistCutoff',[-2000,100],'position','end');
mapMTToKin(randKinTracksInliers(1:100),EB3TracksInliers,1000,'distType','normalDistAndAngle','kinDistCutoff',[-2000,100],'position','end');
toc;

%% Stats
bundleStatistics(MD,'kinBundle',{kinTracksInliers(1:100),randKinTracksInliers(1:100)},'kinBundleName',{'Inlier','RandomInlier'},'plotName','mappedEnd','mappedMTField','associatedMTP1');

%% Select the "N Best" and "N worst" kin
mappedNBDiff=zeros(1,length(kinTracksInliers(1:100)));
for k=1:length(mappedNBDiff)
    mappedNBDiff(k)=length(kinTracksInliers(k).associatedMTP1)-length(randKinTracksInliers(k).associatedMTP1);
end

[diffVal,sortedKinIdx]=sort(mappedNBDiff);

N=5;
capturingKinIdx=sortedKinIdx((end-N):end)
nonCapturingKinIdx=sortedKinIdx(1:N)

%%
tic;
outputDirBundle=[MD.outputDirectory_ filesep 'Kin' filesep 'bundles'];
outputDirAmira=[outputDirBundle filesep 'kinMappingCount' filesep];
amiraWriteTracks([outputDirAmira filesep 'kin_.am'],kinTracksInliers(1:100),'cumulativeOnly',true,'edgeProp',{{'mappedNumber',mappedNBDiff}})
toc


%%
bundleStatistics(MD,'kinBundle',{kinTracksInliers(capturingKinIdx),randKinTracksInliers(capturingKinIdx)},'kinBundleName',{'Inlier','RandomInlier'},'plotName','goodKin','mappedMTField','associatedMTP1');
bundleStatistics(MD,'kinBundle',{kinTracksInliers(nonCapturingKinIdx),randKinTracksInliers(nonCapturingKinIdx)},'kinBundleName',{'Inlier','RandomInlier'},'plotName','badKin','mappedMTField','associatedMTP1');



%% project kin and random 
processProjCapturingKin=projectKinAndRandom(MD,processAddSpindleRef,processDetectPoles,processRandKinAndISO,capturingKinIdx,'name','capturingKin');
processProjNonCapturingKin=projectKinAndRandom(MD,processAddSpindleRef,processDetectPoles,processRandKinAndISO,nonCapturingKinIdx,'name','nonCapturingKin');

%%
tic;
overlayProjTracksList(MD,processProjCapturingKin,kinTracksInliers(capturingKinIdx),kinTracksISOInliers(capturingKinIdx),randKinTracksInliers(capturingKinIdx),randKinTracksISOInliers(capturingKinIdx),P1)
overlayProjTracksList(MD,processProjNonCapturingKin,kinTracksInliers(nonCapturingKinIdx),kinTracksISOInliers(nonCapturingKinIdx),randKinTracksInliers(nonCapturingKinIdx),randKinTracksISOInliers(nonCapturingKinIdx),P1)
toc;
%% SNIPPET AND CODE TEMP TRASH
% It turns out that the response to bundling study is very dissimetric, the
% first few kin show strong differences but the not the overall responses.
% Let s test various part of the kinetochore.

% kinIndxBoundariesLow=1:100:length(kinTracksInliers);
% kinIndxBoundariesHigh=100:100:length(kinTracksInliers);
% 
% parfor b=1:length(kinIndxBoundariesHigh)
%     MTBundling(MD,kinTracksInliers,EB3TracksInliers,randKinTracksInlier, 'kinRange',kinIndxBoundariesLow(b):kinIndxBoundariesHigh(b),'name',['kin-' num2str(kinIndxBoundariesLow(b)) '-' num2str(kinIndxBoundariesHigh(b))]);
% end
