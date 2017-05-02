% Understanding why our basic MT Bundling does not show differences in the 
% between random and kin using local visualization. 
% Possible: better mapping or randomization.
% 1 - recompute the best/worst kin 
% 2 - single kin computation for debug purposes
% 3 - selected kin computation for generalization 

%load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/movieData.mat')
MD=MovieData.loadMatFile('C:\Users\Philippe\project-local\externBetzig\analysis\adavid\smallSample\prometaphase\earlyCell1_12\movieData.mat')
MDDecon=MovieData.loadMatFile('C:\Users\Philippe\project-local\externBetzig\analysis\adavid\smallSample\prometaphase\earlyCell1_12\deconv\analysis\movieData.mat');
MDZ=MovieData.loadMatFile('C:\Users\Philippe\project-local\externBetzig\analysis\adavid\smallSample\prometaphase\earlyCell1_12\movieDataPixelSize.mat');

%% load spindle ref data
processDetectPoles=MD.getPackage(1).getProcess(1);
processAddSpindleRef=MD.getPackage(1).getProcess(2);
processRandKinAndISO=MD.getPackage(1).getProcess(3);

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

%%%%%% Summarizing naive randomization issues (select good kin bad kin on
%%%%%% 100 first kin)
%% Compute random data (OPT)
randomDist=20;
processRandKinAndISO=ExternalProcess(MD,'rawProj');
randomizeKinSpindleRef(MD,randomDist,'processDetectPoles',MD.getPackage(1).getProcess(1),'process',processRandKinAndISO);
MD.getPackage(1).setProcess(3,processRandKinAndISO);

%% load randomization data
tmp=load(processRandKinAndISO.outFilePaths_{1});
randKinTracks=tmp.randKinTracks;
randKinTracksISO=tmp.randKinTracksISO;
randKinTracksInliers=randKinTracks(logical(arrayfun(@(t) t.inliers(1),kinTracks)));
randKinTracksISOInliers=randKinTracksISO(logical(arrayfun(@(t) t.inliers(1),kinTracks)));

poleMovieInfo=load(processDetectPoles.outFilePaths_{1}); poleMovieInfo=poleMovieInfo.poleMovieInfo;
[poleRefsISO,P1,P2]=buildSpindleRef(poleMovieInfo,1);

%% Mapping
tic;
mapMTToKin(kinTracksInliers(1:100),EB3TracksInliers,1000,'distType','normalDistAndAngle','kinDistCutoff',[-2000,100],'position','end');
mapMTToKin(randKinTracksInliers(1:100),EB3TracksInliers,1000,'distType','normalDistAndAngle','kinDistCutoff',[-2000,100],'position','end');
toc;
bundleStatistics(MD,'kinBundle',{kinTracksInliers(1:100),randKinTracksInliers(1:100)},'kinBundleName',{'Inlier','RandomInlier'},'plotName','mappedEnd','mappedMTField','associatedMTP1');

%% Select the "N Best" and "N worst" kin
mappedNBDiff=zeros(1,length(kinTracksInliers(1:100)));
for k=1:length(mappedNBDiff)
    mappedNBDiff(k)=length(kinTracksInliers(k).associatedMTP1)-length(randKinTracksInliers(k).associatedMTP1);
end
[diffVal,sortedKinIdx]=sort(mappedNBDiff);

N=10;
capturingKinIdx=sortedKinIdx((end-N):end)
nonCapturingKinIdx=sortedKinIdx(1:N)

%% check bundling
bundleStatistics(MD,'kinBundle',{kinTracksInliers(capturingKinIdx),randKinTracksInliers(capturingKinIdx)},'kinBundleName',{'Inlier','RandomInlier'},'plotName','goodKin','mappedMTField','associatedMTP1');
bundleStatistics(MD,'kinBundle',{kinTracksInliers(nonCapturingKinIdx),randKinTracksInliers(nonCapturingKinIdx)},'kinBundleName',{'Inlier','RandomInlier'},'plotName','badKin','mappedMTField','associatedMTP1');

%% Render results
tic;
outputDirBundle=[MD.outputDirectory_ filesep 'Kin' filesep 'bundles'];
outputDirAmira=[outputDirBundle filesep 'kinMappingCount' filesep];
amiraWriteTracks([outputDirAmira filesep 'kin_.am'],kinTracksInliers(1:100),'cumulativeOnly',true,'edgeProp',{{'mappedNumber',mappedNBDiff}})
toc

%% project kin and random ROI on random
processProjCapturingKin=projectKinAndRandom(MD,processAddSpindleRef,processDetectPoles,processRandKinAndISO,capturingKinIdx,'name','capturingKin','showRand',true);
processProjNonCapturingKin=projectKinAndRandom(MD,processAddSpindleRef,processDetectPoles,processRandKinAndISO,nonCapturingKinIdx,'name','nonCapturingKin','showRand',true);

%% overlaying 
tic;
overlayProjTracksList(MD,processProjCapturingKin,kinTracksInliers(capturingKinIdx),kinTracksISOInliers(capturingKinIdx),randKinTracksInliers(capturingKinIdx),randKinTracksISOInliers(capturingKinIdx),P1)
overlayProjTracksList(MD,processProjNonCapturingKin,kinTracksInliers(nonCapturingKinIdx),kinTracksISOInliers(nonCapturingKinIdx),randKinTracksInliers(nonCapturingKinIdx),randKinTracksISOInliers(nonCapturingKinIdx),P1)
toc;

%%%%%% Projection/debug on a single kin
%% example kinIdx and associated references 
kinIdx=32;
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
%% projection on deconv
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

%%%%%%%%%% New randomization based trained on a few example
%%
kinTest=[32 44 11 19];
tic;
processProjPyramid=projectKinAndRandom(MD,processAddSpindleRef,processDetectPoles,processRandKinAndISO,kinTest,... 
                'name','testPyramid','showRand',true);
toc;
%%
mapMTToKin(kinTracksInliers(kinTest),EB3TracksInliers,1000,'distType','normalDistAndAngle','kinDistCutoff',[-2000,100],'position','end');
mapMTToKin(randKinTracksInliers(kinTest),EB3TracksInliers,1000,'distType','normalDistAndAngle','kinDistCutoff',[-2000,100],'position','end');
%%
overlayProjTracksList(MD,processProjPyramid,kinTracksInliers(kinTest),kinTracksISOInliers(kinTest),randKinTracksInliers(kinTest),randKinTracksISOInliers(kinTest),P1)




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
