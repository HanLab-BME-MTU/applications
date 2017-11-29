% Testing our new, manifold-based mapping.

% load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/movieData.mat')
MD=MovieData.loadMatFile('C:\Users\Philippe\project-local\externBetzig\analysis\adavid\smallSample\prometaphase\earlyCell1_12\movieData.mat')

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

%% load randomization data
tmp=load(processRandKinAndISO.outFilePaths_{1});
randKinTracks=tmp.randKinTracks;
randKinTracksISO=tmp.randKinTracksISO;
randKinTracksInliers=randKinTracks(logical(arrayfun(@(t) t.inliers(1),kinTracks)));
randKinTracksISOInliers=randKinTracksISO(logical(arrayfun(@(t) t.inliers(1),kinTracks)));

%% Projection for visualization
% Build ref
kinIdx=32;
track=kinTracksISOInliers(kinIdx);

refs=buildRefsFromTracks([P1,P2],track);
processSingleProj=ExternalProcess(MD,'rawProj');
project1D(  MD,[P1,track],'dynPoligonREF',[refs(1).applyBase(P1,[]) refs(1).applyBase(track,[])],'FoF',refs(1), ...
    'name',['bundle-P1-kin-' num2str(kinIdx)],'processSingleProj',processSingleProj);

%%
tic;
mapMTToKin(kinTracksInliers(kinIdx),EB3TracksInliers,1000,'distType','normalDistAndAngle','kinDistCutoff',[-2000,100],'position','end');
toc;

mappedMT=kinTracksInliers(kinIdx).associatedMTP1;
mappedRef=refs(1).applyBase(EB3TracksISO([mappedMT.index]),[]);

%%
track=kinTracksISOInliers(kinIdx);
tic;
manifold=[P1 track];
manifVector= manifold(2).getAddCoord(manifold(1).getMultCoord(-1));
manifNorm=(manifVector.x.^2 + manifVector.y.^2 + manifVector.z.^2).^0.5;
subManifold=[track.getAddCoord(manifVector.getMultCoord(-20./manifNorm)) track.getAddCoord(manifVector.getMultCoord(10./manifNorm))];
mappedTracksP1=mapTracksTo1DManifold(subManifold,EB3TracksISO,1000/MD.pixelSize_,'position','end');
mappedTracksLessP1=mapTracksTo1DManifold(subManifold,EB3TracksISO,500/MD.pixelSize_,'position','end');
toc;

mappedTracksP1Ref=refs(1).applyBase(EB3TracksISO([mappedTracksP1.index]),[]);
mappedTracksLessP1Ref=refs(1).applyBase(EB3TracksISO([mappedTracksLessP1.index]),[]);
overlayProjTracksMovie(processSingleProj,'tracks',[mappedRef; mappedTracksP1Ref;  ],'colorIndx',[  3*ones(1,length(mappedRef)) 1*ones(1,length(mappedTracksP1Ref)) ],'colormap',myColormap,'name','Tube10TubeTubeAngle');

%% test wrapping function
tic;
buildFiberManifoldAndMapMT(P1,track,EB3TracksISO,10)
toc;
%%
mappedTracksP1=track.associatedMT;
mappedTracksP1Ref=refs(1).applyBase(EB3TracksISO([mappedTracksP1.index]),[]);
overlayProjTracksMovie(processSingleProj,'tracks',[ mappedRef; mappedTracksP1Ref; mappedTracksLessP1Ref  ],'colorIndx',[ 3*ones(1,length(mappedRef)) 1*ones(1,length(mappedTracksP1Ref)) 2*ones(1,length(mappedTracksLessP1Ref)) ],'colormap',myColormap,'name','newmap');

%%
tic;
mapMTToKin(kinTracksInliers(1:100),EB3TracksInliers,1000,'distType','normalDistAndAngle','kinDistCutoff',[-2000,100],'position','end');
toc;

%%
mapMTToKin(randKinTracksInliers(1:100),EB3TracksInliers,1000,'distType','normalDistAndAngle','kinDistCutoff',[-2000,100],'position','end');
bundleStatistics(MD,'kinBundle',{kinTracksInliers(1:100),randKinTracksInliers(1:100)},'kinBundleName',{'Inlier','RandomInlier'},'plotName','mappedEnd','mappedMTField','associatedMTP1');

%%
tic;
buildFiberManifoldAndMapMT(P1,kinTracksISOInliers(1:100),EB3TracksISOInliers,10)
toc;

%%
buildFiberManifoldAndMapMT(P1,randKinTracksISOInliers(1:100),EB3TracksISOInliers,10)
bundleStatistics(MD,'kinBundle',{kinTracksISOInliers(1:100),randKinTracksISOInliers(1:100)},'kinBundleName',{'Inlier','RandomInlier'},'plotName','mappedEnd','mappedMTField','associatedMT');

%%
tic;
mappedTracksP1Orig=mapTracksTo1DManifold(subManifold,EB3TracksISO,1000/MD.pixelSize_,'position','end');
    toc;
tic;    
mappedTracksP1Test=mapTracksTo1DManifold(subManifold,EB3TracksISO,1000/MD.pixelSize_,'position','end','distType','normalDistPseudOptimized');
toc;

mappedTracksP1OrigRef=refs(1).applyBase(EB3TracksISO([mappedTracksP1Orig.index]),[]);
mappedTracksP1TestRef=refs(1).applyBase(EB3TracksISO([mappedTracksP1Test.index]),[]);
overlayProjTracksMovie(processSingleProj,'tracks',[mappedTracksP1Ref; mappedTracksP1OrigRef; mappedTracksP1TestRef],'colorIndx',[  2*ones(1,length(mappedTracksP1Ref))  1*ones(1,length(mappedTracksP1OrigRef)) 3*ones(1,length(mappedTracksP1TestRef)) ],'colormap',myColormap,'name','testNewApproach');


%% Test different distance values
kinTest=[2 5 10];

procProjectSelectKin=projectKinAndRandom(MD,processAddSpindleRef, ...
    processDetectPoles,processUniformRandom,kinTest,'name','kinTest-distMeasured','showRand',false);
%%

for dist=[8]
buildFiberManifoldAndMapMT(P1,kinTracksISOInliers(kinTest),EB3TracksISOInliers,dist,'kinDistCutoff',[-20,20]);
bundleStatistics(MD,'kinBundle',{kinTracksISOInliers(kinTest)}, ... 
                'kinBundleName',{'Inlier','RandomInlier'}, ... 
                'plotName','mappedEnd','mappedMTField','associatedMT');
            
overlayProjTracksList(MD,procProjectSelectKin,kinTracksISOInliers(kinTest),kinTracksISOInliers(kinTest),kinTracksISOInliers(kinTest),kinTracksISOInliers(kinTest),P1,'mappedTrackField','associatedMT','name',['dist-' num2str(dist)]);
end 
% overlayProjTracksList(MD,[procProjectSelectKin(1) procProjectSelectKin(1)],kinTracksISOInliers(kinTest),kinTracksISOInliers(kinTest),randKinTracksISOInliers(kinTest),randKinTracksISOInliers(kinTest),P1,'mappedTrackField','associatedMT','name','10')
