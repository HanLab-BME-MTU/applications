function package=add1DManifoldDetectorPackage(MD)
%% Ideally will built a Package to be run letter
%% In practice, for now, we run the whole script through to estimate scores.

% Build a Fake externalProcess for tracking process
processTrackEB3=ExternalProcess(MD,'Tracking');
outputDirProj=[MD.outputDirectory_ filesep 'EB3' filesep 'track' filesep  ];
processTrackEB3.setOutFilePaths([outputDirProj filesep 'tracksLabRef.mat']);

processTrackKin=ExternalProcess(MD,'Tracking');
outputDirProj=[MD.outputDirectory_ filesep 'Kin' filesep 'track' filesep  ];
processTrackKin.setOutFilePaths({[outputDirProj  filesep 'tracksLabRef.mat']})

processDetectPoles=ExternalProcess(MD,'detectPoles',@(p) detectPoles(p.getOwner(),'process',p,'isoOutput',true));
processDetectPoles.run();

processBuildRef=ExternalProcess(MD,'buildRefsAndROI',@(p) buildRefsFromTracks(processDetectPoles,processDetectPoles,'buildROI',true,'process',p));

tmp=load(processTrackKin.outFilePaths_
mappedTracks=mapTracksTo1DManifold(ROIs{1,2},tmp.tracksLabRef,0,'position','start','distType','vertexDistOtsu');


kinTracksInliers=

kinTest=1:100;
maxRandomDist=20;
processUniformRandomMC=ExternalProcess(MD,'randomizeTracks');
tic;
[randTracksCell]=randomizeTracksMC(MD,maxRandomDist,'randomType','manifoldAntispace','tracks',kinTracksISOInliers(kinTest),'process',processUniformRandomMC,'simuNumber',200);
toc;




processProj=ExternalProcess(MD,'projection',@(p) projectFrameOfRef(MD,ref,'process',p));





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
overlayProjTracksList(MD,procProjectSelectKin,kinTracksISOInliers(kinTest),kinTracksISOInliers(kinTest),kinTracksISOInliers(kinTest),kinTracksISOInliers(kinTest),P1,'mappedTrackField','associatedMT','name','mapped-randunif')
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

%%
load('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/Kin/randomized/MC-1000-unif-randKinTracks-mapped.mat')

%%

%%
tic
for kinTest=1:150
    
    mappingDist=10;
    track=kinTracksISOInliers(kinTest);
    refKP1=buildRefsFromTracks(P1,track);
    refP1P2=buildRefsFromTracks(P1,P2);

    ROI=[P1 track];
    processProj=ExternalProcess(MD,'proj');
    myColormap=uint8( ...
        [[0 0 255]; ... % blue "all" tracks
        [0 255 0]; ... % green mapped tracks
        [255 0 100]; ... % kinetochore tracks
        [255 0 200]; ... % kinetochore tracks
        ]);
%     project1D(MD,ROI,'dynPoligonREF',refKP1.applyBase([ROI],''),'FoF',refKP1, ...
%         'name',['rand-' num2str(kinTest)],'channelRender','grayRed','processSingleProj',processProj,'intMinPrctil',[1 95],'intMaxPrctil',[100 100],'fringeWidth',30,'insetFringeWidth',mappingDist);
   project1D(MD,ROI,'dynPoligonREF',refP1P2.applyBase([P2 ROI],''),'FoF',refP1P2, ...
        'name',['rand-' num2str(kinTest)],'channelRender','grayRed','processSingleProj',processProj,'intMinPrctil',[20 99],'intMaxPrctil',[99.9 100],'fringeWidth',30,'insetFringeWidth',mappingDist);
%     [randTracksCell]=randomizeTracksMC(MD,maxRandomDist,'randomType','manifoldAntispace','tracks',kinTracksISOInliers(kinTest),'process',processUniformRandomMC,'simuNumber',200);
%     
%     buildFiberManifoldAndMapMT(P1,[kinTracksISOInliers(kinTest)' [randTracksCell{:}]],EB3TracksISOInliers,10,'kinDistCutoff',[-20,20]);
%     overlayProjTracksMovie(processProj,'tracks',[refKP1.applyBase(track.associatedMT,[])] ,'colorIndx',[ones(size(track.associatedMT))],'colormap',myColormap,'name','test');
% 
%     bundleStatistics(MD,'kinBundle',[{kinTracksISOInliers(kinTest)} {[randTracksCell{:}]}],'plotName',['manifMC-kin-' num2str(kinTest)],'mappedMTField','associatedMT');
end

toc;
%%


tic;
% for sIdx=1:length(randTracksCell)
 buildFiberManifoldAndMapMT(P1,[kinTracksISOInliers(kinTest)' [randTracksCell{:}]],EB3TracksISOInliers,10,'kinDistCutoff',[-20,20]);
% end
toc;

% tic;
%     procFolder=[processUniformRandomMC.getOwner().outputDirectory_  filesep 'Kin' filesep 'randomized' filesep];
%     mkdirRobust(procFolder);
%     save([procFolder 'MC-' '-' 'unif' '-2-randKinTracks-mapped.mat'],'randTracksCell', '-v7.3');
% toc;
% 
% 
% tic;
% randTracksCellTruncate=cell(1,length(randTracksCell));
% for sIdx=1:length(randTracksCell)
%     randTracksCellTruncate{sIdx}=randTracksCell{sIdx};
% end
%
tic
bundleStatistics(MD,'kinBundle',[{kinTracksISOInliers(kinTest)} {[randTracksCell{:}]}],'plotName','manifMC-100b','mappedMTField','associatedMT');
toc;