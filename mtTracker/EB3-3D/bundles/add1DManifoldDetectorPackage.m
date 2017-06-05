function add1DManifoldDetectorPackage(MD)
%% Ideally will built a Package to be run later.
%% In practice, for now, we run the whole script through to estimate scores.

% Build a Fake externalProcess for tracking process
processTrackEB3=ExternalProcess(MD,'Tracking');
outputDirProj=[MD.outputDirectory_ filesep 'EB3' filesep 'track' filesep  ];
processTrackEB3.setOutFilePaths({[outputDirProj filesep 'tracksLabRef.mat']});

processTrackKin=ExternalProcess(MD,'Tracking');
outputDirProj=[MD.outputDirectory_ filesep 'Kin' filesep 'track' filesep  ];
processTrackKin.setOutFilePaths({[outputDirProj  filesep 'tracksLabRef.mat']})

processDetectPoles=ExternalProcess(MD,'detectPoles',@(p) detectPoles(p.getOwner(),'process',p,'isoOutput',true));
processDetectPoles.run();

processBuildRef=ExternalProcess(MD,'buildRefsAndROI',@(p) buildRefsFromTracks(processDetectPoles,processDetectPoles,'buildROI',true,'process',p));
processBuildRef.run();

kinTracksISO=load(processTrackKin.outFilePaths_{1}); kinTracksISO=kinTracksISO.tracksLabRef;
ROIs=load(processBuildRef.outFilePaths_{2}); ROIs=ROIs.ROI;
kinTracksISOInliers=mapTracksTo1DManifold(ROIs{1,2},kinTracksISO,0,'position','start','distType','vertexDistOtsu');


tmp=load(processTrackEB3.outFilePaths_{1}); EB3TracksISO=tmp.tracksLabRef;
EB3TracksInliers=mapTracksTo1DManifold(ROIs{1,2},EB3TracksISO,0,'position','start','distType','vertexDistOtsu');
%% Compute scores
kinTest=1:100;
maxRandomDist=20;
mappingDist=10;
processManifoldAntispace=ExternalProcess(MD,'randomizeTracks');
zScores=zeros(2,length(kinTest));
tmp=load(processDetectPoles.outFilePaths_{1});
P1=tmp.tracks(1);
P2=tmp.tracks(2);
for testIdx=1:length(kinTest)
    tIdx=kinTest(testIdx);
    [randTracksCell]=randomizeTracksMC(MD,maxRandomDist,'randomType','manifoldAntispace','tracks',kinTracksISOInliers(tIdx),'process',processManifoldAntispace,'simuNumber',200);
    buildFiberManifoldAndMapMT(P1,[kinTracksISOInliers(tIdx) randTracksCell{:}],EB3TracksInliers,mappingDist,'kinDistCutoff',[-20,20],'mappedTracksField','associatedMTP1');
    buildFiberManifoldAndMapMT(P2,[kinTracksISOInliers(tIdx) randTracksCell{:}],EB3TracksInliers,mappingDist,'kinDistCutoff',[-20,20],'mappedTracksField','associatedMTP2');
    [~,hFigMapped,zScores(1,testIdx)]=bundleStatistics(MD,'kinBundle',[{kinTracksISOInliers(tIdx)} {[randTracksCell{:}]}],'plotName',['manifMC-P1-' num2str(testIdx)],'mappedMTField','associatedMTP1');
    close(hFigMapped);
    [~,hFigMapped,zScores(2,testIdx)]=bundleStatistics(MD,'kinBundle',[{kinTracksISOInliers(tIdx)} {[randTracksCell{:}]}],'plotName',['manifMC-P2-' num2str(testIdx)],'mappedMTField','associatedMTP2');
    close(hFigMapped);
end

[sorterScore,indx]=sort(zScores(:));
[poleIdx,kinIdx]=ind2sub(size(zScores),indx);
save([MD.outputDirectory_ filesep 'Kin' filesep 'bundles' filesep 'bundleStats.mat'],'zScores')

%% Display the N best and worst scores
N=10;
tic
myColormap=uint8( ...
    [[0 255 255]; ... % blue "all" tracks
    [0 255 0]; ... % green mapped tracks
    [255 0 100]; ... % kinetochore tracks
    [255 0 200]; ... % kinetochore tracks
    ]);
for scoreIdx=[1:5 (length(sorterScore)-N):length(sorterScore)]
    tIdx=kinIdx(scoreIdx);
    track=kinTracksISOInliers(tIdx);
    refP1P2=buildRefsFromTracks(P1,P2);
    if(poleIdx(scoreIdx)==1)
        ROI=[P1 track];
        processProj=ExternalProcess(MD,'projP1');
        project1D(MD,ROI,'dynPoligonREF',refP1P2.applyBase([P2 ROI],''),'FoF',refP1P2, ...
            'name',['rand-P1-' num2str(tIdx) '-S-' num2str(sorterScore(scoreIdx))],'channelRender','grayRed','processSingleProj',processProj,'intMinPrctil',[20 98],'intMaxPrctil',[99.9 100],'fringeWidth',30,'insetFringeWidth',mappingDist);
        overlayProjTracksMovie(processProj,'tracks',[refP1P2.applyBase(track.associatedMTP1,[])] ,'colorIndx',[ones(size(track.associatedMTP1))],'colormap',myColormap,'name','test');
        
        refKP=buildRefsFromTracks(P1,track);
        project1D(MD,ROI,'dynPoligonREF',refKP.applyBase([P2 ROI],''),'FoF',refKP, ...
            'name',['PK-rand-P1-' num2str(tIdx) '-S-' num2str(sorterScore(scoreIdx))],'channelRender','grayRed','processSingleProj',processProj,'intMinPrctil',[20 98],'intMaxPrctil',[99.9 100],'fringeWidth',30,'insetFringeWidth',mappingDist);
        overlayProjTracksMovie(processProj,'tracks',[refKP.applyBase(track.associatedMTP1,[])] ,'colorIndx',[ones(size(track.associatedMTP1))],'colormap',myColormap,'name','test');
    end
    
    if(poleIdx(scoreIdx)==2)
        ROI=[P2 track];
        processProj=ExternalProcess(MD,'projP2');
        project1D(MD,ROI,'dynPoligonREF',refP1P2.applyBase([P2 ROI],''),'FoF',refP1P2, ...
            'name',['rand-P2-' num2str(tIdx) '-S-' num2str(sorterScore(scoreIdx))],'channelRender','grayRed','processSingleProj',processProj,'intMinPrctil',[20 98],'intMaxPrctil',[99.9 100],'fringeWidth',30,'insetFringeWidth',mappingDist);
        overlayProjTracksMovie(processProj,'tracks',[refP1P2.applyBase(track.associatedMTP2,[])] ,'colorIndx',[ones(size(track.associatedMTP2))],'colormap',myColormap,'name','test');

       refKP=buildRefsFromTracks(P2,track);
        project1D(MD,ROI,'dynPoligonREF',refKP.applyBase([P2 ROI],''),'FoF',refKP, ...
            'name',['PK-rand-P2-' num2str(tIdx) '-S-' num2str(sorterScore(scoreIdx))],'channelRender','grayRed','processSingleProj',processProj,'intMinPrctil',[20 98],'intMaxPrctil',[99.9 100],'fringeWidth',30,'insetFringeWidth',mappingDist);
        overlayProjTracksMovie(processProj,'tracks',[refKP.applyBase(track.associatedMTP1,[])] ,'colorIndx',[ones(size(track.associatedMTP1))],'colormap',myColormap,'name','test');
    
    end
end


