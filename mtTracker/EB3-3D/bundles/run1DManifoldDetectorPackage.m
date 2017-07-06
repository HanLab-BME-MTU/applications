function package=run1DManifoldDetectorPackage(MD,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('package',[]);
ip.parse(varargin{:});
p=ip.Results;

%% Ideally will built a Package to be run later.
%% However, external process does not have a simple way yet to edit the parameter of its associated funtion. 
%% In practice, for now, we run the whole script through to estimate scores.

% Process type placeholdes
MD.setPackage(333,GenericPackage({ ... 
    PointSourceDetectionProcess3D(MD,[MD.outputDirectory_ filesep 'EB3']),...
    TrackingProcess(MD,[MD.outputDirectory_ filesep 'EB3']),...
    PointSourceDetectionProcess3D(MD,[MD.outputDirectory_ filesep 'KT']),...
    TrackingProcess(MD,[MD.outputDirectory_ filesep 'KT']),...  
    ExternalProcess(MD,'detectPoles'),...
    ExternalProcess(MD,'buildRefsAndROI'),...
    ExternalProcess(MD,'manifoldScoring')...
  }));


if(~isempty(p.package)&&(~isempty(p.package.getProcess(1))))
   processDetectEB3=p.package.getProcess(1);
else
    processDetectEB3=PointSourceDetectionProcess3D(MD, [MD.outputDirectory_ filesep 'EB3'],UTrackPackage3D.getDefaultDetectionParams(MD,[MD.outputDirectory_ filesep 'EB3']));
    MD.addProcess(processDetectEB3);
    funParams = processDetectEB3.funParams_;
    funParams.showAll=true;
    funParams.alpha=0.05;
    funParams.filterSigma=[1.4396 1.4396;1.2913 1.2913  ];
    funParams.WindowSize={[],[]};
    funParams.algorithmType= {'pointSourceLM'  'pointSourceLM'}
    funParams.ConfRadius={[],[]};       
    processDetectEB3.setPara(funParams);
    paramsIn.ChannelIndex=1;
    paramsIn.isoCoord=true;
    processDetectEB3.run(paramsIn);
end

MD.getPackage(333).setProcess(1,processDetectEB3);

if(~isempty(p.package)&&(~isempty(p.package.getProcess(2))))
    processTrackEB3=p.package.getProcess(2);
else
    processTrackEB3=TrackingProcess(MD, [MD.outputDirectory_ filesep 'EB3'],UTrackPackage3D.getDefaultTrackingParams(MD, [MD.outputDirectory_ filesep 'EB3']));
    MD.addProcess(processTrackEB3);    
    funParams = processTrackEB3.funParams_;
    [costMatrices,gapCloseParam,kalmanFunctions,probDim]=plusTipCometTracker3DParam(MD);
    funParams.gapCloseParam=gapCloseParam;
    funParams.costMatrices=costMatrices;
    funParams.kalmanFunctions=kalmanFunctions;
    funParams.probDim=probDim;
    processTrackEB3.setPara(funParams);
    paramsIn.ChannelIndex=1;
    paramsIn.DetProcessIndex=processDetectEB3.getIndex();
    processTrackEB3.run(paramsIn);
end

MD.getPackage(333).setProcess(2,processTrackEB3);

if(~isempty(p.package)&&(~isempty(p.package.getProcess(3))))
   processDetectKT=p.package.getProcess(3);
else
    processDetectKT=PointSourceDetectionProcess3D(MD,[MD.outputDirectory_ filesep 'KT'],UTrackPackage3D.getDefaultDetectionParams(MD, [MD.outputDirectory_ filesep 'KT']));
    MD.addProcess(processDetectKT);
    funParams = processDetectKT.funParams_;
    funParams.showAll=true;
    funParams.alpha=0.05;
    funParams.filterSigma=[1.6 1.6;1.5 1.5  ];
    funParams.WindowSize={[],[]};
    funParams.algorithmType= {'pointSourceFit'  'pointSourceFit'}
    funParams.ConfRadius={[],[]};       
    processDetectKT.setPara(funParams);
    paramsIn.ChannelIndex=2;
    paramsIn.isoCoord=true;
    processDetectKT.run(paramsIn);    
end

MD.getPackage(333).setProcess(3,processDetectKT);
                    
if(~isempty(p.package)&&(~isempty(p.package.getProcess(4))))
    processTrackKT=p.package.getProcess(4);
else
    processTrackKT=TrackingProcess(MD, [MD.outputDirectory_ filesep 'KT'],UTrackPackage3D.getDefaultTrackingParams(MD,[MD.outputDirectory_ filesep 'KT']));
    MD.addProcess(processTrackKT);    
    funParams = processTrackKT.funParams_;
    [gapCloseParam,costMatrices,kalmanFunctions,probDim,verbose]=kinTrackingParam();
    funParams.gapCloseParam=gapCloseParam;
    funParams.costMatrices=costMatrices;
    funParams.kalmanFunctions=kalmanFunctions;
    funParams.probDim=probDim;
    processTrackKT.setPara(funParams);
    paramsIn.ChannelIndex=2;
    paramsIn.DetProcessIndex=processDetectKT.getIndex();
    processTrackKT.run(paramsIn);
end

MD.getPackage(333).setProcess(4,processTrackKT);


% Build a Fake externalProcess for tracking process
% processTrackEB3=ExternalProcess(MD,'Tracking');
% outputDirProj=[MD.outputDirectory_ filesep 'EB3' filesep 'track' filesep  ];
% processTrackEB3.setOutFilePaths({[outputDirProj filesep 'tracksLabRef.mat']});
% 
% processTrackKT=ExternalProcess(MD,'Tracking');
% outputDirProj=[MD.outputDirectory_ filesep 'Kin' filesep 'track' filesep  ];
% processTrackKT.setOutFilePaths({[outputDirProj  filesep 'tracksLabRef.mat']})

if(~isempty(p.package)&&(~isempty(p.package.getProcess(5))))
    processDetectPoles=p.package.getProcess(5);
else
    processDetectPoles=ExternalProcess(MD,'detectPoles',@(p) detectPoles(p.getOwner(),'process',p,'isoOutput',true));
    processDetectPoles.run();
end

MD.getPackage(333).setProcess(5,processDetectPoles);


if(~isempty(p.package)&&(~isempty(p.package.getProcess(6))))
    processBuildRef=p.package.getProcess(6);
else
    processBuildRef=ExternalProcess(MD,'buildRefsAndROI',@(p) buildRefsFromTracks(processDetectPoles,processDetectPoles,'buildROI',true,'process',p));
    processBuildRef.run();
end

MD.getPackage(333).setProcess(6,processBuildRef);


%% Load tracks, convert and select inliers
ROIs=load(processBuildRef.outFilePaths_{2}); ROIs=ROIs.ROI;
tmp=load(processTrackKT.outFilePaths_{2}); kinTracksISO=TracksHandle(tmp.tracksFinal);
kinTracksISOInliers=mapTracksTo1DManifold(ROIs{1,2},kinTracksISO,0,'position','start','distType','vertexDistOtsu');

tmp=load(processTrackEB3.outFilePaths_{1}); EB3TracksISO=TracksHandle(tmp.tracksFinal);
EB3TracksInliers=mapTracksTo1DManifold(ROIs{1,2},EB3TracksISO,0,'position','start','distType','vertexDistOtsu');

clear tmp;

tmp=load(processDetectPoles.outFilePaths_{1});
P1=tmp.tracks(1);
P2=tmp.tracks(2);
%% Compute scores
if(~isempty(p.package)&&(~isempty(p.package.getProcess(7))))
    processScoring=p.package.getProcess(7);
else
    lftThreshold=20;
    kinTest=find([kinTracksISOInliers.lifetime]>lftThreshold);
    kinTest=[1 100];
    maxRandomDist=20;
    mappingDist=10;
    processManifoldAntispace=ExternalProcess(MD,'randomizeTracks');
    zScores=nan(2,length(kinTracksISOInliers));

    for testIdx=1:length(kinTest)
        tIdx=kinTest(testIdx);
        [randTracksCell]=randomizeTracksMC(MD,maxRandomDist,'randomType','manifoldAntispace','tracks',kinTracksISOInliers(tIdx),'process',processManifoldAntispace,'simuNumber',200);
        buildFiberManifoldAndMapMT(P1,[kinTracksISOInliers(tIdx) randTracksCell{:}],EB3TracksInliers,mappingDist,'kinDistCutoff',[-20,20],'mappedTracksField','associatedMTP1');
        buildFiberManifoldAndMapMT(P2,[kinTracksISOInliers(tIdx) randTracksCell{:}],EB3TracksInliers,mappingDist,'kinDistCutoff',[-20,20],'mappedTracksField','associatedMTP2');
        [~,hFigMapped,zScores(1,tIdx)]=bundleStatistics(MD,'kinBundle',[{kinTracksISOInliers(tIdx)} {[randTracksCell{:}]}],'plotName',['manifMC-P1-tr' num2str(tIdx)],'mappedMTField','associatedMTP1');
        close(hFigMapped);
        [~,hFigMapped,zScores(2,tIdx)]=bundleStatistics(MD,'kinBundle',[{kinTracksISOInliers(tIdx)} {[randTracksCell{:}]}],'plotName',['manifMC-P2-tr' num2str(tIdx)],'mappedMTField','associatedMTP2');
        close(hFigMapped);
    end
    %%
    processScoring=ExternalProcess(MD,'manifoldScoring');
    mkdirRobust([MD.outputDirectory_ filesep 'Kin' filesep 'bundles' ]);
    save([MD.outputDirectory_ filesep 'Kin' filesep 'bundles' filesep 'bundleStats.mat'],'zScores','mappingDist','kinTracksISOInliers')
    processScoring.setOutFilePaths({[MD.outputDirectory_ filesep 'Kin' filesep 'bundles' filesep 'bundleStats.mat']});
end

MD.getPackage(333).setProcess(7,processScoring);


%% Display the N best and worst scores
tmp=load(processScoring.outFilePaths_{1});
zScores=tmp.zScores;
mappingDist=tmp.mappingDist;  

kinTracksISOInliers=tmp.kinTracksISOInliers;
%%
[sortedScore,indx]=sort(zScores(:));

% suppressing the unprocessed manifold
indx(isnan(sortedScore(:)))=[];
sortedScore(isnan(sortedScore(:)))=[];
[poleIdx,kinIdx]=ind2sub(size(zScores),indx);
%%
N=1;
tic
myColormap=uint8( ...
    [[0 255 255]; ... % blue "all" tracks
    [0 255 0]; ... % green mapped tracks
    [255 0 100]; ... % kinetochore tracks
    [255 0 200]; ... % kinetochore tracks
    ]);

for scoreIdx=[1:N (length(sortedScore)-N):length(sortedScore)]
    tIdx=kinIdx(scoreIdx);
    track=kinTracksISOInliers(tIdx);
    refP1P2=buildRefsFromTracks(P1,P2);
    if(poleIdx(scoreIdx)==1)
        ROI=[P1 track];
        processProj=ExternalProcess(MD,'projP1');
        project1D(MD,ROI,'dynPoligonREF',refP1P2.applyBase([P2 ROI],''),'FoF',refP1P2, ...
            'name',['SP-P1-' num2str(tIdx) '-S-' num2str(sortedScore(scoreIdx))], ...
            'channelRender','grayRed','processSingleProj',processProj, ... 
            'intMinPrctil',[20 98],'intMaxPrctil',[99.9 100],'fringeWidth',30,'insetFringeWidth',mappingDist);
        overlayProjTracksMovie(processProj,'tracks',[refP1P2.applyBase([track.associatedMTP1; track],[]) ] , ... 
            'colorIndx',[ones(size(track.associatedMTP1)); 3],'colormap',myColormap,'name','test');
        
        refKP=buildRefsFromTracks(P1,track);
        project1D(MD,ROI,'dynPoligonREF',refKP.applyBase([ROI],''),'FoF',refKP, ...
            'name',['PK-P1-' num2str(tIdx) '-S-' num2str(sortedScore(scoreIdx))], ... 
            'channelRender','grayRed','processSingleProj',processProj, ... 
            'intMinPrctil',[20 98],'intMaxPrctil',[99.9 100],'fringeWidth',30,'insetFringeWidth',mappingDist);
        overlayProjTracksMovie(processProj,'tracks',refKP.applyBase([track.associatedMTP1; track],[]), ... 
            'colorIndx',[ones(size(track.associatedMTP1)); 3],'colormap',myColormap,'name','tracks');
    end
    
    if(poleIdx(scoreIdx)==2)
        ROI=[P2 track];
        processProj=ExternalProcess(MD,'projP2');
        project1D(MD,ROI,'dynPoligonREF',refP1P2.applyBase([P1 ROI],''),'FoF',refP1P2, ...
            'name',['SP-P2-' num2str(tIdx) '-S-' num2str(sortedScore(scoreIdx))], ...
            'channelRender','grayRed','processSingleProj',processProj, ...
            'intMinPrctil',[20 98],'intMaxPrctil',[99.9 100],'fringeWidth',30,'insetFringeWidth',mappingDist);
        overlayProjTracksMovie(processProj,'tracks',refP1P2.applyBase([track.associatedMTP2; track],[]), ...
            'colorIndx',[ones(size(track.associatedMTP2)); 3],'colormap',myColormap,'name','test');

       refKP=buildRefsFromTracks(P2,track);
        project1D(MD,ROI,'dynPoligonREF',refKP.applyBase([ROI],''),'FoF',refKP, ...
            'name',['PK-P2-' num2str(tIdx) '-S-' num2str(sortedScore(scoreIdx))], ...
            'channelRender','grayRed','processSingleProj',processProj, ...
            'intMinPrctil',[20 98],'intMaxPrctil',[99.9 100],'fringeWidth',30,'insetFringeWidth',mappingDist);
        overlayProjTracksMovie(processProj,'tracks',refKP.applyBase([track.associatedMTP2; track],[]) , ... 
            'colorIndx',[ones(size(track.associatedMTP2)); 3],'colormap',myColormap,'name','tracks');
    
    end
end

package=GenericPackage({processDetectEB3    processTrackEB3     processDetectKT processTrackKT ... 
                        processDetectPoles  processBuildRef     processScoring  });

