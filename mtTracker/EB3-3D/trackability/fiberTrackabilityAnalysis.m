function [overlayCell]=fiberTrackabilityAnalysis(MD,varargin)
    ip = inputParser;
    ip.CaseSensitive = false;
    ip.KeepUnmatched = true;
    ip.addParameter('package',[]);
    ip.addParameter('debug',[]);
    ip.addParameter('packPID',1000);
    ip.addParameter('forceRunIdx',[]);
    ip.addParameter('dynROIData',[]);   % Side loading projection
    ip.addParameter('KT',[]);
    ip.addParameter('printManifCount',[]);
    ip.parse(varargin{:});
    p=ip.Results;

% Process type placeholdes
packPID=p.packPID;
projPackPID=1003;
overlayPackPID=1004;
% We save the temp processes at packID+1, that way if process get killed
% halfway through, or other issues, the previous logs are not overidden but
% the last computation location are still available.

% Process type placeholder (they should be defined by tag and then
% retrieved below)
% Process type placeholdes
templatePackage=GenericPackage({ ... 
    PointSourceDetectionProcess3D(MD,[MD.outputDirectory_ filesep 'EB3']),...
    TrackingProcess(MD,[MD.outputDirectory_ filesep 'EB3']),...
    PointSourceDetectionProcess3D(MD,[MD.outputDirectory_ filesep 'KT']),...
    TrackingProcess(MD,[MD.outputDirectory_ filesep 'KT']),...  
    ExternalProcess(MD,'detectPoles'),...
    ExternalProcess(MD,'buildRefsAndROI'),...
    ExternalProcess(MD,'project1D'),...   
    ExternalProcess(MD,'manifoldScoring'),...
    ExternalProcess(MD,'builDynROIs'),...
    ExternalProcess(MD,'builDynROIs'),...
    ExternalProcess(MD,'overlayProjDetectionMovie'),...  
    ExternalProcess(MD,'builDynROIs'),...
  },[],'name_','fiberTrackabilityAnalysisBackup');


packPIDTMP=packPID+1;
%packs=MD.searchPackageName('fiberTrackabilityAnalysis');
package=p.package;
% if(isempty(package)&&~isempty(packs))
%     package=packs(end-1);
% end

% if(~isempty(packs))
%     MD.setPackage(packPIDTMP,packs(end-1));
% else
MD.setPackage(packPIDTMP,templatePackage);
% end

lpid=0;
lpid=lpid+1;                    
if(GenericPackage.processExist(package,lpid)&&(~any(lpid==p.forceRunIdx)))
   processDetectEB3=package.getProcess(lpid);
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
MD.getPackage(packPIDTMP).setProcess(lpid,processDetectEB3);


lpid=lpid+1;                    
if(GenericPackage.processExist(package,lpid)&&(~any(lpid==p.forceRunIdx)))
    processTrackEB3=package.getProcess(lpid);
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
MD.getPackage(packPIDTMP).setProcess(lpid,processTrackEB3);

lpid=lpid+1;                    
if(GenericPackage.processExist(package,lpid)&&(~any(lpid==p.forceRunIdx)))
   processDetectKT=package.getProcess(lpid);
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
MD.getPackage(packPIDTMP).setProcess(lpid,processDetectKT);


lpid=lpid+1;                    
if(GenericPackage.processExist(package,lpid)&&(~any(lpid==p.forceRunIdx)))
    processTrackKT=package.getProcess(lpid);
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
MD.getPackage(packPIDTMP).setProcess(lpid,processTrackKT);

lpid=lpid+1;
if(GenericPackage.processExist(package,lpid)&&(~any(lpid==p.forceRunIdx)))
    processDetectPoles=package.getProcess(lpid);
else
    processDetectPoles=ExternalProcess(MD,'detectPoles',@(p) detectPoles(p.getOwner(),'process',p,'isoOutput',true));
    processDetectPoles.run();
end

MD.getPackage(packPIDTMP).setProcess(lpid,processDetectPoles);

lpid=lpid+1;
if(GenericPackage.processExist(package,lpid)&&(~any(lpid==p.forceRunIdx)))
    processBuildRef=package.getProcess(lpid);
else
    processBuildRef=ExternalProcess(MD,'buildRefsAndROI',@(p) buildRefsFromTracks(processDetectPoles,processDetectPoles,'buildROI',true,'process',p));
    processBuildRef.run();
end
MD.getPackage(packPIDTMP).setProcess(lpid,processBuildRef);

%%
tmp=load(processDetectPoles.outFilePaths_{1});
P1=tmp.tracks(1);
P2=tmp.tracks(2);
fiducials=[P1 P2];

% inlier tracks
refs=load(processBuildRef.outFilePaths_{1}); refs=refs.refs;
ROIs=load(processBuildRef.outFilePaths_{2}); ROIs=ROIs.ROI;


%% Mapping tracks
tmp=load(processTrackKT.outFilePaths_{2}); kinTracksISO=TracksHandle(tmp.tracksFinal);
kinTracksISOInliers=mapTracksTo1DManifold(ROIs{1,2},kinTracksISO,0,'position','start','distType','vertexDistOtsu');
tmp=load(processTrackEB3.outFilePaths_{1}); EB3TracksISO=TracksHandle(tmp.tracksFinal);
EB3TracksInliers=mapTracksTo1DManifold(ROIs{1,2},EB3TracksISO,0,'position','start','distType','vertexDistOtsu');

%%
tmp=load(processDetectEB3.outFilePaths_{1}); detection=tmp.movieInfo;
oDetections=Detections(detection);

%% Compute scores
lpid=lpid+1;
if(GenericPackage.processExist(package,lpid)&&(~any(lpid==p.forceRunIdx)))
    processScoring=package.getProcess(lpid);
else
    lftThreshold=min(10,MD.nFrames_-1);
    %    lftThreshold=20;

    kinTest=find([kinTracksISOInliers.lifetime]>lftThreshold);
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
MD.getPackage(packPIDTMP).setProcess(lpid,processScoring);


%%
lpid=lpid+1;
if(GenericPackage.processExist(package,lpid)&&(~any(lpid==p.forceRunIdx)))
    processProjSpindleRef=package.getProcess(lpid);
else
    disp('Producing Spindle Ref MIPs');
    processProjSpindleRef=ExternalProcess(MD,'project1D');
    project1D(  MD,[P1,P2],'FoF',refs(1,2), ...
        'name','CroppedSpindleRef','channelRender','grayRed', 'insetFringeWidth',80, ...
        'processSingleProj',processProjSpindleRef, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);
end
MD.getPackage(packPIDTMP).setProcess(lpid,processProjSpindleRef);


% spindle frame of reference for detection  
tmp=load(processScoring.outFilePaths_{1});
zScores=tmp.zScores;
mappingDist=tmp.mappingDist;
kinTracksISOInliers=tmp.kinTracksISOInliers;
[sortedScore,indx]=sort(zScores(:));

if(isempty(p.KT))
    % suppressing the unprocessed manifold
    indx(isnan(sortedScore(:)))=[];
    sortedScore(isnan(sortedScore(:)))=[];
    [poleIdx,kinIdx]=ind2sub(size(zScores),indx);
    N=p.printManifCount;
    projKin=kinIdx([1:N (length(sortedScore)-N+1):length(sortedScore)]);
    projPole=poleIdx([1:N (length(sortedScore)-N+1):length(sortedScore)]);
else
    selectedKT=p.KT(ismember(p.KT,1:length(kinTracksISOInliers)));
    projKin=[selectedKT selectedKT];
    projPole=[ones(size(selectedKT)) 2*ones(size(selectedKT))];
end

disp('Produce fiber ROI associated with labRef');
lpid=lpid+1;
if(GenericPackage.processExist(package,lpid)&&(~any(lpid==p.forceRunIdx)))
    dynROIData=load(package.getProcess(lpid).outFilePaths_{1});
    dynROIData=dynROIData.dynROIData;
else
    dynROIData=cell(2,length(kinTracksISOInliers));
end

for KTIndex=1:length(projKin)
    tIdx=projKin(KTIndex);
    pIdx=projPole(KTIndex);
    if(isempty(dynROIData{pIdx,tIdx}))
        tic;
        track=kinTracksISOInliers(projKin(KTIndex));
        PTrack=fiducials(projPole(KTIndex));
        mappedTracks={track.associatedMTP1 track.associatedMTP2};
        ROI=[PTrack track];
        ref=FrameOfRef().genCanonicalRef(MD.nFrames_);
        ldynROIData.ref=ref;
        ldynROIData.ROI=ROI;
        ldynROIData.mappedTracks=mappedTracks;
        dynROIData{pIdx,tIdx}=ldynROIData;
        toc;
    end
end
processDynROI=ExternalProcess(MD,'builDynROIs');
processDynROI.setProcessTag('builDynROIs_KTPole_Lab');
mkdirRobust([MD.outputDirectory_ filesep 'dynROI' filesep 'builDynROIs_KTPole_Lab' ]);
save([MD.outputDirectory_ filesep 'dynROI' filesep 'builDynROIs_KTPole_Lab' filesep 'dynROIData.mat'],'dynROIData');
processDynROI.setOutFilePaths({[MD.outputDirectory_ filesep 'dynROI' filesep 'builDynROIs_KTPole_Lab' filesep 'dynROIData.mat']});
MD.getPackage(packPIDTMP).setProcess(lpid,processDynROI);

%% The part below is way to computationally intensive to be run routinely



disp('Producing KT fiber MIPs');
lpid=lpid+1;
tic
if(isempty(p.dynROIData))
    if(GenericPackage.processExist(package,lpid)&&(~any(lpid==p.forceRunIdx)))
        dynROIData=load(package.getProcess(lpid).outFilePaths_{1});    
        dynROIData=dynROIData.dynROIData;
    else
        dynROIData=cell(2,length(kinTracksISOInliers));
    end
else
    dynROIData=p.dynROIData;
end
for KTIndex=1:length(projKin)
    tIdx=projKin(KTIndex);
    pIdx=projPole(KTIndex);
    if(isempty(dynROIData{pIdx,tIdx}))
        tic;
        track=kinTracksISOInliers(projKin(KTIndex));
        PTrack=fiducials(projPole(KTIndex));
        mappedTracks={track.associatedMTP1 track.associatedMTP2};
        ROI=[PTrack track];
        processProj=ExternalProcess(MD,'dynROIProj');
        refKP=buildRefsFromTracks(PTrack,track);
        project1D(MD,ROI,'dynPoligonREF',refKP.applyBase([ROI],''),'FoF',refKP, ...
            'name',['PK-P' num2str(projPole(KTIndex)) '-' num2str(tIdx)], ...
            'channelRender','grayRed','processSingleProj',processProj, ...
            'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99],'fringeWidth',30,'insetFringeWidth',mappingDist);
        ldynROIData.processProj=processProj;
        ldynROIData.ref=refKP;
        ldynROIData.ROI=ROI;
        ldynROIData.mappedTracks=mappedTracks;
        dynROIData{pIdx,tIdx}=ldynROIData;
        toc;
    end
end
processDynROI=ExternalProcess(MD,'builDynROIs');
processDynROI.setProcessTag('self');
mkdirRobust([MD.outputDirectory_ filesep 'dynROI' filesep 'self']);
save([MD.outputDirectory_ filesep 'dynROI' filesep 'self' filesep 'dynROIData.mat'],'dynROIData');
processDynROI.setOutFilePaths({[MD.outputDirectory_ filesep 'dynROI' filesep 'self' filesep 'dynROIData.mat']});
MD.getPackage(packPIDTMP).setProcess(lpid,processDynROI);

disp('Computing trackability');
[maxSpeed,~,densities]=estimateTrackability(oDetections,1/MD.timeInterval_,'debugMode','');


disp('Overlay Detection over fiber ROI');
% systematically, validating measure accross scale. 
tic;
lpid=lpid+1;
if(GenericPackage.processExist(package,lpid)&&(~any(lpid==p.forceRunIdx)))
    overlayCell=load(package.getProcess(lpid).outFilePaths_{1});    
    overlayCell=overlayCell.overlayCell;
else
    overlayCell=cell(2,length(kinTracksISOInliers));
end

for KTIndex=1:length(projKin)
    tIdx=projKin(KTIndex);
    pIdx=projPole(KTIndex);
    if(isempty(overlayCell{pIdx,tIdx}))
        processSelectedDetectionOverlay=ExternalProcess(MD,'overlayProjDetectionMovie');
        myColormap=255*flipud(jet(256));
        [selectedDetection,indices]=mapDetectionsTo1DManifold(dynROIData{pIdx,tIdx}.ROI,oDetections,mappingDist);
        colorIndex=cellfun(@(m,i) ceil(254*mat2gray(m(i),[0.5,0.8]))+1,maxSpeed,indices,'unif',0);
        overlayProjDetectionMovie(dynROIData{pIdx,tIdx}.processProj,'detections', dynROIData{pIdx,tIdx}.ref.applyBase(selectedDetection,''), ... 
            'colorIndx',colorIndex, ...
            'process',processSelectedDetectionOverlay, ...
            'colormap',myColormap,'name','trackability')
        overlayCell{pIdx,tIdx}=processSelectedDetectionOverlay;
    end
end;
processOverlay=ExternalProcess(MD,'overlays');
processOverlay.setProcessTag('self');
mkClrDir([MD.outputDirectory_ filesep 'overlays']);
save([MD.outputDirectory_ filesep 'overlays' filesep 'overlays.mat'],'overlayCell');
processOverlay.setOutFilePaths({[MD.outputDirectory_ filesep 'overlays' filesep 'overlays.mat']});
toc;
MD.getPackage(packPIDTMP).setProcess(lpid,processOverlay);

tic
disp('Produce fiber ROI in Lab frame of ref.');
lpid=lpid+1;
if(GenericPackage.processExist(package,lpid)&&(~any(lpid==p.forceRunIdx)))
    dynROIData=load(package.getProcess(lpid).outFilePaths_{1});
    dynROIData=dynROIData.dynROIData;
else
    dynROIData=cell(2,length(kinTracksISOInliers));
end

for KTIndex=1:length(projKin)
    tIdx=projKin(KTIndex);
    pIdx=projPole(KTIndex);
    if(isempty(dynROIData{pIdx,tIdx}))
        tic;
        track=kinTracksISOInliers(projKin(KTIndex));
        PTrack=fiducials(projPole(KTIndex));
        mappedTracks={track.associatedMTP1 track.associatedMTP2};
        ROI=[PTrack track];
        processProj=ExternalProcess(MD,'dynROIProj');
        processVolMask=ExternalProcess(MD,'volMask');
        spindleRef=buildRefsFromTracks(fiducials,fiducials);
        ref=FrameOfRef().genCanonicalRef(MD.nFrames_);
        projectDynROI(MD,ROI, ...
            'name',['Lab' num2str(projPole(KTIndex)) '-' num2str(tIdx)], ...
            'channelRender','grayRed','processRenderer',processProj, ...
            'processMaskVolume',processVolMask,'crop','full', ... 
            'intMinPrctil',[20 98],'intMaxPrctil',[100 100],'fringeWidth',30,'insetFringeWidth',mappingDist);
        ldynROIData.processProj=processProj;    
        ldynROIData.ref=ref;
        ldynROIData.ROI=ROI;
        ldynROIData.mappedTracks=mappedTracks;
        dynROIData{pIdx,tIdx}=ldynROIData;
        toc;
    end
end
processDynROI=ExternalProcess(MD,'builDynROIs');
processDynROI.setProcessTag('P1P2');
mkdirRobust([MD.outputDirectory_ filesep 'dynROI' filesep 'P1P2' ]);
save([MD.outputDirectory_ filesep 'dynROI' filesep 'P1P2' filesep 'dynROIData.mat'],'dynROIData');
processDynROI.setOutFilePaths({[MD.outputDirectory_ filesep 'dynROI' filesep 'P1P2' filesep 'dynROIData.mat']});
MD.getPackage(packPIDTMP).setProcess(lpid,processDynROI);

completePackage=GenericPackage(MD.getPackage(packPIDTMP).processes_,[],'name_','fiberTrackabilityAnalysis');
MD.setPackage(packPID,completePackage);

