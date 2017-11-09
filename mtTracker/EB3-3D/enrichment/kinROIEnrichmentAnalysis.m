function [dynROIData,overlayCell]=kinROIEnrichmentAnalysis(MD,varargin)
    ip = inputParser;
    ip.CaseSensitive = false;
    ip.KeepUnmatched = true;
    ip.addParameter('package',[]);
    ip.addParameter('debug',[]);
    ip.addParameter('packPID',850);
    ip.addParameter('dynROIData',[]);   % Side loading projection
    ip.addParameter('KT',[]);
    ip.addParameter('printManifCount',[]);
    ip.parse(varargin{:});
    p=ip.Results;

% Process type placeholdes
packPID=p.packPID;
packPIDTMP=packPID+1;
projPackPID=packPID+2;
overlayPackPID=packPID+3;
% We save the temp processes at packID+1, that way if process get killed
% halfway through, or other issues, the previous logs are not overidden but
% the last computation location are still available.

% Process type placeholder (they should be defined by tag and then
% retrieved below)
% Process type placeholdes
MD.setPackage(packPIDTMP,GenericPackage({ ... 
    PointSourceDetectionProcess3D(MD,[MD.outputDirectory_ filesep 'EB3']),...
    TrackingProcess(MD,[MD.outputDirectory_ filesep 'EB3']),...
    PointSourceDetectionProcess3D(MD,[MD.outputDirectory_ filesep 'KT']),...
    TrackingProcess(MD,[MD.outputDirectory_ filesep 'KT']),...  
    ExternalProcess(MD,'detectPoles'),...
    ExternalProcess(MD,'buildRefsAndROI'),...
    ExternalProcess(MD,'project1D'),...   
    ExternalProcess(MD,'manifoldScoring'),...
    ExternalProcess(MD,'project1D'),...
    ExternalProcess(MD,'builDynROIs'),...
    ExternalProcess(MD,'overlayProjDetectionMovie'),...    
  },[],'fiberTrackabilityAnalysis'));


lpid=0;
lpid=lpid+1;                    
if(GenericPackage.processExist(p.package,lpid))
   processDetectEB3=p.package.getProcess(lpid);
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
if(GenericPackage.processExist(p.package,lpid))
    processTrackEB3=p.package.getProcess(lpid);
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
if(GenericPackage.processExist(p.package,lpid))
   processDetectKT=p.package.getProcess(lpid);
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
if(GenericPackage.processExist(p.package,lpid))
    processTrackKT=p.package.getProcess(lpid);
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
if(GenericPackage.processExist(p.package,lpid))
    processDetectPoles=p.package.getProcess(lpid);
else
    processDetectPoles=ExternalProcess(MD,'detectPoles',@(p) detectPoles(p.getOwner(),'process',p,'isoOutput',true));
    processDetectPoles.run();
end

MD.getPackage(packPIDTMP).setProcess(lpid,processDetectPoles);

lpid=lpid+1;
if(GenericPackage.processExist(p.package,lpid))
    processBuildRef=p.package.getProcess(lpid);
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
if(GenericPackage.processExist(p.package,lpid))
    processScoring=p.package.getProcess(lpid);
else
    lftThreshold=4;
    %    lftThreshold=2;

    kinTest=find([kinTracksISOInliers.lifetime]>lftThreshold);
    maxRandomDist=20;
    mappingDist=10;
    processManifoldAntispace=ExternalProcess(MD,'randomizeTracks');
    zScores=nan(2,length(kinTracksISOInliers));

    for testIdx=1:length(kinTest)
        tIdx=kinTest(testIdx);
        [randTracksCell]=randomizeTracksMC(MD,maxRandomDist,'randomType','manifoldAntispace',...
            'tracks',kinTracksISOInliers(tIdx),'process',processManifoldAntispace,'simuNumber',200);
        buildFiberManifoldAndDetections(P1,[kinTracksISOInliers(tIdx) randTracksCell{:}],EB3TracksInliers,mappingDist,'kinDistCutoff',[-20,20],'mappedTracksField','associatedMTP1');
        buildFiberManifoldAndDetections(P2,[kinTracksISOInliers(tIdx) randTracksCell{:}],EB3TracksInliers,mappingDist,'kinDistCutoff',[-20,20],'mappedTracksField','associatedMTP2');
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
if(GenericPackage.processExist(p.package,lpid))
    processProjSpindleRef=p.package.getProcess(lpid);
else
    processProjSpindleRef=ExternalProcess(MD,'project1D');
    project1D(  MD,[P1,P2],'FoF',refs(1,2), ...
        'name','CroppedSpindleRef','channelRender','grayRed', 'insetFringeWidth',80, ...
        'processSingleProj',processProjSpindleRef, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);
end
MD.getPackage(packPIDTMP).setProcess(lpid,processProjSpindleRef);


%% Display the N best and worst scores

% spindle frame of reference for detection    cellfun(@(p) disp(p.name_),MD.getPackage(150).processes_)
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

disp('Producing MIPs');
tic
if(isempty(p.dynROIData))
    dynROIData=cell(2,length(kinTracksISOInliers));
else
    dynROIData=p.dynROIData;
end

for KTIndex=1:length(projKin)
    tIdx=projKin(KTIndex);
    pIdx=projPole(KTIndex);
    if(isempty(dynROIData{pIdx,tIdx}))
        track=kinTracksISOInliers(projKin(KTIndex));
        PTrack=fiducials(projPole(KTIndex));
        mappedTracks={track.associatedMTP1 track.associatedMTP2};
        mappedTracksCell{KTIndex}=mappedTracks(projPole(KTIndex));
        %refP1P2=buildRefsFromTracks(P1,P2);
        ROI=[PTrack track];
        processProj=ExternalProcess(MD,'dynROIProj');
        refKP=buildRefsFromTracks(PTrack,track);
        project1D(MD,ROI,'dynPoligonREF',refKP.applyBase([ROI],''),'FoF',refKP, ...
            'name',['PK-P' num2str(projPole(KTIndex)) '-' num2str(tIdx)], ...
            'channelRender','grayRed','processSingleProj',processProj, ...
            'intMinPrctil',[20 98],'intMaxPrctil',[100 100],'fringeWidth',30,'insetFringeWidth',mappingDist);
        ldynROIData.processProjCell=processProj;
        ldynROIData.refsCell=refKP;
        ldynROIData.ROIsCell=ROI;
        ldynROIData.mappedTracksCell=mappedTracksCell;
        dynROIData{projPole(KTIndex),tIdx}=ldynROIData;
        processDynROI=ExternalProcess(MD,'builDynROIs');
        save([MD.outputDirectory_ filesep 'dynROI' filesep 'dynROIData.mat'],'dynROIData');
        processDynROI.setOutFilePaths([MD.outputDirectory_ filesep 'dynROI' filesep 'dynROIData.mat']);
    end
end
MD.setPackage(packPIDTMP).setProcess(lpid,processDynROI);

% systematically, validating measure accross scale. 

disp('Overlay hard-code ROI selected detection');
tic;
lpid=lpid+1;
overlayCell=cell(2,length(projKin));
for KTIndex=1:length(projKin)
    tIdx=projKin(KTIndex);
    pIdx=projPole(KTIndex);

    processSelectedDetectionOverlay=ExternalProcess(MD,'overlayProjDetectionMovie');
    myColormap=255*jet(256);
    selectedDetection=mapDetectionsTo1DManifold(dynROIData{pIdx,tIdx}.ROIsCell,oDetections,mappingDist);
    overlayProjDetectionMovie(dynROIData{pIdx,tIdx}.processProjCell,'detections', dynROIData{pIdx,tIdx}.refsCell.applyBase(selectedDetection,''), ... 
        'colorIndx',arrayfun(@(d) ceil(100*ones(length(d.xCoord))), selectedDetection,'unif',0), ...
        'process',processSelectedDetectionOverlay, ...
        'colormap',myColormap,'name','Red')
    overlayCell{pIdx,tIdx}=processSelectedDetectionOverlay;
end;
toc;
MD.getPackage(packPIDTMP).setProcess(lpid,processSelectedDetectionOverlay);

% disp('Build MIP')
% tic;
% printProcMIPArray(overlayCell,    ... 
%                    [fileparts(fileparts(overlayCell.outFilePaths_{4})) filesep 'MIPArray' filesep], ...
%                     'MIPIndex',2,'MIPSize',400,'maxHeight',1500,'maxWidth',1000);
% toc;

MD.setPackage(packPID,MD.getPackage(packPIDTMP));

