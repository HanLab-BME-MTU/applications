function [processOverlayCells]=trackKT(MD,varargin)
    ip = inputParser;
    ip.CaseSensitive = false;
    ip.KeepUnmatched = true;
    ip.addParameter('package',[]);
    ip.addParameter('debug',false);
    ip.addParameter('packPID',1600);
    ip.addParameter('name','');
    ip.addParameter('dynROIView','');
    ip.addParameter('forceRunIdx',[]);
    ip.addParameter('dynROIData',[]);   % Side loading projection
    ip.addParameter('KT',[]);
    ip.addParameter('printManifCount',[]);
    ip.parse(varargin{:});
    p=ip.Results;

% Process type placeholdes
packPID=p.packPID;
% We save the temp processes at packID+1, that way if process get killed
% halfway through, or other issues, the previous logs are not overidden but
% the last computation location are still available.

% Process type placeholder (they should be defined by tag and then
% retrieved below)
% Process type placeholdes
templatePackage=GenericPackage({ ... 
    PointSourceDetectionProcess3D(MD,[MD.outputDirectory_ filesep 'KT']),...
    PointSourceDetectionProcess3D(MD,[MD.outputDirectory_ filesep 'KT' filesep 'spindleRef']),...
    TrackingProcess(MD,[MD.outputDirectory_ filesep 'KT']),...  
    ExternalProcess(MD,'projDynROI'),...  
    ExternalProcess(MD,'overlayProjDetectionMovie'),...  
    ExternalProcess(MD,'overlayProjacksMovie')... 
    },[],'name_',['trackKT' p.name '_backup']);


packPIDTMP=packPID+1;
package=p.package;

MD.setPackage(packPIDTMP,templatePackage);

pack=MD.searchPackageName('buildAndProjectSpindleRef');
processDetectPoles=pack.getProcess(1);
processBuildRef=pack.getProcess(2);
processProj=pack.getProcess(3);

%% loading fiducials and references (should be dynROI with associated ref really...)
tmp=load(processDetectPoles.outFilePaths_{1});
P1=tmp.tracks(1);
P2=tmp.tracks(2);
fiducials=[P1 P2];
refs=load(processBuildRef.outFilePaths_{1}); refs=refs.refs;
ROIs=load(processBuildRef.outFilePaths_{2}); ROIs=ROIs.ROI;

lpid=0;

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

%% register tracks in the spindle frame of reference
tmp=load(processDetectKT.outFilePaths_{2}); detection=tmp.movieInfo;
oDetections=Detections(detection);
movieInfoOrig=detection;

oDetectionsP1P2=refs(1,2).applyBase(oDetections,'');

%% Create fack process for detections, shifted to deal with U-track BS.
detCopy=copy(oDetectionsP1P2);
detCopy.addOffset(100,100,100);
movieInfo=detCopy.getStruct();
outputFolder=fullfile(fileparts(processDetectKT.outFilePaths_{2}),'spindleRef'); mkdirRobust(outputFolder);
outputFilePath=fullfile(outputFolder,'channel_2.mat'); 
save(outputFilePath,'movieInfo');
processDetectKTRef=PointSourceDetectionProcess3D(MD, outputFolder,UTrackPackage3D.getDefaultDetectionParams(MD,outputFolder));
processDetectKTRef.setOutFilePaths({[] ,outputFilePath});
MD.addProcess(processDetectKTRef);  

lpid=lpid+1;                    
MD.getPackage(packPIDTMP).setProcess(lpid,processDetectKTRef);

%%
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
    paramsIn.DetProcessIndex=processDetectKTRef.getIndex();
    processTrackKT.run(paramsIn);
end
MD.getPackage(packPIDTMP).setProcess(lpid,processTrackKT);

%% Setup track back in the normal referential
tmp=load(processTrackKT.outFilePaths_{2}); 
tracksFinalStripped=rmfield(tmp.tracksFinal,'tracksCoordAmpCG');
kinTracksISO=TracksHandle(tracksFinalStripped,movieInfoOrig); 
tracksFinal=kinTracksISO.getStruct()';
save(processTrackKT.outFilePaths_{2},'tracksFinal');

%% Mapping tracks
kinTracksISOInliers=mapTracksTo1DManifold(ROIs{1,2},kinTracksISO,0,'position','start','distType','vertexDistOtsu');

if(p.debug)
%% Display CH1 only
lpid=lpid+1;     
if(~isempty(p.package)&&(length(p.package.processes_)>=lpid)&&(~isempty(p.package.getProcess(lpid))))
    processProjSpindleRef=p.package.getProcess(lpid);
else
    processPro=ExternalProcess(MD,'projectDynROI');
    processProjSpindleRef=ExternalProcess(MD,'projectDynROI');
    projectDynROI(  MD,[P1,P2],'FoF',refs(1,2),'renderedChannel',[2], ...
        'name','CroppedSpindleRef','channelRender','grayRed', 'insetFringeWidth',80, ...
        'processSingleProj',processPro,'processRenderer',processProjSpindleRef, ...
        'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);
end
MD.getPackage(packPIDTMP).setProcess(lpid,processProjSpindleRef);

%% Overlay data
trs=kinTracksISO;
colorIndx=arrayfun(@(d) zeros(length(d.xCoord(:,1)),1),oDetections,'unif',0); 
tracksColorIndx=randperm(length(trs));
for fIdx=1:trs.numTimePoints()
    for tIdx=1:length(trs)
        if(trs(tIdx).tracksFeatIndxCG(trs(tIdx).f==fIdx)>0)
            colorIndx{fIdx}(trs(tIdx).tracksFeatIndxCG(trs(tIdx).f==fIdx))=tracksColorIndx(tIdx);
        end
    end
end

colorIndx=cellfun(@(c) ceil(255*mat2gray(c,[1,length(trs)]))+1,colorIndx,'unif',0);
colorIndx=cellfun(@(c) 100*ones(size(c)),colorIndx,'unif',0);

lpid=lpid+1;
if(GenericPackage.processExist(package,lpid)&&(~any(lpid==p.forceRunIdx)))
    processAllDetectOverlay=package.getProcess(lpid);
else
    disp('Overlay KT detections');tic;
    processAllDetectOverlay=ExternalProcess(MD,'overlayProjDetectionMovie');
    myColormap=255*jet(256);
    overlayProjDetectionMovie(processProjSpindleRef,'detections', oDetectionsP1P2 , ... 
        'colorIndx',colorIndx, ...
        'colormap',myColormap,'name',['allDets' p.name],'process',processAllDetectOverlay);
    toc;
end
MD.getPackage(packPIDTMP).setProcess(lpid,processAllDetectOverlay);


tic;
myColormap=255*jet(256);
lpid=lpid+1;
if(GenericPackage.processExist(package,lpid)&&(~any(lpid==p.forceRunIdx)))
    process=package.getProcess(lpid);
else
    disp('Overlay KT tracks');
    process=ExternalProcess(MD,'overlayProjTracksMovie');
    overlayProjTracksMovie(processAllDetectOverlay,'tracks', refs(1,2).applyBase(kinTracksISO,''), ... 
        'colorIndx',ceil(255*mat2gray([kinTracksISO.lifetime]',[1 150]))+1,'dragonTail',10,'colormap',myColormap,'name',['allTracks' p.name],'process',process);
end
MD.getPackage(packPIDTMP).setProcess(lpid,process);

%% Isolating the first 3 KTs that do not have a full lifetime
lft=[kinTracksISOInliers.lifetime];
cutKT=(lft<120);
KTIindices=find(cutKT);
KTIindices=KTIindices(1:3);
processDynROIPack=p.dynROIView;
if(isempty(processDynROIPack))
    processDynROICells=cell(1,length(KTIindices));
    for KTIdx=1:length(KTIindices)
        tr=kinTracksISOInliers(KTIindices(KTIdx));

        processProj=ExternalProcess(MD,'projectDynROI');
        projectDynROI(MD,[],tr,'FoF',refs(1,2), ...
            'name',['SingleTrack-' num2str(KTIindices(KTIdx)) p.name],'renderedChannel',2,'channelRender','grayRed', 'fringeWidth',20,'insetFringeWidth',20, ...
            'processRenderer',processProj, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);
        processDynROICells{KTIdx}=processProj;
    end
    processDynROIPack=GenericPackage(processDynROICells,[],'name_','dynROIView');
    MD.addPackage(processDynROIPack);
end

if(nargout>0)
    processOverlayCells=cell(1,length(KTIindices));
    for KTIdx=1:length(KTIindices)
        processProj=processDynROIPack.getProcess(KTIdx);
        process=ExternalProcess(MD,'overlayProjTracksMovie');
        processAllDetectOverlay=ExternalProcess(MD,'overlayProjDetectionMovie');
        overlayProjDetectionMovie(processProj,'detections', oDetectionsP1P2 , ... 
            'colorIndx',colorIndx, ...
            'colormap',myColormap,'name',['localDets' p.name],'process',processAllDetectOverlay);
        overlayProjTracksMovie(processAllDetectOverlay,'tracks', kinTracksISOInliers , ... 
            'colorIndx',ceil(255*mat2gray([kinTracksISOInliers.lifetime]',[1 150]))+1,'dragonTail',10,'colormap',myColormap,'name','localTrack','process',process);
        processOverlayCells{KTIdx}=process;
    end
end
end
completePackage=GenericPackage(MD.getPackage(packPIDTMP).processes_,[],'name_',['trackKT' p.name]);
MD.setPackage(packPID,completePackage);

