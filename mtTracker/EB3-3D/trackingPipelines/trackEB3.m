function [overlayCell]=trackEB3(MD,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('package',[]);
ip.addParameter('debug',false);
ip.addParameter('packPID',1700);
ip.addParameter('name','');
ip.addParameter('dynROIView','');
ip.addParameter('forceRunIdx',[]);
ip.addParameter('dynROIData',[]);   % Side loading projection
ip.addParameter('KT',[]);
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
        PointSourceDetectionProcess3D(MD,[MD.outputDirectory_ filesep 'EB3']),...
        PointSourceDetectionProcess3D(MD,[MD.outputDirectory_ filesep 'EB3' filesep 'spindleRef']),...
        TrackingProcess(MD,[MD.outputDirectory_ filesep 'EB3']),...  
        ExternalProcess(MD,'projDynROI'),...  
        ExternalProcess(MD,'overlayProjDetectionMovie'),...  
        ExternalProcess(MD,'overlayProjacksMovie')... 
        },[],'name_',['trackEB3' p.name '_backup']);


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

%% register tracks in the spindle frame of reference
tmp=load(processDetectEB3.outFilePaths_{1}); detection=tmp.movieInfo;
oDetections=Detections(detection);
movieInfoOrig=detection;

oDetectionsP1P2=refs(1,2).applyBase(oDetections,'');

%% Create fack process for detections, shifted to deal with U-track BS.
detCopy=copy(oDetectionsP1P2);
detCopy.addOffset(100,100,100);
movieInfo=detCopy.getStruct();
outputFolder=fullfile(fileparts(processDetectEB3.outFilePaths_{1}),'spindleRef'); mkdirRobust(outputFolder);
outputFilePath=fullfile(outputFolder,'channel_1.mat'); 
save(outputFilePath,'movieInfo');
processDetectEB3Ref=PointSourceDetectionProcess3D(MD, outputFolder,UTrackPackage3D.getDefaultDetectionParams(MD,outputFolder));
processDetectEB3Ref.setOutFilePaths({outputFilePath,[]});
MD.addProcess(processDetectEB3Ref);  

lpid=lpid+1;                    
MD.getPackage(packPIDTMP).setProcess(lpid,processDetectEB3Ref);

lpid=lpid+1;                    
if(GenericPackage.processExist(package,lpid)&&(~any(lpid==p.forceRunIdx)))
    processTrackEB3=package.getProcess(lpid);
else
    processTrackEB3=TrackingProcess(MD, [MD.outputDirectory_ filesep 'EB3'],UTrackPackage3D.getDefaultTrackingParams(MD, [MD.outputDirectory_ filesep 'EB3']));
    MD.addProcess(processTrackEB3);    
    funParams = processTrackEB3.funParams_;
    [costMatrices,gapCloseParam,kalmanFunctions,probDim]=plusTipCometTracker3DParamDebug(MD);
    funParams.gapCloseParam=gapCloseParam;
    funParams.costMatrices=costMatrices;
    funParams.kalmanFunctions=kalmanFunctions;
    funParams.probDim=probDim;
    processTrackEB3.setPara(funParams);
    paramsIn.ChannelIndex=1;
    paramsIn.DetProcessIndex=processDetectEB3Ref.getIndex();
    processTrackEB3.run(paramsIn);
end
MD.getPackage(packPIDTMP).setProcess(lpid,processTrackEB3);

%% Setup track back in the normal referential
tmp=load(processTrackEB3.outFilePaths_{1}); 
tracksFinalStripped=rmfield(tmp.tracksFinal,'tracksCoordAmpCG');
EB3TracksISO=TracksHandle(tracksFinalStripped,movieInfoOrig); 
tracksFinal=EB3TracksISO.getStruct()';
save(processTrackEB3.outFilePaths_{1},'tracksFinal');

%% Mapping tracks
EB3TracksInliers=mapTracksTo1DManifold(ROIs{1,2},EB3TracksISO,80,'position','start','distType','euclideanDist');


%% Tracks in the lab FoR.
%% Tracks in the lab FoR.
amiraWriteTracks([fileparts(processTrackEB3.outFilePaths_{1}) filesep 'Amira' filesep 'AmiraTrackLabRefDT10' filesep 'tracksLabRef.am'],EB3TracksInliers,'dragonTail',10);
amiraWriteTracks([fileparts(processTrackEB3.outFilePaths_{1}) filesep 'Amira' filesep 'AmiraTrackLabRef' filesep 'tracksLabRef.am'],EB3TracksInliers);
amiraWriteTracks([fileparts(processTrackEB3.outFilePaths_{1}) filesep 'Amira' filesep  'AmiraTrackLabRef20plus' filesep 'trackLabRef20plus.am'],EB3TracksInliers([EB3TracksInliers.lifetime]>20));

%%
tmp=load(processDetectEB3.outFilePaths_{1}); detection=tmp.movieInfo;
oDetections=Detections(detection);

if(p.debug)
    %% Display CH1 only
    lpid=lpid+1;     
    if(~isempty(p.package)&&(length(p.package.processes_)>=lpid)&&(~isempty(p.package.getProcess(lpid))))
        processProjSpindleRef=p.package.getProcess(lpid);
    else
        processPro=ExternalProcess(MD,'projectDynROI');
        processProjSpindleRef=ExternalProcess(MD,'projectDynROI');
        projectDynROI(  MD,[P1,P2],'FoF',refs(1,2),'renderedChannel',[1 2], ...
            'name','CroppedSpindleRef','channelRender','grayRed', 'insetFringeWidth',80, ...
            'processSingleProj',processPro,'processRenderer',processProjSpindleRef, ...
            'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);
    end
    MD.getPackage(packPIDTMP).setProcess(lpid,processProjSpindleRef);

    %% Overlay data
    trs=EB3TracksISO;
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
        overlayProjTracksMovie(processProjSpindleRef,'tracks', refs(1,2).applyBase(EB3TracksISO,''), ... 
            'colorIndx',ceil(255*mat2gray([EB3TracksISO.lifetime]',[1 20]))+1,'dragonTail',10,'colormap',myColormap,'name',['allTracks' p.name],'process',process);
    end
    MD.getPackage(packPIDTMP).setProcess(lpid,process);
end



completePackage=GenericPackage(MD.getPackage(packPIDTMP).processes_,[],'name_','trackEB3');
MD.setPackage(packPID,completePackage);