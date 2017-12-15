function [processOverlayCells]=trackGeneric(MD,varargin)
    ip = inputParser;
    ip.CaseSensitive = false;
    ip.KeepUnmatched = true;
    ip.addParameter('package',[]);
    ip.addParameter('debug',false);
    ip.addParameter('createROIs',false);
    ip.addParameter('packPID',10);
    ip.addParameter('name','');
    ip.addParameter('dynROIView',[]);
    ip.addParameter('forceRunIdx',[]);
    %ip.addParameter('dynROIData',[]);   % Side loading projection
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
if(p.createROIs)
    ROIName='trackingROI'
    processROIs=ExternalProcess(MD,'createROIs');
    [~,mask]=createROI(MD);
    MD=crop3D(MD,mask,'name',ROIName);
end

templatePackage=GenericPackage({ ... 
    PointSourceDetectionProcess3D(MD,[MD.outputDirectory_ ]),...
    TrackingProcess(MD,[MD.outputDirectory_ ]),...  
    ExternalProcess(MD,'projDynROI'),...  
    ExternalProcess(MD,'overlayProjDetectionMovie'),...  
    ExternalProcess(MD,'overlayProjacksMovie')... 
    },[],'name_',['trackGeneric' p.name '_backup']);


packPIDTMP=packPID+1;
package=p.package;

MD.setPackage(packPIDTMP,templatePackage);


lpid=0;

lpid=lpid+1;                    
if(GenericPackage.processExist(package,lpid)&&(~any(lpid==p.forceRunIdx)))
 processDetect=package.getProcess(lpid);
else
    processDetect=PointSourceDetectionProcess3D(MD,[MD.outputDirectory_ filesep ],UTrackPackage3D.getDefaultDetectionParams(MD, [MD.outputDirectory_ filesep ]));
    MD.addProcess(processDetect);
    funParams = processDetect.funParams_;
    funParams.showAll=true;
    funParams.alpha=0.05;
    funParams.filterSigma=[1.8; .75];
    funParams.WindowSize={[]};
    funParams.ClearMaskBorder=false;
    funParams.algorithmType= {'pointSourceFit'}
    funParams.ConfRadius={[]};       
    processDetect.setPara(funParams);
    paramsIn.ChannelIndex=1;
    paramsIn.isoCoord=true;
    processDetect.run(paramsIn);    
end
MD.getPackage(packPIDTMP).setProcess(lpid,processDetect);

disp('loading detections')
tmp=load(processDetect.outFilePaths_{1}); detection=tmp.movieInfo;
oDetections=Detections(detection);


disp('Computing trackability');
[maxSpeed,~,densities]=estimateTrackability(oDetections,1/MD.timeInterval_,'debugMode','');
histogram(vertcat(maxSpeed{:}));
xlabel('Max. speed in pixel');

lpid=lpid+1;                    
if(GenericPackage.processExist(package,lpid)&&(~any(lpid==p.forceRunIdx)))
    processTrack=package.getProcess(lpid);
else
    processTrack=TrackingProcess(MD, [MD.outputDirectory_],UTrackPackage3D.getDefaultTrackingParams(MD,[MD.outputDirectory_ ]));
    MD.addProcess(processTrack);    
    funParams = processTrack.funParams_;
    [gapCloseParam,costMatrices,kalmanFunctions,probDim,verbose]=genericTrackingParam();
    funParams.gapCloseParam=gapCloseParam;
    funParams.costMatrices=costMatrices;
    funParams.kalmanFunctions=kalmanFunctions;
    funParams.probDim=probDim;
    processTrack.setPara(funParams);
    paramsIn.ChannelIndex=1;
    paramsIn.DetProcessIndex=processDetect.getIndex();
    processTrack.run(paramsIn);
end
MD.getPackage(packPIDTMP).setProcess(lpid,processTrack);

%% Setup track back in the normal referential
disp('loading tracks')
tmp=load(processTrack.outFilePaths_{1}); 
tracks=TracksHandle(tmp.tracksFinal); 

%% Tracks in the lab FoR.
disp('Amira write')
% amiraWriteTracks([fileparts(processTrack.outFilePaths_{1}) filesep 'Amira' filesep 'AmiraTrackLabRefDT10' filesep 'tracksLabRef.am'],tracks,'dragonTail',10);
% amiraWriteTracks([fileparts(processTrack.outFilePaths_{1}) filesep 'Amira' filesep 'AmiraTrackLabRef' filesep 'tracksLabRef.am'],tracks);
% amiraWriteTracks([fileparts(processTrack.outFilePaths_{1}) filesep 'Amira' filesep  'AmiraTrackLabRef20plus' filesep 'trackLabRef20plus.am'],tracks([tracks.lifetime]>20));

if(p.debug)
%% Display CH1 only
lpid=lpid+1;     
if(GenericPackage.processExist(p.package,lpid))
    processRenderer=p.package.getProcess(lpid);
else
    processProj=ExternalProcess(MD,'projectDynROI');
    processRenderer=ExternalProcess(MD,'projectDynROI');
    projectDynROI(  MD,'renderedChannel',[1], ...
        'name','FullCell','channelRender','grayRed', 'insetFringeWidth',80, ...
        'processSingleProj',processProj,'processRenderer',processRenderer, ...
        'intMinPrctil',[1],'intMaxPrctil',[100]);
end
MD.getPackage(packPIDTMP).setProcess(lpid,processRenderer);

%% Build detection color according to tracks ID
trs=tracks;
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
%colorIndx=cellfun(@(m) ceil(254*mat2gray(m,prctile(vertcat(maxSpeed{:}),[0.05 0.95])))+1,maxSpeed,'unif',0);

lpid=lpid+1;
if(GenericPackage.processExist(package,lpid)&&(~any(lpid==p.forceRunIdx)))
    processAllDetectOverlay=package.getProcess(lpid);
else
    disp('Overlay detections');tic;
    processAllDetectOverlay=ExternalProcess(MD,'overlayProjDetectionMovie');
    myColormap=255*jet(256);
    overlayProjDetectionMovie(processRenderer,'detections', oDetections , ... 
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
    disp('Overlay tracks');
    process=ExternalProcess(MD,'overlayProjTracksMovie');
    overlayProjTracksMovie(processAllDetectOverlay,'tracks', tracks, ... 
        'colorIndx',ceil(255*mat2gray([tracks.lifetime]',[1 MD.nFrames_]))+1,'dragonTail',10,'colormap',myColormap,'name',['allTracks' p.name],'process',process);
end
MD.getPackage(packPIDTMP).setProcess(lpid,process);

%% Isolating the first 3 tracks that do not have a full lifetime
if(nargout>0)
    processDynROIPack=p.dynROIView;


    if(isempty(processDynROIPack))
        KTIindices=1;
        processDynROICells=cell(1,length(KTIindices));
        for KTIdx=1:length(KTIindices)
            tr=tracks(KTIindices(KTIdx));
            tr=tr.copy();
            tr.x(:)=tr.x(1);
            tr.y(:)=tr.y(1);
            tr.z(:)=tr.z(1);
            ref=FrameOfRef().setOriginFromTrack(tr).genCanonicalBase();
            processProj=ExternalProcess(MD,'projectDynROI');
            processRendererROI=ExternalProcess(MD,'projectDynROI');
            projectDynROI(MD,[],ref.applyBase(tr,''),'FoF',ref, ...
                'name',['SingleTrack-' num2str(KTIindices(KTIdx)) p.name],'renderedChannel',1,'channelRender','grayRed', 'fringeWidth',30,'insetFringeWidth',30, ...
                'processSingleProj',processProj,'processRenderer',processRendererROI, 'intMinPrctil',[1],'intMaxPrctil',[100]);
            processDynROICells{KTIdx}=processRendererROI;
        end
        processDynROIPack=GenericPackage(processDynROICells,[],'name_','dynROIView');
        MD.addPackage(processDynROIPack);
    end

    processOverlayCells=cell(1,length(processDynROIPack.processes_));
    for KTIdx=1:length(processDynROIPack.processes_)
        processRendererROI=processDynROIPack.getProcess(KTIdx);
        
        processAllDetectOverlay=ExternalProcess(MD,'overlayProjDetectionMovie');
        colorIndx=arrayfun(@(d) ceil(255*mat2gray(d.yCoord(:,1)))+1,oDetections,'unif',0);
        colorIndx=cellfun(@(s) ceil(255*mat2gray(s,[0,3]))+1,maxSpeed,'unif',0);
        overlayProjDetectionMovie(processRendererROI,'detections', oDetections , ... 
            'colorIndx',colorIndx, 'colormap',myColormap, ... 
            'name',['localDets' p.name],'process',processAllDetectOverlay,'radius',2);
        colorIndx=ceil(255*mat2gray([tracks.lifetime]',[1 MD.nFrames_]))+1;

        process=ExternalProcess(MD,'overlayProjTracksMovie');
        colorIndx=ceil(255*mat2gray(randi(length(tracks),1,length(tracks))))+1;
        overlayProjTracksMovie(processAllDetectOverlay,'tracks', tracks , ... 
             'colorIndx',colorIndx,'dragonTail',20,'colormap',myColormap,'name','localTrack','process',process);
        processOverlayCells{KTIdx}=process;
    end
end
end
completePackage=GenericPackage(MD.getPackage(packPIDTMP).processes_,[],'name_',['trackGeneric' p.name]);
MD.setPackage(packPID,completePackage);

