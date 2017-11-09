function testProjDynROIMetaphasePlate(MD,varargin)
    ip = inputParser;
    ip.CaseSensitive = false;
    ip.KeepUnmatched = true;
    ip.addParameter('package',[]);
    ip.addParameter('debug',[]);
    ip.addParameter('dynROI','');
    ip.addParameter('packPID',500);
    ip.parse(varargin{:});
    p=ip.Results;

% Process type placeholdes
packPID=p.packPID;
% We save the temp processes at packID+1, that way if process get killed
% halfway through, or other issues, the previous logs are not overidden but
% the last computation location are still available.

% Process type placeholder (they should be defined by tag and then
% retrieved below)
packPIDTMP=packPID+1;
MD.setPackage(packPIDTMP,GenericPackage({ ... 
    ExternalProcess(MD,'project1D'), ...
    ExternalProcess(MD,'detectPoles'),...
    ExternalProcess(MD,'buildRefsAndROI'),...
    ExternalProcess(MD,'project1D'), ... 
    PointSourceDetectionProcess3D(MD,[MD.outputDirectory_ filesep 'KT']),...
    TrackingProcess(MD,[MD.outputDirectory_ filesep 'KT']),...  
    ExternalProcess(MD,'project1D'), ...  
    ComputeMIPProcess(MD), ...
    ExternalProcess(MD,'overlayProjDetectionMovie'), ...    
    ExternalProcess(MD,'overlayProjTrackMovie')
    },'','name_','testProjDynROIMetaphasePlate'));

lpid=0;

lpid=lpid+1;
if(~isempty(p.package)&&(~isempty(p.package.getProcess(lpid))))
    processProj=p.package.getProcess(lpid);
else
    processProj=ExternalProcess(MD,'project1D');
    project1D(  MD, ...
        'name','fullMIPLabFrame','channelRender','grayRed', ...
        'processSingleProj',processProj, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);
end
MD.getPackage(packPIDTMP).setProcess(lpid,processProj);

lpid=lpid+1;
if(~isempty(p.package)&&(~isempty(p.package.getProcess(lpid))))
    processDetectPoles=p.package.getProcess(lpid);
else
    processDetectPoles=ExternalProcess(MD,'detectPoles',@(p) detectPoles(p.getOwner(),'process',p,'isoOutput',true));
    processDetectPoles.run();
end

MD.getPackage(packPIDTMP).setProcess(lpid,processDetectPoles);

lpid=lpid+1;
if GenericPackage.processExist(p.package,lpid)
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

%% 
% inlier tracks
refs=load(processBuildRef.outFilePaths_{1}); refs=refs.refs;
ROIs=load(processBuildRef.outFilePaths_{2}); ROIs=ROIs.ROI;

%%
lpid=lpid+1;      

if GenericPackage.processExist(p.package,lpid)
    processProjSpindleRef=p.package.getProcess(lpid);
else
    processProjSpindleRef=ExternalProcess(MD,'project1D');
    disp('test no  GPU')
    tic;
    project1D(  MD,[P1,P2],'FoF',refs(1,2), ...
        'name','CroppedSpindleRef','channelRender','grayRed', 'insetFringeWidth',80, ...
        'processSingleProj',processProjSpindleRef, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99],'useGPU',false);
    toc
    disp('test GPU')
    tic;
    project1D(  MD,[P1,P2],'FoF',refs(1,2), ...
        'name','CroppedSpindleRef','channelRender','grayRed', 'insetFringeWidth',80, ...
        'processSingleProj',processProjSpindleRef, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99],'useGPU',true);
    toc
end
MD.getPackage(packPIDTMP).setProcess(lpid,processProjSpindleRef);

lpid=lpid+1;
if GenericPackage.processExist(p.package,lpid)
   processDetectKT=p.package.getProcess(lpid);
else
    processDetectKT=PointSourceDetectionProcess3D(MD,[MD.outputDirectory_ filesep 'KT'],UTrackPackage3D.getDefaultDetectionParams(MD, [MD.outputDirectory_ filesep 'KT']));
    MD.addProcess(processDetectKT);
    funParams = processDetectKT.funParams_;
    funParams.showAll=true;
    funParams.alpha=0.05;
    funParams.frameRange=[1 5];
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
if(~isempty(p.package)&&(length(p.package.processes_)>=lpid)&&(~isempty(p.package.getProcess(lpid))))
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

% selectedDetectionP1P2=refs(1,2).applyBase(selectedDetection,'');
% selectedDetectionP1P2.addSphericalCoord();overlayProjDetectionMovie

tmp=load(processTrackKT.outFilePaths_{2}); kinTracksISO=TracksHandle(tmp.tracksFinal);
tmp=load(processDetectKT.outFilePaths_{2}); kinTracksISODetect=Detections(tmp.movieInfo);

oDetections=Detections(kinTracksISODetect);
oDetectionsP1P2=refs(1,2).applyBase(oDetections,'');
oDetectionsP2P1=refs(2,1).applyBase(oDetections,'');
oDetectionsP1P2.addSphericalCoord();
oDetectionsP2P1.addSphericalCoord();
elevations=arrayfun(@(d,D) min(d.elevation,D.elevation),oDetectionsP1P2,oDetectionsP2P1,'unif',0);

trEndElevation=arrayfun(@(t) elevations{t.f(end)}(t.tracksFeatIndxCG(end)),kinTracksISO);
interpolarKin=kinTracksISO(trEndElevation>pi/8);
if(isempty(interpolarKin))
    return;
end
lpid=lpid+1;
if(~isempty(p.package)&&(length(p.package.processes_)>=lpid)&&(~isempty(p.package.getProcess(lpid))))
    processProjMetPlate=p.package.getProcess(lpid);
else
    disp('Produce hard-code ROI');
    tic;
    processProjMetPlate=ExternalProcess(MD,'project1D');
    project1D(  MD,'FoF',refs(1,2), 'processFrame',1:5,'dynPoligonREF',[refs(1,2).applyBase(interpolarKin,'')], ...
        'name','metaphasePlate','channelRender','grayRed', 'fringeWidth',0, ...
        'processSingleProj',processProjMetPlate, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);
    toc;
end
MD.getPackage(packPIDTMP).setProcess(lpid,processProjMetPlate);

lpid=lpid+1;
if(~isempty(p.package)&&(length(p.package.processes_)>=lpid)&&(~isempty(p.package.getProcess(lpid))))
    processDynROIProj=p.package.getProcess(lpid);
    processRenderer=p.package.getProcess(lpid+1);
else
    disp('Produce hard-code ROI');
    tic;
    %processDynROIProj=ComputeMIPProcess(MD);
    processDynROIProj=projectDynROIProcess(MD);
    processRenderer=ExternalProcess(MD,'processRenderer');
    dynROI=copy(refs(1,2).applyBase(interpolarKin,''));
    offset=50;
    dynROIExtended=[dynROI.copy().addOffset(offset,offset,0); dynROI.copy().addOffset(-offset,-offset,0)];
    projectDynROI(  MD,'FoF',refs(1,2),'dynPoligonREF',dynROIExtended, ...
        'name','metaphasePlateDynROI','channelRender','grayRed','rawTIFF',true, 'fringeWidth',0, ...
        'processSingleProj',processDynROIProj, 'processRenderer',processRenderer,  'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);
    toc;
    MD.addProcess(processDynROIProj);
end
MD.getPackage(packPIDTMP).setProcess(lpid,processDynROIProj);
lpid=lpid+1;
MD.getPackage(packPIDTMP).setProcess(lpid,processRenderer);



disp('Overlay hard-code ROI selected detection on lab ref');
tic;
lpid=lpid+1;
processSelectedDetectionOverlay=ExternalProcess(MD,'overlayProjDetectionMovie');
myColormap=255*jet(256);
overlayProjDetectionMovie(processProjMetPlate,'detections',  refs(1,2).applyBase(oDetections,''), ... 
          'colorIndx',arrayfun(@(d) ceil(254*mat2gray(d.xCoord(:,1)))+1, oDetections,'unif',0), ...
          'colormap',myColormap,'process',processSelectedDetectionOverlay,'name','kin')
MD.getPackage(packPIDTMP).setProcess(lpid,processSelectedDetectionOverlay);
toc;

disp('Overlay hard-code ROI tracks');
tic;
lpid=lpid+1;
processOverlay=ExternalProcess(MD,'overlayProjTracksMovie');
overlayProjTracksMovie(processProjMetPlate,'tracks', refs(1,2).applyBase([interpolarKin],'') , ...
                        'name','tracks','process',processOverlay);
% MD.getPackage(packPIDTMP).setProcess(lpid,processOverlay);
toc;



MD.setPackage(packPID,MD.getPackage(packPIDTMP));

