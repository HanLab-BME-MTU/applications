function astralTracksAnalysis(MD,varargin)
    ip = inputParser;
    ip.CaseSensitive = false;
    ip.KeepUnmatched = true;
    ip.addParameter('package',[]);
    ip.addParameter('debug',[]);
    ip.addParameter('dynROI','');
    ip.addParameter('packPID',150);
    ip.parse(varargin{:});
    p=ip.Results;

% Process type placeholdes
packPID=p.packPID;
% We save the temp processes at packID+1, that way if process get killed
% halfway through, or other issues, the previous logs are not overidden but
% the last computation location are still available.
package=p.package;
% Process type placeholder (they should be defined by tag and then
% retrieved below)
packPIDTMP=packPID+1;
MD.setPackage(packPIDTMP,GenericPackage({ ... 
    ExternalProcess(MD,'project1D'), ...
    PointSourceDetectionProcess3D(MD,[MD.outputDirectory_ filesep 'EB3']),...
    ExternalProcess(MD,'detectPoles'),...
    ExternalProcess(MD,'buildRefsAndROI'),...
    ExternalProcess(MD,'project1D'), ... 
    TrackingProcess(MD,[MD.outputDirectory_ filesep 'EB3']),...   
    ExternalProcess(MD,'project1D'), ...   
    ExternalProcess(MD,'project1D'), ...    
    ExternalProcess(MD,'overlayProjDetectionMovie'), ...
    ExternalProcess(MD,'overlayProjTracksMovie'), ...  
    ExternalProcess(MD,'overlayProjDetectionMovie'), ...
    ExternalProcess(MD,'overlayProjTracksMovie'), ...    
    ExternalProcess(MD,'astralTracksAnalysis') ...
    }));

lpid=0;

lpid=lpid+1;
if(GenericPackage.processExist(package,lpid))
    processProj=p.package.getProcess(lpid);
else
    processProj=ExternalProcess(MD,'project1D');
    project1D(  MD, ...
        'name','fullMIPLabFrame','channelRender','grayRed', ...
        'processSingleProj',processProj, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);
end
MD.getPackage(packPIDTMP).setProcess(lpid,processProj);


lpid=lpid+1;
if(GenericPackage.processExist(package,lpid))
 processDetectEB3=p.package.getProcess(lpid);
else
    processDetectEB3=PointSourceDetectionProcess3D(MD, [MD.outputDirectory_ filesep 'EB3'],UTrackPackage3D.getDefaultDetectionParams(MD,[MD.outputDirectory_ filesep 'EB3']));
    MD.addProcess(processDetectEB3);
    funParams = processDetectEB3.funParams_;
    funParams.showAll=true;
    funParams.alpha=0.05;
    funParams.filterSigma=[1.4396 1.4396;1.2913 1.2913  ];
    funParams.WindowSize={[],[]};
    funParams.algorithmType= {'pointSourceLM'  'pointSourceLM'};
    funParams.ConfRadius={[],[]};       
    processDetectEB3.setPara(funParams);
    paramsIn.ChannelIndex=1;
    paramsIn.isoCoord=true;
    processDetectEB3.run(paramsIn);
end
MD.getPackage(packPIDTMP).setProcess(lpid,processDetectEB3);

lpid=lpid+1;
if(GenericPackage.processExist(package,lpid))
    processDetectPoles=p.package.getProcess(lpid);
else
    processDetectPoles=ExternalProcess(MD,'detectPoles',@(p) detectPoles(p.getOwner(),'process',p,'isoOutput',true));
    processDetectPoles.run();
end
MD.getPackage(packPIDTMP).setProcess(lpid,processDetectPoles);

lpid=lpid+1;
if(GenericPackage.processExist(package,lpid))
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

%% Build static tracks from P1 P2
statP2=copy(P2);
statP2.x(:)=P1.x-P1.x(1)+P2.x(1);
statP2.y(:)=P1.y-P1.y(1)+P2.y(1);
statP2.z(:)=P1.z-P1.z(1)+P2.z(1);

statRefP1=buildRefsFromTracks(P1,statP2,'buildROI',true);

statP1=copy(P1);
statP1.x(:)=P1.x(1);
statP1.y(:)=P1.y(1);
statP1.z(:)=P1.z(1);


statP2=copy(P2);
statP2.x(:)=P2.x(1);
statP2.y(:)=P2.y(1);
statP2.z(:)=P2.z(1);

statRef=buildRefsFromTracks(statP1,statP2,'buildROI',true);


%% 
% inlier tracks
refs=load(processBuildRef.outFilePaths_{1}); refs=refs.refs;
ROIs=load(processBuildRef.outFilePaths_{2}); ROIs=ROIs.ROI;

%%
lpid=lpid+1;
if(GenericPackage.processExist(package,lpid))
    processProjSpindleRef=p.package.getProcess(lpid);
else
    processProjSpindleRef=ExternalProcess(MD,'project1D');
    project1D(  MD,[P1,P2],'FoF',refs(1,2), ...
        'name','CroppedSpindleRef','channelRender','grayRed', 'insetFringeWidth',80, ...
        'processSingleProj',processProjSpindleRef, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);
end
MD.getPackage(packPIDTMP).setProcess(lpid,processProjSpindleRef);

% spindle frame of reference for detection    cellfun(@(p) disp(p.name_),MD.getPackage(150).processes_)
 
tmp=load(processDetectEB3.outFilePaths_{1}); detection=tmp.movieInfo;
oDetections=Detections(detection);
oDetectionsP1P2=refs(1,2).applyBase(oDetections,'');
oDetectionsP2P1=refs(2,1).applyBase(oDetections,'');
oDetectionsP1P2.addSphericalCoord();
oDetectionsP2P1.addSphericalCoord();
elevations=arrayfun(@(d,D) min(d.elevation,D.elevation),oDetectionsP1P2,oDetectionsP2P1,'unif',0);

%% selecting tracks from a single pole
if(~isempty(p.dynROI))
    switch p.dynROI
    case 'hardcoded'
    case 'astral'
        elevationCutoffComp=[-pi/2 0];
        elevationCutoffROI=[-pi/2 -pi/10];
        xCutoffComp=[-20 20];
        xCutoffROI=[-12 12];
        rhoCutoffComp=[0 100];
        rhoROI=[8 100];
        selectedDetection=copy(oDetections);
        for i=1:length(selectedDetection)
            d=oDetectionsP1P2(i);
            selectedDetectionIdx=(d.elevation<elevationCutoffComp(2))&(d.elevation>elevationCutoffComp(1));
            selectedDetectionIdx=selectedDetectionIdx&(d.rho<rhoCutoffComp(2))&(d.rho>rhoCutoffComp(1));
            selectedDetectionIdx=selectedDetectionIdx&((d.xCoord(:,1))<xCutoffComp(2))&((d.xCoord(:,1))>xCutoffComp(1));
            selectedDetection(i).selectIdx(selectedDetectionIdx);
        end
    case 'interpolar'
        elevationCutoffComp=[pi/8 pi/2+1];
        elevationCutoffROI=[pi/8 pi/2+1];
        xCutoffComp=[-10 10];
        xCutoffROI=[-6 6];
%         rhoCutoffComp=[0 100];
        rhoROI=[0 500];
        selectedDetection=copy(oDetections);
        for i=1:length(selectedDetection)
            d=oDetectionsP1P2(i);
            selectedDetectionIdx=(elevations{i}<elevationCutoffComp(2))&(elevations{i}>elevationCutoffComp(1));
            selectedDetectionIdx=selectedDetectionIdx&((d.xCoord(:,1))<xCutoffComp(2))&((d.xCoord(:,1))>xCutoffComp(1));
            selectedDetection(i).selectIdx(selectedDetectionIdx);
        end
    otherwise
        error('Unknown ROI spec');
    end


    detCopy=copy(selectedDetection);
    detCopy.addOffset(100,100,100);
    movieInfo=detCopy.getStruct();

    % save data in a process to feed to u-track
    outputFilePath=fullfile(fileparts(processDetectEB3.outFilePaths_{1}),p.dynROI); mkdirRobust(outputFilePath);
    outputFilePath=fullfile(outputFilePath,'channel_1.mat'); 

    processDetectEB3=PointSourceDetectionProcess3D(MD, [MD.outputDirectory_ filesep 'EB3'],UTrackPackage3D.getDefaultDetectionParams(MD,[MD.outputDirectory_ filesep 'EB3' p.dynROI]));
    save(outputFilePath,'movieInfo');
    processDetectEB3.setOutFilePaths({outputFilePath});
    MD.addProcess(processDetectEB3);    

end;

lpid=lpid+1;    

if(GenericPackage.processExist(package,lpid))
    processTrackEB3=p.package.getProcess(lpid);
else
    processTrackEB3=TrackingProcess(MD, [MD.outputDirectory_ filesep 'EB3'],UTrackPackage3D.getDefaultTrackingParams(MD, [MD.outputDirectory_ filesep 'EB3' p.dynROI]));
    MD.addProcess(processTrackEB3);    
    funParams = processTrackEB3.funParams_;
    [costMatrices,gapCloseParam,kalmanFunctions,probDim]=astralPlusTipCometTracker3DParam(MD);
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

tmp=load(processTrackEB3.outFilePaths_{1}); EB3TracksISO=TracksHandle(tmp.tracksFinal);
EB3TracksISO.addOffset(-100,-100,-100);
EB3TracksInliers=mapTracksTo1DManifold(ROIs{1,2},EB3TracksISO,0,'position','start','distType','vertexDistOtsu');
% detect tracks that finish in the astral region
trEndElevation=arrayfun(@(t) elevations{t.f(end)}(t.tracksFeatIndxCG(end)),EB3TracksInliers);
astralTracks=EB3TracksInliers(trEndElevation<pi/8);


if(~isempty(p.dynROI))
    selectedDetectionP1P2=refs(1,2).applyBase(selectedDetection,'');
    selectedDetectionP1P2.addSphericalCoord();
    trEndElevation=arrayfun(@(t) selectedDetectionP1P2(t.f(end)).elevation(t.tracksFeatIndxCG(end)),EB3TracksISO);
    trEndRho=arrayfun(@(t) selectedDetectionP1P2(t.f(end)).rho(t.tracksFeatIndxCG(end)),EB3TracksISO);
    trEndX=arrayfun(@(t) selectedDetectionP1P2(t.f(end)).xCoord(t.tracksFeatIndxCG(end),1),EB3TracksISO);
    EB3TracksISODisplayOnly=EB3TracksISO((trEndElevation<elevationCutoffROI(2))    ...
                                        &(trEndElevation>elevationCutoffROI(1))    ...
                                        &(trEndRho<rhoROI(2))                   ...
                                        &(trEndRho>rhoROI(1)                    ...
                                        &(trEndX<xCutoffROI(2))                 ...
                                        &(trEndX>xCutoffROI(1))));
    EB3TracksISO=EB3TracksISO((trEndRho>rhoROI(1)));
    astralTracks=EB3TracksISO;
end


lpid=lpid+1;
if(~isempty(p.dynROI))
    if(GenericPackage.processExist(package,lpid))
        processProjStatRef=p.package.getProcess(lpid);
    else
        disp('Produce hard-code ROI stat lab ref');
        tic;
        processProjStatRef=ExternalProcess(MD,'project1D');
        project1D(  MD,'FoF',statRef, 'dynPoligonREF',[statRef.applyBase([P1 ; EB3TracksISODisplayOnly],''); ], ...
            'name',['dynROI-' p.dynROI '-labRef'],'channelRender','grayRed', 'fringeWidth',1, ...
            'processSingleProj',processProjStatRef, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);
    %      project1D(  MD,'dynPoligonREF',[P1 ; EB3TracksISODisplayOnly], ...
    toc;
end
MD.getPackage(packPIDTMP).setProcess(lpid,processProjStatRef);
end

% lpid=lpid+1;
% if(~isempty(p.package)&&(length(p.package.processes_)>=lpid)&&(~isempty(p.package.getProcess(lpid))))
%     processProjStatRefP1=p.package.getProcess(lpid);
% else
%     disp('Produce hard-code ROI stat P1 ref');
%     tic;
%     processProjStatRefP1=ExternalProcess(MD,'project1D');
%     project1D(  MD,'FoF',statRefP1, 'dynPoligonREF',[statRefP1.applyBase([P1 ; EB3TracksISODisplayOnly],''); ], ...
%         'name','dynROI-hardcoded-P1ref','channelRender','grayRed', 'fringeWidth',1, ...
%         'processSingleProj',processProjStatRefP1, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);
%     toc;
% end
% MD.getPackage(packPIDTMP).setProcess(lpid,processProjStatRefP1);
% 
lpid=lpid+1;
if(GenericPackage.processExist(package,lpid))
    processProjdynROI=p.package.getProcess(lpid);
else
    disp('Produce hard-code ROI');
    tic;
    processProjdynROI=ExternalProcess(MD,'project1D');
    project1D(  MD,'FoF',refs(1,2), 'dynPoligonREF',[refs(1,2).applyBase([P1 ; EB3TracksISODisplayOnly],''); ], ...
        'name',['dynROI-' p.dynROI '-refs'],'channelRender','grayRed', 'fringeWidth',1,                              ...
        'processSingleProj',processProjdynROI, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);
    toc;
end
MD.getPackage(packPIDTMP).setProcess(lpid,processProjdynROI);

% disp('Overlay hard-code ROI all detection');
% tic;
% lpid=lpid+1;
% processAllDetectOverlay=ExternalProcess(MD,'overlayProjDetectionMovie');
% myColormap=255*jet(256);
% overlayProjDetectionMovie(processProjdynROI,'detections', oDetectionsP1P2 , ... 
%             'colorIndx',arrayfun(@(d,s) (1+255*double(ismember(d.zCoord(:,1),s.zCoord(:,1)))) , oDetections,selectedDetection,'unif',0), ...
%             'colormap',myColormap,'name','allDetection','process',processAllDetectOverlay)
% MD.getPackage(packPIDTMP).setProcess(lpid,processAllDetectOverlay);
% toc;

% 
% disp('Overlay hard-code ROI detection on full cell');
% tic;
% lpid=lpid+1;
% processOverlay=ExternalProcess(MD,'overlayProjDetectionMovie');
% myColormap=255*jet(256);
% overlayProjDetectionMovie(processProj,'detections', selectedDetection, ... 
%             'colorIndx',arrayfun(@(d) ceil(254*mat2gray(d.xCoord(:,1)))+1, selectedDetection,'unif',0), ...
%             'colormap',myColormap,'name','selectedDetection')
% MD.getPackage(packPIDTMP).setProcess(lpid,processOverlay);
% toc;

disp('Overlay hard-code ROI selected detection on lab ref');
tic;
lpid=lpid+1;
processSelectedDetectionOverlaySRLab=ExternalProcess(MD,'overlayProjDetectionMovie');
myColormap=255*jet(256);
overlayProjDetectionMovie(processProjStatRef,'detections',  statRef.applyBase(selectedDetection,''), ... 
            'colorIndx',arrayfun(@(d) ceil(254*mat2gray(d.xCoord(:,1)))+1, selectedDetection,'unif',0), ...
            'process',processSelectedDetectionOverlaySRLab, ...
            'colormap',myColormap,'name','selectedDetection')
MD.getPackage(packPIDTMP).setProcess(lpid,processSelectedDetectionOverlaySRLab);
toc;

disp('Overlay hard-code ROI tracks on lab ref');
tic;
lpid=lpid+1;
processOverlaySRLab=ExternalProcess(MD,'overlayProjTracksMovie');
overlayProjTracksMovie(processSelectedDetectionOverlaySRLab,'tracks', [statRef.applyBase([P1; EB3TracksISO],'')] , ...
    'colorIndx',[255; 1*ones(size(EB3TracksISO))],'colormap',myColormap,'name','tracks','process',processOverlaySRLab);
MD.getPackage(packPIDTMP).setProcess(lpid,processOverlaySRLab);
toc;

% disp('Overlay hard-code ROI selected detection on P1 ref');
% tic;
% lpid=lpid+1;
% processSelectedDetectionOverlaySRP1=ExternalProcess(MD,'overlayProjDetectionMovie');
% myColormap=255*jet(256);
% overlayProjDetectionMovie(processProjStatRefP1,'detections',  statRefP1.applyBase(selectedDetection,''), ... 
%             'colorIndx',arrayfun(@(d) ceil(254*mat2gray(d.xCoord(:,1)))+1, selectedDetection,'unif',0), ...
%             'process',processSelectedDetectionOverlaySRP1, ...
%             'colormap',myColormap,'name','selectedDetection')
% MD.getPackage(packPIDTMP).setProcess(lpid,processSelectedDetectionOverlaySRP1);
% toc;
% 
% disp('Overlay hard-code ROI tracks on P1 ref');
% tic;
% lpid=lpid+1;
% processOverlaySRP1=ExternalProcess(MD,'overlayProjTracksMovie');
% overlayProjTracksMovie(processSelectedDetectionOverlaySRP1,'tracks', [statRefP1.applyBase([P1; EB3TracksISO],'')] , ...
%     'colorIndx',[255; 1*ones(size(EB3TracksISO))],'colormap',myColormap,'name','tracks','process',processOverlaySRP1);
% MD.getPackage(packPIDTMP).setProcess(lpid,processOverlaySRP1);
% toc;


disp('Overlay hard-code ROI selected detection');
tic;
lpid=lpid+1;
processSelectedDetectionOverlay=ExternalProcess(MD,'overlayProjDetectionMovie');
myColormap=255*jet(256);
overlayProjDetectionMovie(processProjdynROI,'detections',  refs(1,2).applyBase(selectedDetection,''), ... 
            'colorIndx',arrayfun(@(d) ceil(254*mat2gray(d.xCoord(:,1)))+1, selectedDetection,'unif',0), ...
            'process',processSelectedDetectionOverlay, ...
            'colormap',myColormap,'name','selectedDetection')
MD.getPackage(packPIDTMP).setProcess(lpid,processSelectedDetectionOverlay);
toc;

disp('Overlay hard-code ROI tracks');
tic;
lpid=lpid+1;
processOverlay=ExternalProcess(MD,'overlayProjTracksMovie');
overlayProjTracksMovie(processSelectedDetectionOverlay,'tracks', [refs(1,2).applyBase([P1; EB3TracksISO],'')] , ...
    'colorIndx',[255; 1*ones(size(EB3TracksISO))],'colormap',myColormap,'name','tracks','process',processOverlay);
MD.getPackage(packPIDTMP).setProcess(lpid,processOverlay);
toc;

disp('Build MIP')
tic;
printProcMIPArray({ processOverlaySRLab,processOverlay,    ...
                    processOverlaySRLab,processOverlay ...
                    },    ... 
                   [fileparts(fileparts(processOverlaySRLab.outFilePaths_{4})) filesep 'MIPArray' filesep], ...
                    'MIPIndex',[2,2,3,3],'MIPSize',400,'maxHeight',1500,'maxWidth',1000);
toc;

% disp('Overlay full MIP');
% tic;
% lpid=lpid+1;
% processOverlay=ExternalProcess(MD,'overlayProjTracksMovie');
% overlayProjTracksMovie(processProj,'tracks',[P1 ;P2 ;astralTracks;], ... 
%     'colorIndx',[1 ;1; 2*ones(size(astralTracks))],'colormap',myColormap,'name','astralEB3','process',processOverlay);
% MD.getPackage(packPIDTMP).setProcess(lpid,processOverlay);
% toc;

clear tmp;

% TODO Length
st=MD.timeInterval_;
sx=0.1;sz=0.1;
distance=arrayfun(@(t) sum(nansum(( [sx*t.x(1:end-1);sx*t.y(1:end-1);sz*t.z(1:end-1)]- ...
  [sx*t.x(2:end);  sx*t.y(2:end);  sz*t.z(2:end)  ]).^2).^0.5) , astralTracks);

speedStd=arrayfun(@(t) nanstd(nansum(( [sx*t.x(1:end-1);sx*t.y(1:end-1);sz*t.z(1:end-1)]- ...
  [sx*t.x(2:end);  sx*t.y(2:end);  sz*t.z(2:end)  ]).^2).^0.5/st) , astralTracks);

lpid=lpid+1;
processAstralTracks=ExternalProcess(MD,'astralTracksAnalysis');
mkdirRobust([MD.outputDirectory_ filesep 'astralTracks' filesep  p.dynROI]);
save([MD.outputDirectory_ filesep 'astralTracks' filesep p.dynROI filesep 'trackLengthStat.mat'],'distance','speedStd')
processAstralTracks.setOutFilePaths({[MD.outputDirectory_ filesep 'astralTracks' filesep p.dynROI filesep 'trackLengthStat.mat']});
MD.getPackage(packPIDTMP).setProcess(lpid,processAstralTracks);

MD.setPackage(packPID,MD.getPackage(packPIDTMP));

