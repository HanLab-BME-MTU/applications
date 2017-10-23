function tracks3DAP2(MD,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('package',[]);
ip.addParameter('createROIs',false);
ip.addParameter('packID',100);
ip.addParameter('name',[]);
ip.addParameter('forceRun',[]);
ip.parse(varargin{:});
p=ip.Results;

forceRun=p.forceRun;

% Process type placeholdes
packPID=p.packID;
% We save the temp processes at packID+1, that way if process get killed
% halfway through, or other issues, the previous logs are not overidden but
% the last computation location are still available.

nFrames=min(300,MD.nFrames_);

% Process type placeholder
packPIDTMP=packPID+1;

if(p.createROIs)
    processROIs=ExternalProcess(MD,'createROIs');
    [~,mask]=createROI(MD);
    MD=crop3D(MD,mask);
end

MD.setPackage(packPIDTMP,GenericPackage({ ... 
    ComputeMIPProcess(MD,[MD.outputDirectory_ filesep 'AP2']),...
    PointSourceDetectionProcess3D(MD,[MD.outputDirectory_ filesep 'AP2']),...
    TrackingProcess(MD,[MD.outputDirectory_ filesep 'AP2']),...
    ExternalProcess(MD,'rawProj'), ...
    ExternalProcess(MD,'trackOverlay')
  }));

lpid=0;

lpid=lpid+1;
if(~isempty(p.package)&&(~isempty(p.package.getProcess(lpid))))
    processMIP=p.package.getProcess(lpid);
    if(~isempty(forceRun)&&forceRun(lpid))
        processMIP.run();
    end
else
    processMIP=ComputeMIPProcess(MD);
    MD.addProcess(processMIP);
    processMIP.run();
end
MD.getPackage(packPIDTMP).setProcess(lpid,processMIP);

lpid=lpid+1;
if(~isempty(p.package)&&(~isempty(p.package.getProcess(lpid))))
    processDetectAP2=p.package.getProcess(lpid);
    if(~isempty(forceRun)&&forceRun(lpid))
        processDetectAP2.run();
    end
else
    processDetectAP2=PointSourceDetectionProcess3D(MD,[MD.outputDirectory_ filesep 'AP2'],UTrackPackage3D.getDefaultDetectionParams(MD, [MD.outputDirectory_ filesep 'AP2']));
    MD.addProcess(processDetectAP2);
    funParams = processDetectAP2.funParams_;
    funParams.showAll=true;
    funParams.alpha=0.05;
    funParams.frameRange=[1 nFrames];    
    funParams.fitMixtures = false;
    funParams.MaxMixtures = 5;   
    funParams.filterSigma=[1.6; 1.5  ];
    funParams.WindowSize={[]};
    funParams.algorithmType= {'pointSourceFit'};
    funParams.ConfRadius={[]};       
    processDetectAP2.setPara(funParams);
    paramsIn.ChannelIndex=1;
    paramsIn.isoCoord=true;
    processDetectAP2.run(paramsIn);    
end
MD.getPackage(packPIDTMP).setProcess(lpid,processDetectAP2);

lpid=lpid+1;
if(~isempty(p.package)&&(~isempty(p.package.getProcess(lpid))))
    processTrack=p.package.getProcess(lpid);
    if(~isempty(forceRun)&&forceRun(lpid))
        processTrack.run();
    end
else
    processTrack=TrackingProcess(MD, [MD.outputDirectory_ filesep 'AP2'],UTrackPackage3D.getDefaultTrackingParams(MD,[MD.outputDirectory_ filesep 'AP2']));
    MD.addProcess(processTrack);    
    funParams = processTrack.funParams_;
    newFunParams=AP2TrackingParam();
    F=fields(newFunParams); for fIdx=1:length(F) funParams.(F{fIdx})=newFunParams.(F{fIdx}); end;
    processTrack.setPara(funParams);
    paramsIn.ChannelIndex=1;
    paramsIn.DetProcessIndex=processDetectAP2.getIndex();
    processTrack.run(paramsIn);
end
MD.getPackage(packPIDTMP).setProcess(lpid,processTrack);


tic
lpid=lpid+1;
if(~isempty(p.package)&&(~isempty(p.package.getProcess(lpid))))
    processProj=p.package.getProcess(lpid);
    if(~isempty(forceRun)&&forceRun(lpid))
        processProj.run();
    end
else
    disp('MIP')
    processProj=ExternalProcess(MD,'rawProj');
    project1D(  MD, ...
        'name','fullMIPNoManifold','channelRender','grayRed', ...
        'processSingleProj',processProj,'processFrame',1:nFrames, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);
end
MD.getPackage(packPIDTMP).setProcess(lpid,processProj);


lpid=lpid+1;
tic;
disp('Loading and Converting tracks')
tracksFinal=processTrack.loadChannelOutput(1);
tracks=TracksHandle(tracksFinal);
toc;

disp('Overlay')
tic;
NColor=256;
cMap=255*jet(NColor);
cIdx=arrayfun(@(t) t.z(1),tracks);
cIdx=ceil((NColor-1)*mat2gray(cIdx))+1;
processOverlay=ExternalProcess(MD,['track-depth' p.name]);
overlayProjTracksMovie(processProj,'tracks',tracks, ... 
            'colorIndx',cIdx,'colormap',cMap,'name',processOverlay.name_,'process',processOverlay);
toc;
MD.getPackage(packPIDTMP).setProcess(lpid,processOverlay);

% Save tmp package at final memory location.
MD.setPackage(packPID,MD.getPackage(packPIDTMP));
