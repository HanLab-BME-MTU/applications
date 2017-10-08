function tracks3DAP2(MD,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('package',[]);
ip.addParameter('createROIs',false);
ip.parse(varargin{:});
p=ip.Results;

% Process type placeholdes
packPID=100;
% We save the temp processes at packID+1, that way if process get killed
% halfway through, or other issues, the previous logs are not overidden but
% the last computation location are still available.

nFrames=min(200,MD.nFrames_);

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
else
    processMIP=ComputeMIPProcess(MD);
    MD.addProcess(processMIP);
    processMIP.run();    
end
MD.getPackage(packPIDTMP).setProcess(lpid,processMIP);

lpid=lpid+1;
if(~isempty(p.package)&&(~isempty(p.package.getProcess(lpid))))
   processDetectAP2=p.package.getProcess(lpid);
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
else
    processTrack=TrackingProcess(MD, [MD.outputDirectory_ filesep 'AP2'],UTrackPackage3D.getDefaultTrackingParams(MD,[MD.outputDirectory_ filesep 'AP2']));
    MD.addProcess(processTrack);    
    funParams = processTrack.funParams_;
    [gapCloseParam,costMatrices,kalmanFunctions,probDim,verbose]=kinTrackingParam();
    funParams.gapCloseParam=gapCloseParam;
    funParams.costMatrices=costMatrices;
    funParams.kalmanFunctions=kalmanFunctions;
    funParams.probDim=probDim;
    processTrack.setPara(funParams);
    paramsIn.ChannelIndex=1;
    paramsIn.DetProcessIndex=processDetectAP2.getIndex();
    processTrack.run(paramsIn);
end
MD.getPackage(packPIDTMP).setProcess(lpid,processTrack);


tic
disp('MIP')
lpid=lpid+1;
processProj=ExternalProcess(MD,'rawProj');
project1D(  MD, ...
            'name','fullMIPNoManifold','channelRender','grayRed', ...
            'processSingleProj',processProj,'processFrame',1:nFrames, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]); 
toc;
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
processOverlay=ExternalProcess(MD,'trackOverlay');
overlayProjTracksMovie(processProj,'tracks',tracks, ... 
            'colorIndx',cIdx,'colormap',cMap,'name','track-depth','process',processOverlay);
toc;
MD.getPackage(packPIDTMP).setProcess(lpid,processOverlay);

% Save tmp package at final memory location.
MD.setPackage(packPID,MD.getPackage(packPIDTMP));
