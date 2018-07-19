function [processOverlayCells]=testImaris(MD,varargin)
    ip = inputParser;
    ip.CaseSensitive = false;
    ip.KeepUnmatched = true;
    ip.addParameter('package',[]);
    ip.addParameter('debug',false);
    ip.addParameter('name','');
    ip.addParameter('iceConnector',[])
    ip.addParameter('dynROIView',[]);
    ip.addParameter('forceRunIdx',[]);
    %ip.addParameter('dynROIData',[]);   % Side loading projection
    ip.parse(varargin{:});
    p=ip.Results;


templatePackage=GenericPackage({ ... 
    PointSourceDetectionProcess3D(MD,[MD.outputDirectory_ ]),...
    TrackingProcess(MD,[MD.outputDirectory_ ])  
    },[],'name_',['testImaris' p.name '_backup']);


package=p.package;

MD.addPackage(templatePackage);

iceConnector=p.iceConnector;
if(isempty(iceConnector))
    try
        iceConn = movieViewerImaris(MD,'UseImarisFileReader',true);
    catch
        iceConn = movieViewerImaris(MD);
    end
end


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
    funParams.algorithmType= {'pointSourceFit'}
    funParams.ConfRadius={[]};       
    processDetect.setPara(funParams);
    paramsIn.ChannelIndex=2;
    paramsIn.isoCoord=true;
    processDetect.run(paramsIn);    
end
templatePackage.setProcess(lpid,processDetect);

disp('loading detections')
tmp=load(processDetect.outFilePaths_{2}); detection=tmp.movieInfo;
oDetections=Detections(detection);


disp('Computing trackability');
% [maxSpeed,~,densities]=estimateTrackability(oDetections,1/MD.timeInterval_,'debugMode','');


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
    paramsIn.ChannelIndex=2;
    paramsIn.DetProcessIndex=processDetect.getIndex();
    processTrack.run(paramsIn);
end
templatePackage.setProcess(lpid,processTrack);

%% Setup track back in the normal referential
disp('loading tracks')
tmp=load(processTrack.outFilePaths_{2}); 
tracks=TracksHandle(tmp.tracksFinal); 
imarisShowTracks(tracks,iceConnector);

completePackage=GenericPackage(templatePackage.processes_,[],'name_',['testImaris' p.name]);
MD.addPackage(completePackage);
