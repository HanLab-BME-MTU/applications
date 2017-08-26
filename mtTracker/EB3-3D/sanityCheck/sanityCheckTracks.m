function sanityCheckTracks(MD,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('package',[]);
ip.parse(varargin{:});
p=ip.Results;

% Process type placeholdes
packPID=100;
% We save the temp processes at packID+1, that way if process get killed
% halfway through, or other issues, the previous logs are not overidden but
% the last computation location are still available.


% Process type placeholdes
packPIDTMP=packPID+1;
MD.setPackage(packPIDTMP,GenericPackage({ ... 
    PointSourceDetectionProcess3D(MD,[MD.outputDirectory_ filesep 'EB3']),...
    TrackingProcess(MD,[MD.outputDirectory_ filesep 'EB3']),...
    ExternalProcess(MD,'detectPoles'),...
    ExternalProcess(MD,'buildRefsAndROI'),...
    ExternalProcess(MD,'sanityCheckTracks')
  }));

lpid=1;
if(~isempty(p.package)&&(~isempty(p.package.getProcess(lpid))))
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
if(~isempty(p.package)&&(~isempty(p.package.getProcess(lpid))))
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
if(~isempty(p.package)&&(~isempty(p.package.getProcess(lpid))))
    processDetectPoles=p.package.getProcess(lpid);
else
    processDetectPoles=ExternalProcess(MD,'detectPoles',@(p) detectPoles(p.getOwner(),'process',p,'isoOutput',true));
    processDetectPoles.run();
end

MD.getPackage(packPIDTMP).setProcess(lpid,processDetectPoles);

lpid=lpid+1;
if(~isempty(p.package)&&(~isempty(p.package.getProcess(lpid))))
    processBuildRef=p.package.getProcess(lpid);
else
    processBuildRef=ExternalProcess(MD,'buildRefsAndROI',@(p) buildRefsFromTracks(processDetectPoles,processDetectPoles,'buildROI',true,'process',p));
    processBuildRef.run();
end

MD.getPackage(packPIDTMP).setProcess(lpid,processBuildRef);



%% Load tracks, convert and select inliers
ROIs=load(processBuildRef.outFilePaths_{2}); ROIs=ROIs.ROI;
tmp=load(processTrackEB3.outFilePaths_{1}); EB3TracksISO=TracksHandle(tmp.tracksFinal);
EB3TracksInliers=mapTracksTo1DManifold(ROIs{1,2},EB3TracksISO,0,'position','start','distType','vertexDistOtsu');

clear tmp;

% speed
st=MD.timeInterval_;
sx=0.1;sz=0.1;
avgspeed=arrayfun(@(t) mean(nansum(( [sx*t.x(1:end-1);sx*t.y(1:end-1);sz*t.z(1:end-1)]- ...
                                  [sx*t.x(2:end);  sx*t.y(2:end);  sz*t.z(2:end)  ]).^2).^0.5/st) ,EB3TracksInliers);
medspeed=arrayfun(@(t) median(nansum(( [sx*t.x(1:end-1);sx*t.y(1:end-1);sz*t.z(1:end-1)]- ...
                                  [sx*t.x(2:end);  sx*t.y(2:end);  sz*t.z(2:end)  ]).^2).^0.5/st) ,EB3TracksInliers);


lpid=lpid+1;
processEnrichVsElev=ExternalProcess(MD,'sanityCheckTracks');
mkdirRobust([MD.outputDirectory_ filesep 'sanity' filesep ]);
save([MD.outputDirectory_ filesep 'sanity' filesep 'growthRate-density.mat'],'avgspeed','medspeed')
processEnrichVsElev.setOutFilePaths({[MD.outputDirectory_ filesep 'sanity' filesep 'growthRate-density.mat']});
MD.getPackage(packPIDTMP).setProcess(lpid,processEnrichVsElev);

MD.setPackage(packPID,MD.getPackage(packPIDTMP));

