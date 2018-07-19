% A simple script to test the use of multiple detection and tracking processes for MPT on dual channel data
% the results should be the same as scriptEB3DetectAndTracks and
% ScriptKinetochoresDetectAndTracks.

% Starting from 

% MATLAB-based testing/performance suite for Utrack3D
% Andrew R. Jamieson 2017
% Test Utrack3D Package
%-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Preconditions
%-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Init MD in a new folder 
%%
MD=MDOrig.addAnalysisFolder('C:\Users\Philippe\project-local\externBetzig\analysis\adavid\smallSample\prometaphase\','C:\Users\Philippe\project-local\externBetzig\analysis\package\trackProcess');
MD.reset();

%% Load MD
MD=MovieData.load('C:\Users\Philippe\project-local\externBetzig\analysis\package\trackProcess\earlyCell1_12\movieData.mat');

%%

% ---------------------------------------------------
% ---------------------------------------------------
% Generate Package
Package_ = UTrackPackage3D(MD);
% Set-up analysis infrastructure via command line interace
MD.addPackage(Package_);
% Specify processes to include in test.
stepNames = Package_.getProcessClassNames;
disp('===============================');
disp('Available Package Process Steps');
disp('===============================');
disp(stepNames');
steps2Test = [1, 2, 3];
assert(length(Package_.processes_) >= length(steps2Test));
assert(length(Package_.processes_) >= max(steps2Test));

disp('Selected Package Process Steps');
for i=steps2Test
    disp(['Step ' num2str(i) ': ' stepNames{i}]);
end

%%% BROWNIAN MOTION: (CHANNEL 2: KINETOCHORE)
   
%% Step 1: Detection 3D
disp('===================================================================');
disp('Running (1st) Utrack 3D Process');
disp('===================================================================');
iPack = 1;
step_ = 1;
if isempty(MD.getPackage(iPack).processes_{step_})
    MD.getPackage(iPack).createDefaultProcess(step_);
end
funParams = MD.getPackage(iPack).processes_{step_}.funParams_
funParams.showAll=true;
funParams.Alpha=0.05;
funParams.filterSigma=[1.6 1.6;1.5 1.5  ];
funParams.WindowSize={[],[]};
funParams.algorithmType= {'pointSourceFit'  'pointSourceFit'}
funParams.ConfRadius={[],[]};       
MD.getPackage(iPack).getProcess(step_).setPara(funParams);
paramsIn.ChannelIndex=2;
MD.getPackage(iPack).processes_{step_}.run(paramsIn);
%% Estimed scales: 1.7087       1.521

%% Compare results with script based detection 
%% Due to the stochastic aspect of the scale, it does not really work out.
detectProcess=MD.getPackage(iPack).processes_{step_}.loadChannelOutput(2);
size(detectProcess(1).xCoord)

%%
outputDirDetect=[MDOrig.outputDirectory_ filesep 'Kin'  filesep 'detection' filesep];
tmp=load([outputDirDetect 'detectionLabRef.mat']);
detectionScript=tmp.detectionsLabRef;
size(detectionScript(1).xCoord)

assert(detectionScript(1).xCoord==detectProcess(1).xCoord)
%%
[detectionsLabRef,~]=detectEB3(MDOrig,'type','pointSourceFit','showAll',false,'channel',2,'scales',[1.6 1.5]);
%Estimed scales: 1.552      1.5818
size(detectionsLabRef(1).xCoord)
assert(all(all(detectProcess(1).xCoord==detectionsLabRef(1).xCoord)))


%% Step 2: Tracking 3D
disp('===================================================================');
disp('Running (2nd) Utrack 3D Process');
disp('===================================================================');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Step 2: Tracking 3D');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
step_ = 2;
if isempty(MD.getPackage(iPack).processes_{step_})
    MD.getPackage(iPack).createDefaultProcess(step_);
end
funParams = MD.getPackage(iPack).processes_{step_}.funParams_
[gapCloseParam,costMatrices,kalmanFunctions,probDim,verbose]=kinTrackingParam();
funParams.gapCloseParam=gapCloseParam;
funParams.costMatrices=costMatrices;
funParams.kalmanFunctions=kalmanFunctions;
funParams.probDim=probDim;

% MD.getPackage(iPack).getProcess(step_).setPara(funParams);
paramsIn.ChannelIndex=2
MD.getPackage(iPack).processes_{step_}.run(paramsIn);

%% Step 3: Post Tracking Processing Motion Analysis 3D
disp('===================================================================');
disp('Running (3rd) Utrack 3D Process');
disp('===================================================================');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Step 3: Post Tracking Processing Motion Analysis 3D');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
step_ = 3;
if isempty(MD.getPackage(iPack).processes_{step_})
    MD.getPackage(iPack).createDefaultProcess(step_);
end
funParams = MD.getPackage(iPack).processes_{step_}.funParams_
% MD.getPackage(iPack).getProcess(step_).setPara(funParams);
paramsIn.ChannelIndex=2;
MD.getPackage(iPack).processes_{step_}.run(paramsIn);

%% GUI: Package and movieViwer with Overlays
disp('===================================================================');
disp('Running (GUI Output Display Generation) Utrack 3D ');
disp('===================================================================');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Checking Basic GUI functionlality');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
% Load GUI (with Results)
h = MD.getPackage(iPack).GUI(MD);


%%% +TIP TRACKING: (CHANNEL 2: KINETOCHORE)

%% Step 1: Detection 3D
disp('===================================================================');
disp('Running (1st) Utrack 3D Process');
disp('===================================================================');
iPack = 2;
step_ = 1;
if isempty(MD.getPackage(iPack).processes_{step_})
    MD.getPackage(iPack).createDefaultProcess(step_);
end
funParams = MD.getPackage(iPack).processes_{step_}.funParams_
funParams.showAll=true;
funParams.alpha=0.05;
funParams.filterSigma=[1.4396 1.4396;1.2913 1.2913  ];
funParams.WindowSize={[],[]};
funParams.algorithmType= {'pointSourceLM'  'pointSourceLM'}
funParams.ConfRadius={[],[]};       
MD.getPackage(iPack).getProcess(step_).setPara(funParams);
paramsIn.ChannelIndex=1;
paramsIn.isoCoord=true;
MD.getPackage(iPack).processes_{step_}.run(paramsIn);

%
outputDirDetect=[MD.outputDirectory_ filesep 'EB3'  filesep 'detection' filesep];
tmp=load([outputDirDetect 'detectionLabRef.mat']);
detectionScript=tmp.detectionsLabRef;
size(detectionScript(20).xCoord)
detectionScript(20).zCoord(1,:)
%
detectProcess=MD.getPackage(2).processes_{1}.loadChannelOutput(1);
size(detectProcess(20).xCoord)
detectProcess(20).zCoord(1,:)

%%
% [detectionsLabRef,~]=detectEB3(MDOrig,'type','pointSourceAutoSigmaLM','showAll',false,'channel',1,'scales',[1.4396 1.2913]);
% size(detectionsLabRef(1).xCoord)
% assert(all(all(detectProcess(1).xCoord==detectionsLabRef(1).xCoord)))

%% Step 2: Tracking 3D
disp('===================================================================');
disp('Running (2nd) Utrack 3D Process');
disp('===================================================================');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Step 2: Tracking 3D');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
step_ = 2;
if isempty(MD.getPackage(iPack).processes_{step_})
    MD.getPackage(iPack).createDefaultProcess(step_);
end
funParams = MD.getPackage(iPack).processes_{step_}.funParams_;
[costMatrices,gapCloseParam,kalmanFunctions,probDim]=plusTipCometTracker3DParam(MD);
funParams.gapCloseParam=gapCloseParam;
funParams.costMatrices=costMatrices;
funParams.kalmanFunctions=kalmanFunctions;
funParams.probDim=probDim;
MD.getPackage(iPack).processes_{step_}.setPara(funParams);

paramsIn.ChannelIndex=1;
paramsIn.DetProcessIndex=MD.getPackage(iPack).getProcess(1).getIndex();
MD.getPackage(iPack).processes_{step_}.run(paramsIn);

% Compare against script
%
outputDirTrack=[MD.outputDirectory_ filesep 'EB3' filesep 'track' filesep ];
plusTipTrackerLegacyRoot=[outputDirTrack filesep 'plustiptrackerio' filesep];
tmp=load([plusTipTrackerLegacyRoot filesep 'track' filesep 'trackResults.mat']);
tracksFinal=tmp.tracksFinal;
size(tracksFinal)

%
EB3TrackOutput=MD.getPackage(iPack).processes_{2}.loadChannelOutput(1);
size(EB3TrackOutput)


