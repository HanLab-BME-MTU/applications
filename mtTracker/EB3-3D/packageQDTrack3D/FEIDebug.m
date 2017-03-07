%%% Description
% This script facilitates debugging from the classic data sent by FEI


%%
[gapCloseParam,costMatrices,kalmanFunctions,...
    probDim,verbose]=QDTrackerParam();
% watch_KF_iter=0;

saveResults.dir =  '.'; %directory where to save input and output
saveResults.filename = 'trackResults.mat'; %name of file where input and output are saved


load('gapCloseParam.mat');
load('detection.mat');
load('costMatricesLink.mat');
load('costMatricesGap.mat');

%% %% USER INPUT 
%costMatricesGap.parameters.brownStdMult=repmat(costMatricesGap.parameters.brownStdMult,gapCloseParam.timeWindow,1);
for i=1:length(detection) detection(i).yCoord(detection(i).yCoord(:,1)==0,1)=1;  end;
for i=1:length(detection) detection(i).xCoord(detection(i).xCoord(:,1)==0,1)=1;  end;

costMatrices(1)=costMatricesLink;
costMatrices(2)=costMatricesGap;


[tracksFinal,kalmanInfoLink,errFlag] = ...
    trackCloseGapsKalmanSparse(detection, ...
    costMatrices,gapCloseParam,kalmanFunctions,...
    probDim,saveResults,verbose);
%%
tracks=TracksHandle(tracksFinal);


if(runAmiraRendering)
    load([outputDir filesep 'tracksHandle.mat']);
    s=[MD.pixelSize_ MD.pixelSize_ MD.pixelSize_ MD.timeInterval_];
    
    amiraWriteTracks([outputDir filesep 'AmiraTracks' filesep 'tracking_' detectionMethod '.am'],tracks,'scales',s);
end
