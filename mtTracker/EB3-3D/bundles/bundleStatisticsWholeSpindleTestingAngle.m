function [fiberCell,randFiberCell]=bundleStatisticsWholeSpindleTestingAngle(MD,varargin)
% Swipe different angle used to consider capute Detection, output the
% bundle, save them if a process is input. 
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addParameter('printAll',false, @islogical);
ip.addParameter('process',[]);
ip.parse(MD,varargin{:});
p=ip.Results;

[kinTracksPlus,EB3tracksPlus]=addSpindleRef(MD);
% inlier index
EB3tracksInliers=EB3tracksPlus(logical(arrayfun(@(eb) eb.inliers(1),EB3tracksPlus)));
inliersKin=logical(arrayfun(@(k) k.inliers(1),kinTracksPlus));
kinTracksInliers=kinTracksPlus(inliersKin);

% Randomize pixel domain Kinetochore and create the associted sphercial coordinates.
randomDist=10; % in pixel

outputDirTrack=[MD.outputDirectory_ filesep 'Kin' filesep 'track' filesep ];
kinTrackData=load([outputDirTrack  filesep 'tracksLabRef.mat']);
kinTracks=kinTrackData.tracksLabRef;
[randKinTracks]=randomizeKinetochore(kinTracks,randomDist);

% Translate these changes in the detection structure and associated polar
% coordiante
% Load associated data
outputDirDetect=[MD.outputDirectory_ filesep 'Kin'  filesep 'detection' filesep]
tmp=load([outputDirDetect 'detectionLabRef.mat']);
detectionsLabRef=tmp.detectionsLabRef;

dataIsotropy=[MD.pixelSize_ MD.pixelSize_ MD.pixelSize_];
EB3poleDetectionMethod=['simplex_scale_003'];
outputDirPoleDetect=[MD.outputDirectory_ filesep 'EB3' filesep 'poles' filesep EB3poleDetectionMethod filesep];
poleData=load([outputDirPoleDetect filesep 'poleDetection.mat']);
poleMovieInfo=poleData.poleMovieInfo;

[~,kinSphericalCoord,inliers]=tracks2detections(randKinTracks,detectionsLabRef,poleMovieInfo,dataIsotropy)

% Rebuild the augmented kin
[randKinTracksPlus]=addSpindleRef(MD,'kinTracks',randKinTracks,'kinSphericalCoord',kinSphericalCoord,'kinInliers',inliers);
randKinTracksInlier=randKinTracksPlus(inliersKin);

angleRange=0.02:0.01:0.1;

fiberCell=cell(1,length(angleRange));
randFiberCell=cell(1,length(angleRange));

parfor angleDistIdx=1:length(angleRange)
    
    angleDist=angleRange(angleDistIdx);
    
    %% Estimate bundle in KinPole axis
    captureDetection(MD,'kinTracks',kinTracksInliers,'EB3tracks',EB3tracksInliers,'name',['inlier_angle_' num2str(angleDist)],'distanceCutOff',angleDist);
    fiber=detectMTAlignment(MD,'name',['inlier_angle_' num2str(angleDist)]);
    fiberCell{angleDistIdx}=fiber;
    
    % Estimate bundle outside the Kin-Plan refencial
    captureDetection(MD,'kinTracks',randKinTracksInlier,'EB3tracks',EB3tracksInliers,'name',['randInlier_angle_' num2str(angleDist)],'distanceCutOff',angleDist);
    randomFiber=detectMTAlignment(MD,'name',['randInlier_angle_' num2str(angleDist)]);
    randFiberCell{angleDistIdx}=randomFiber;
end

displayAngleSwipeResults(fiberCell,randFiberCell)
