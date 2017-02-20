function [kinTracksP1,kinTracksP2]=MTDirectionalBias(MD,varargin)
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

printAll=p.printAll;
%%
[kinTracksPlus,EB3TracksPlus,EB3TracksP,kinTracksP]=addSpindleRef(MD);

% inlier index
inliersEB3=(logical(arrayfun(@(eb) eb.inliers(1),EB3TracksPlus)));
inliersKin=logical(arrayfun(@(k) k.inliers(1),kinTracksPlus));


%%
angleCutoff=0.1;

[kinTracksP1,kinTracksP2]=mapMTApparitionToKin(kinTracksP{1}(inliersKin),kinTracksP{2}(inliersKin),EB3TracksP{1}(inliersEB3),EB3TracksP{2}(inliersEB3),angleCutoff);
if(printAll)
    %% For each kinetochore, plot an Amira file with attached mt
    outputDirAmira=[MD.outputDirectory_ filesep 'Kin' filesep 'directionalBias' filesep  'Amira' filesep];
    for kIdx=1:length(kinTracksP2)
        kinTrack=kinTracksP2(kIdx);
        trackSet=[kinTrack; kinTrack.appearingMT];
        trackType=[1; zeros(length(kinTrack.appearingMT),1)];
        amiraWriteTracks([outputDirAmira filesep 'kin_' num2str(kIdx) '.am'],trackSet,'cumulativeOnly',true,'edgeProp',{{'kinEB3',trackType}})
    end
end
       
%% Estimate bundle in KinPole axis

% Randomize pixel domain Kinetochore and create the associted sphercial coordinates.
% randomDist=10; % in pixel
% 
% outputDirTrack=[MD.outputDirectory_ filesep 'Kin' filesep 'track' filesep ];
% kinTrackData=load([outputDirTrack  filesep 'tracksLabRef.mat']);
% kinTracks=kinTrackData.tracksLabRef;
% [randKinTracks]=randomizeKinetochore(kinTracks,randomDist);
% 
% % Translate these changes in the detection structure and associated polar
% % coordiante
% % Load associated data
% outputDirDetect=[MD.outputDirectory_ filesep 'Kin'  filesep 'detection' filesep]
% tmp=load([outputDirDetect 'detectionLabRef.mat']);
% detectionsLabRef=tmp.detectionsLabRef;
% 
% dataIsotropy=[MD.pixelSize_ MD.pixelSize_ MD.pixelSize_];
% EB3poleDetectionMethod=['simplex_scale_003'];
% outputDirPoleDetect=[MD.outputDirectory_ filesep 'EB3' filesep 'poles' filesep EB3poleDetectionMethod filesep];
% poleData=load([outputDirPoleDetect filesep 'poleDetection.mat']);
% poleMovieInfo=poleData.poleMovieInfo;
% 
% [~,kinSphericalCoord,inliers]=tracks2detections(randKinTracks,detectionsLabRef,poleMovieInfo,dataIsotropy)
% 
% % Rebuild the augmented kin
% [randKinTracksPlus]=addSpindleRef(MD,'kinTracks',randKinTracks,'kinSphericalCoord',kinSphericalCoord,'kinInliers',inliers);
% 
% % Estimate bundle outside the Kin-Plan refencial
% randKinTracksInlier=randKinTracksPlus(inliersKin);
% captureDetection(MD,'kinTracks',randKinTracksInlier,'EB3tracks',EB3tracksInliers,'name',['randInlier_angle_' num2str(angleDist)],'distanceCutOff',angleDist);
% randomFiber=detectMTAlignment(MD,'name',['randInlier_angle_' num2str(angleDist)]);
% randFiberCell=[randFiberCell {randomFiber}];


%displayAngleSwipeResults(fiberCell,randFiberCell)
