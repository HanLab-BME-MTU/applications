function MTDirectionalBias(MD,varargin)
% Swipe different angle used to consider capute Detection, output the
% bundle, save them if a process is input. 
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addParameter('kinBundle',[]);
ip.addParameter('printAll',false, @islogical);
ip.addParameter('angleCutoff',0.1, @isnumeric);
ip.addParameter('plotHandle',[]);
ip.addParameter('process',[]);
ip.parse(MD,varargin{:});
p=ip.Results;

printAll=p.printAll;
%%
[kinTracks,EB3Tracks]=addSpindleRef(MD);

% inlier index
inliersEB3=(logical(arrayfun(@(eb) eb.inliers(1),EB3Tracks)));
inliersKin=logical(arrayfun(@(k) k.inliers(1),kinTracks));

outputDirAmira=[MD.outputDirectory_ filesep 'Kin' filesep 'directionalBias' filesep  'Amira' filesep];

% Randomize pixel domain Kinetochore and create the associted sphercial coordinates.
randomDist=10; % in pixel

outputDirTrack=[MD.outputDirectory_ filesep 'Kin' filesep 'track' filesep ];
kinTrackData=load([outputDirTrack  filesep 'tracksLabRef.mat']);
kinTracksOrig=kinTrackData.tracksLabRef;
[randKinTracks]=randomizeKinetochore(kinTracksOrig,randomDist);

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

% Estimate bundle outside the Kin-Plan refencial
randKinTracksInlier=randKinTracksPlus(inliersKin);


%%
angleCutoffs=p.angleCutoff;
numAngle=length(angleCutoffs);
[Hs]=setupFigure(1,(numAngle),numAngle,'AspectRatio',1,'AxesWidth',4);
for aIdx=1:numAngle
    angleCutoff=angleCutoffs(aIdx);
    [kinTracks]=mapMTApparitionToKin(kinTracks(inliersKin),EB3Tracks(inliersEB3),angleCutoff);
    appearingMTBias(kinTracks);
    % For each kinetochore, plot an Amira file with attached mt
    if(printAll)
        %% For each kinetochore, plot an Amira file with attached mt
        for kIdx=1:length(kinTracks)
            trackSet=[kinTracks(kIdx); kinTracks(kIdx).appearingMTP1];
            trackType=[1; zeros(length(kinTracks(kIdx).appearingMTP1),1)];
            amiraWriteTracks([outputDirAmira filesep 'kin_' num2str(kIdx) '.am'],trackSet,'cumulativeOnly',true,'edgeProp',{{'kinEB3',trackType}})
        end
        for kIdx=1:length(kinTracks)
            trackSet=[kinTracks(kIdx); kinTracks(kIdx).appearingMTP1; kinTracks(kIdx).appearingMTP2];
            trackType=[4; kinTracks(kIdx).MTP1Angle; kinTracks(kIdx).MTP2Angle];
            amiraWriteTracks([outputDirAmira filesep 'kin_' num2str(kIdx) '.am'],trackSet,'cumulativeOnly',true,'edgeProp',{{'kinEB3',trackType}})
        end
        %% For each kinetochore, plot an Amira file with attached mt in the spindle
        % ref
        for kIdx=1:length(kinTracks)
            appearingMT=[kinTracks(kIdx).appearingMTP1];
            if (~isempty(appearingMT))
                trackSet=[kinTracks(kIdx).pole1,[appearingMT.pole1]];
            else
                trackSet=kinTracks(kIdx).pole1;
            end
            trackType=[4; kinTracks(kIdx).MTP1Angle];
            amiraWriteTracks([outputDirAmira filesep 'spindleRefP1' filesep 'kin_' num2str(kIdx) '.am'],trackSet,'cumulativeOnly',true,'edgeProp',{{'kinEB3',trackType}})
        end
    end
    
    %% Estimate bias outside KinPole axis
    [randKinTracksInlier]=mapMTApparitionToKin(randKinTracksInlier,EB3Tracks(inliersEB3),angleCutoff);
    appearingMTBias(randKinTracksInlier);
    displayBiasStat({kinTracks,randKinTracksInlier},'plotHandleArray',Hs(aIdx));
end

%displayAngleSwipeResults(fiberCell,randFiberCell)
