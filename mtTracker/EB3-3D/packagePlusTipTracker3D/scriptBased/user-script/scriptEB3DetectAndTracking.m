%% %% Description
% This script measure and process EB3s trajectories. The resulting tracks
% of each movies are saved in their respective "analysis" in the original
% u-track format, a more use-friendly structured format and an
% Amira-readable format for visualization. If required, stage shift can be
% registered automatically to correct detection coordinates.

%% %% USER INPUT 

%% LOADING MOVIES INFORMATION (TWO OPTIONS)
%% Loading MovieList file
% MovieList file paths:
% Multiple movieList can be loaded by using a collection of
% movieListFileNames:
MLPath='/project/bioinformatics/Danuser_lab/externBetzig/raw/Yuko/Processed data/anaphase2/analysis/';
movieListFileNames={'movieList.mat'};

% An other example of usage with movieList in different folder
%MLPath='/work/gdanuser/proudot/project/EB3-3D-track/packaging/alpha/plusTipTracker3D-alpha-1/data/A1_HeLa_Cells_EB1/'
%movieListFileNames={'prometaphase/analysis/movieList.mat','prometaphase/analysis/movieList.mat'};

% Build the array of MovieList (automatic)
aMovieListArray=cellfun(@(x) MovieList.loadMatFile([MLPath filesep x]), movieListFileNames,'unif',0);
aMovieListArray=aMovieListArray{:};

%% ALTERNATIVELY, an array of movieList object are already loaded can be used (see loadMoviesManagementFile.m).
%  aMovieListArray=[ inter10nML inter33nML inter100nML];

%% PROCESS MANAGEMENT 
%  - 1 to (re-)run the algorithm
%  - 0 to load previously computed results (if any)
runEB3Detection=1;  % Detect EB3 "comets"   
runPoleDetection=1; % Detec spindle Poles
runRegistration=0;  % Register sudden stage shift during acquisition.
runEB3Tracking=1; runAmiraWrite=1;      
runPostProcessing=1;      
% 
% Frame to be processed ([] for whole movie)
processFrame=[];


%% ALGORITHM PARAMETERS

%% EB3 Detection parameters
% The detection method must be chosen depending on the data at hand
% - 'pointSourceAutoSigmaLM' is parameter-free and suitable for raw 
%   (but deskewed) data
% - 'bandPassWatershed' is more suitable for deconvolve data, requires 
%   input threshold parameter <waterThresh>.
detectionMethod='pointSourceAutoSigmaLM';
scales=[];
%scales=[1.4396 1.2913];
scales=[2 1.5];
alpha=0.0005;
waterThresh=80; % Only for bandPassWatershed

%% Pole detection parameters
poleScale=3;  % expected scale of the poles, should not need to be changed.


%% Tracking parameters
% Segments presenting a lifetime below <minTrackLength> are discarded. 
minTrackLength=2; % In frames.
 
% Lower and upper bound for the radius of the sphere considered to link 
% a detection one frame to the next
minRadius=2; % In pixels.
maxRadius=5; % In pixels.
searchRadiusMult=3;
searchRadiusFirstIteration=20;

% In the context of pause or shrinkage detection, <maxFAngle> is the
% maximum angle between the estimated speed vectors at the end
% of a segment and a possible link between two segments. Also the maximum angle
% between the speed estimated the beginning and the end of a segment.
maxFAngle=10; % In degrees.

% In the context of shrinkage detection: <maxBAngle>  defines the area
% considered  around the segment to detect shrinkage event. The angle only
% defines the distance orthoganal distance between a point and track (as if
% the track was a straight line). See Figure 5B in Applegate and al 2012
% for an illustration.
maxBAngle=10; % In degrees.

% Maximum skrinkage factor (multiply with the maximum growth rate) 
maxShrinkFactor=1.5;

% The fluctuation radius <fluctRad> models unexpected fluctuations during
% shrinkage or pauses. The search volume is defined bye area described by
% the parameter above and dilatted by the fluctuation radius.
fluctRad=1.; % In pixels.

% Break non linear tracks into multiple tracks.
breakNonLinearTracks=false;

%% Post processing parameter
sphericalProjectionRadius=5000;

%% DEBUGGING PARAMETER (Could void your warranty)
trackAfterRegistration=0;
% MDIndex=[1:length()];


%% %% Process each movie
for aPoleScale=poleScale
    for k=1:length(aMovieListArray)
        ML=aMovieListArray(k);
        
        %% loop over each different cell in each condition
        for i=1:length(ML.movieDataFile_)           
            MD=MovieData.loadMatFile(ML.movieDataFile_{i});
            disp(['Processing movie: ' MD.getFullPath]);
            dataAnisotropy=[MD.pixelSize_ MD.pixelSize_ MD.pixelSizeZ_];
            
            %% EB3 Detection
            outputDirDetect=[MD.outputDirectory_ filesep 'EB3' filesep 'detection' filesep];mkdir(outputDirDetect);
            if(runEB3Detection)
                [detectionsLabRef,lab]=detectEB3(MD,'type',detectionMethod,'waterThresh',waterThresh,'showAll',false,'scales',scales,'Alpha',alpha);
                save([outputDirDetect filesep 'detectionLabRef.mat'],'detectionsLabRef');
                amiraWriteMovieInfo([outputDirDetect filesep 'Amira' filesep 'amiraVertexLabRef' filesep  'detectionsLabRef.am'], detectionsLabRef,'scales',dataAnisotropy);
                mkdir([outputDirDetect filesep 'detectionMaskLabRef']);
                for tidx=1:length(detectionsLabRef)
                    stackWrite(lab{tidx},[outputDirDetect filesep 'detectionMaskLabRef' filesep 'detect_T_' num2str(tidx,'%05d') '.tif']);
                end
            else
                if (exist([outputDirDetect filesep 'detectionLabRef.mat'], 'file') == 2)
                    tmp=load([outputDirDetect 'detectionLabRef.mat']);
                    detectionsLabRef=tmp.detectionsLabRef;
                else
                    error('Detection files not found. Please run detection beforehand');
                end
            end

            
            %% Registration            
            if(runRegistration)
                registerTranslation3D(MD,'show',true,'warp',true,'computeShift',true);
            end
            
            detectionsStageRef=detectionsLabRef;
            driftFilename=[MD.outputDirectory_ filesep 'regFile' filesep 'driftParameter.mat'];
            if (exist(driftFilename, 'file') == 2)
                driftParameter=load(driftFilename);
                displacements=driftParameter.displacements;
                jumpIdx=driftParameter.jumpIdx;
                parfor i=1:MD.nFrames_
                    jIdx=find((jumpIdx<i));
                    for j=jIdx
                        detectionsStageRef(i).xCoord(:,1)=detectionsStageRef(i).xCoord(:,1)+displacements{j}.T(4,1);
                        detectionsStageRef(i).yCoord(:,1)=detectionsStageRef(i).yCoord(:,1)+displacements{j}.T(4,2);
                        detectionsStageRef(i).zCoord(:,1)=detectionsStageRef(i).zCoord(:,1)+displacements{j}.T(4,3);
                    end
                end
                save([outputDirDetect filesep 'detectionStageRef.mat'],'detectionsStageRef');
                amiraWriteMovieInfo([outputDirDetect filesep 'Amira' filesep 'amiraVertexStageRef' filesep 'detectionsStageRef.am'], detectionsStageRef,'scales',dataAnisotropy);
                tmp=load([outputDirDetect filesep 'detectionStageRef.mat']);
                detectionsStageRef=tmp.detectionsStageRef;
            else
                disp('Movie has not been registered, using the laboratory frame of reference.');
            end
            
            
            
            %% Pole detection
            poleDetectionMethod=['simplex_scale_' num2str(aPoleScale,'%03d')];
            outputDirPoleDetect=[MD.outputDirectory_ filesep 'EB3' filesep 'poles' filesep poleDetectionMethod filesep];mkdir(outputDirPoleDetect);
            if(runPoleDetection)
                % Detect Poles and write results
                [poleMovieInfo]=detectPoles(MD,'showAll',false,'processFrames',processFrame,'scales',(aPoleScale*(dataAnisotropy/dataAnisotropy(1)).^(-1)));
                save([outputDirPoleDetect filesep 'poleDetection.mat'],'poleMovieInfo');
                amiraWriteMovieInfo([outputDirPoleDetect filesep 'amiraVertex' filesep poleDetectionMethod '.am'],poleMovieInfo,'scales',dataAnisotropy);
            
                %% Set detection in the spindle referential
                [dist,sphCoord,poleId,inliers,originProb,minProb,sphCoordBest,detectionsSpindleRef]=poleDist(poleMovieInfo,detectionsStageRef,'anisotropy',dataAnisotropy,'angleRef','poles');
                save([outputDirDetect filesep 'dist.mat'],'dist','minProb','poleId','inliers');
                save([outputDirDetect filesep 'sphericalCoord.mat'],'sphCoordBest');
                save([outputDirDetect filesep 'detectionSpindleRef.mat'],'detectionsSpindleRef');
 
                amiraWriteMovieInfo([outputDirDetect filesep 'Amira' filesep 'amiraVertexLabRef' filesep 'detectionLabRef.am'],detectionsLabRef, ...
                    'scales',dataAnisotropy,'prop',{{'minProb',minProb},{'azimuth',sphCoordBest.azimuth},{'elevation',sphCoordBest.elevation},{'poleId',cellfun(@(x,y) x.*y,inliers,poleId,'unif',0)}});
                amiraWriteMovieInfo([outputDirDetect filesep 'Amira' filesep 'amiraVertexSpindleRef' filesep 'detectionSpindleRef.am'],detectionsSpindleRef, ...
                    'scales',dataAnisotropy,'prop',{{'minProb',minProb},{'azimuth',sphCoordBest.azimuth},{'elevation',sphCoordBest.elevation},{'poleId',cellfun(@(x,y) x.*y,inliers,poleId,'unif',0)}});                
                if (exist([outputDirDetect filesep 'detectionStageRef.mat'], 'file') == 2)
                    amiraWriteMovieInfo([outputDirDetect filesep 'Amira' filesep  'amiraVertexStageRef' filesep 'detectionStageRef.am'],detectionsStageRef, ...
                        'scales',dataAnisotropy,'prop',{{'minProb',minProb},{'azimuth',sphCoordBest.azimuth},{'elevation',sphCoordBest.elevation},{'poleId',cellfun(@(x,y) x.*y,inliers,poleId,'unif',0)}});
                end
            else
                load([outputDirDetect filesep 'detectionSpindleRef.mat']);
                load([outputDirDetect filesep 'sphericalCoord.mat']);
                load([outputDirDetect filesep 'dist.mat']);
            end
            
            
            %% Tracking
            outputDirTrack=[MD.outputDirectory_ filesep 'EB3' filesep 'track' filesep ];
            if(trackAfterRegistration);
                outputDirTrack=[MD.outputDirectory_ filesep 'EB3' filesep 'trackSpindleRef' filesep];
            end;
            mkdir(outputDirTrack);
            if runEB3Tracking
                movieInfo=detectionsStageRef;
                if(trackAfterRegistration);
                    movieInfo=detectionsSpindleRef;
                end
                
                plusTipTrackerLegaxyRoot=[outputDirTrack filesep 'plustiptrackerio' filesep];
                mkdir([plusTipTrackerLegaxyRoot filesep  'feat' filesep]); 
                save([plusTipTrackerLegaxyRoot filesep  'feat' filesep 'movieInfo.mat'],'movieInfo');
                
                projDir=struct();
                projDir.anDir=plusTipTrackerLegaxyRoot;
                projDir.imDir='';                       % tested but useless
                
                % Additional, fixed parameter
                timeRange=[1 MD.nFrames_];
                timeWindow=2;
                plusTipCometTracker3DQD(projDir,timeWindow,...
                    minTrackLength,minRadius,maxRadius,...
                    maxFAngle,maxBAngle,maxShrinkFactor,...
                    fluctRad,timeRange,[],breakNonLinearTracks, ...
                    searchRadiusMult,searchRadiusFirstIteration);
                
                %% Convert tracks final in a user-friendlier format in three versions
                % - lab frame of reference
                % - stage frame of reference 
                % - spindle frame of reference 
                trackFile=load([plusTipTrackerLegaxyRoot filesep 'track' filesep 'trackResults.mat']);
                
                tracksLabRef=TracksHandle(trackFile.tracksFinal,detectionsLabRef);
                save([outputDirTrack filesep  'tracksLabRef.mat'],'tracksLabRef')
                tracksStageRef=TracksHandle(trackFile.tracksFinal,detectionsStageRef);
                save([outputDirTrack filesep 'tracksStageRef.mat'],'tracksStageRef')
                tracksSpindleRef=TracksHandle(trackFile.tracksFinal,detectionsSpindleRef);
                save([outputDirTrack filesep  'tracksSpindleRef.mat'],'tracksSpindleRef')
                
            else
                if (exist([outputDirTrack  filesep 'tracksLabRef.mat'], 'file') == 2)
                    load([outputDirTrack  filesep 'tracksLabRef.mat']);
                    load([outputDirTrack  filesep 'tracksStageRef.mat']);
                    load([outputDirTrack  filesep 'tracksSpindleRef.mat']);
                end
            end
            
            
            if(runAmiraWrite)
                
                %% Tracks in the lab FoR.
                amiraWriteTracks([outputDirTrack filesep 'Amira' filesep 'AmiraTrackLabRef' filesep 'tracksLabRef.am'],tracksLabRef,'MD',MD);
                amiraWriteTracks([outputDirTrack filesep 'Amira' filesep  'AmiraTrackLabRef20plus' filesep 'trackLabRef20plus.am'],tracksLabRef([tracksLabRef.lifetime]>20),'MD',MD);

                %% Tracks in the Stage FoR.
                amiraWriteTracks([outputDirTrack filesep 'Amira' filesep  'AmiraTracksStageRef' filesep 'tracksLabRef.am'],tracksStageRef,'MD',MD);
                amiraWriteTracks([outputDirTrack filesep 'Amira' filesep  'AmiraTrackStageRef20plus' filesep 'trackLabRef20plus.am'],tracksStageRef([tracksStageRef.lifetime]>20),'MD',MD);                
                
                %% Tracks in the spindle FoR. 
                amiraWriteTracks([outputDirTrack filesep 'Amira' filesep  'AmiraTrackSpindleRef' filesep 'trackSpindleRef.am'],tracksSpindleRef,'MD',MD);
                amiraWriteTracks([outputDirTrack filesep 'Amira' filesep  'AmiraTrackSpindleRef20plus' filesep 'trackSpindleRef20plus.am'],tracksSpindleRef([tracksSpindleRef.lifetime]>20),'MD',MD);

            end
            
            if(runPostProcessing)
                [sphericalAzimuth,sphericalElevation,time,trackId,poleId]=sphericalDistribution(tracksStageRef,sphCoordBest.azimuth,sphCoordBest.elevation,sphCoordBest.rho,poleId,sphericalProjectionRadius);
                mkdir([MD.outputDirectory_ filesep 'EB3' filesep 'sphericalProjection' filesep 'radius-' num2str(sphericalProjectionRadius)]);
                save([MD.outputDirectory_ filesep 'EB3' filesep 'sphericalProjection' filesep 'radius-' num2str(sphericalProjectionRadius) filesep 'sphericalProjection.mat'],'sphericalAzimuth','sphericalElevation','time','trackId','poleId');                
            end;
            
        end
    end
end