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
% An other example of usage with movieList in different folder
%MLPath='/work/gdanuser/proudot/project/EB3-3D-track/packaging/alpha/plusTipTracker3D-alpha-1/data/A1_HeLa_Cells_EB1/'
MLPath='/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/smallData/prometaphase/analysis/';
movieListFileNames={'Cell2.mat'};
% Build the array of MovieList (automatic)
aMovieListArray=cellfun(@(x) MovieList.loadMatFile([MLPath filesep x]), movieListFileNames,'unif',0);
aMovieListArray=aMovieListArray{:};

%% ALTERNATIVELY, an array of movieList object are already loaded can be used (see loadMoviesManagementFile.m).
%  aMovieListArray=[ inter10nML inter33nML inter100nML];

%% PROCESS MANAGEMENT 
%  - 1 to (re-)run the algorithm
%  - 0 to load previously computed results (if any)
runEB3Detection=0;  % Detect EB3 "comets"   
runPoleDetection=0; % Detec spindle Poles
runRegistration=0;  % Register sudden stage shift during acquisition.
runEB3Tracking=1; 
runPostProcessing=0;      
% Frame to be processed ([] for whole movie)
processFrame=[];

%% Optional output
printAmiraFile=1;      
printDetectionMask=0;

%% ALGORITHM PARAMETERS

%% EB3 Detection parameters
% The detection method must be chosen depending on the data at hand
% - 'pointSourceAutoSigmaLM' is parameter-free and suitable for raw 
%   (but deskewed) data
% - 'bandPassWatershed' is more suitable for deconvolve data, requires 
%   input threshold parameter <waterThresh>.
detectionMethod='pointSourceAutoSigmaLM';
scales=[];
scales=[1.4396 1.2913];
%scales=[2 1.5];
alpha=0.05;
waterThresh=80; % Only for bandPassWatershed

%% Pole detection parameters
poleScale=3;  % expected scale of the poles, should not need to be changed.


%% Tracking parameters
% Segments presenting a lifetime below <minTrackLength> are discarded. 
minTrackLength=2; % In frames.
 
% Lower and upper bound for the radius of the sphere considered to link 
% a detection one frame to the next
minRadius=2; % In pixels.
maxRadius=8; % In pixels.
searchRadiusMult=3;
searchRadiusFirstIteration=10;

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

timeWindow=2;

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
            dataIsotropy=[MD.pixelSize_ MD.pixelSize_ MD.pixelSize_];
            
            %% EB3 Detection
            outputDirDetect=[MD.outputDirectory_ filesep 'EB3' filesep 'detection' filesep];mkdir(outputDirDetect);
            if(runEB3Detection)
                [detectionsLabRef,lab]=detectEB3(MD,'type',detectionMethod,'waterThresh',waterThresh,'showAll',false,'scales',scales,'Alpha',alpha);
                                
                for fIdx=1:length(detectionsLabRef)
                    detectionsLabRef(fIdx).zCoord(:,1)=detectionsLabRef(fIdx).zCoord(:,1)*MD.pixelSizeZ_/MD.pixelSize_;
                end

                save([outputDirDetect filesep 'detectionLabRef.mat'],'detectionsLabRef');
                
                if(printAmiraFile)
                    amiraWriteMovieInfo([outputDirDetect filesep 'Amira' filesep 'amiraVertexLabRef' filesep  'detectionLabRef.am'], detectionsLabRef,'scales',dataIsotropy);
                end
                
                if(printDetectionMask)
                    mkdir([outputDirDetect filesep 'detectionMaskLabRef']);
                    for tidx=1:length(detectionsLabRef)
                        stackWrite(lab{tidx},[outputDirDetect filesep 'detectionMaskLabRef' filesep 'detect_T_' num2str(tidx,'%05d') '.tif']);
                    end
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
                if(printAmiraFile)
                    amiraWriteMovieInfo([outputDirDetect filesep 'Amira' filesep 'amiraVertexStageRef' filesep 'detectionsStageRef.am'], detectionsStageRef,'scales',dataIsotropy);
                end
                tmp=load([outputDirDetect filesep 'detectionStageRef.mat']);
                detectionsStageRef=tmp.detectionsStageRef;
            else
                disp('Movie has not been registered, using the laboratory frame of reference.');
            end
            
            
            
            %% Pole detection
            poleDetectionMethod=['simplex_scale_' num2str(aPoleScale,'%03d')];
            outputDirPoleDetect=[MD.outputDirectory_ filesep 'EB3' filesep 'poles' filesep poleDetectionMethod filesep];mkdir(outputDirPoleDetect);
            spindleReferential=1;
            if(runPoleDetection)
                % Detect Poles and write results
                [poleMovieInfo]=detectPoles(MD,'showAll',false,'processFrames',processFrame,'scales',(aPoleScale*(dataAnisotropy/dataAnisotropy(1)).^(-1)));
                                % Isotropic coordinate (
                for fIdx=1:length(poleMovieInfo)
                    poleMovieInfo(fIdx).zCoord(:,1)=poleMovieInfo(fIdx).zCoord(:,1)*MD.pixelSizeZ_/MD.pixelSize_;
                end
                save([outputDirPoleDetect filesep 'poleDetection.mat'],'poleMovieInfo');
                if(printAmiraFile)
                    amiraWriteMovieInfo([outputDirPoleDetect filesep 'amiraVertex' filesep poleDetectionMethod '.am'],poleMovieInfo,'scales',dataIsotropy);
                end
                
                
                %% Set detection in the spindle referential
                [dist,sphCoord,poleId,inliers,originProb,minProb,sphCoordBest,detectionsSpindleRef]=poleDist(poleMovieInfo,detectionsStageRef,'anisotropy',dataIsotropy,'angleRef','poles');
                save([outputDirDetect filesep 'dist.mat'],'dist','minProb','poleId','inliers');
                save([outputDirDetect filesep 'sphericalCoord.mat'],'sphCoordBest','sphCoord');
                save([outputDirDetect filesep 'detectionSpindleRef.mat'],'detectionsSpindleRef');
                
                if(printAmiraFile)
                    amiraWriteMovieInfo([outputDirDetect filesep 'Amira' filesep 'amiraVertexLabRef' filesep 'detectionLabRef.am'],detectionsLabRef, ...
                        'scales',dataIsotropy,'prop',{{'minProb',minProb},{'azimuth',sphCoordBest.azimuth},{'elevation',sphCoordBest.elevation},{'poleId',cellfun(@(x,y) x.*y,inliers,poleId,'unif',0)}});
                    amiraWriteMovieInfo([outputDirDetect filesep 'Amira' filesep 'amiraVertexSpindleRef' filesep 'detectionSpindleRef.am'],detectionsSpindleRef, ...
                        'scales',dataIsotropy,'prop',{{'minProb',minProb},{'azimuth',sphCoordBest.azimuth},{'elevation',sphCoordBest.elevation},{'poleId',cellfun(@(x,y) x.*y,inliers,poleId,'unif',0)}});
                    if (exist([outputDirDetect filesep 'detectionStageRef.mat'], 'file') == 2)
                        amiraWriteMovieInfo([outputDirDetect filesep 'Amira' filesep  'amiraVertexStageRef' filesep 'detectionStageRef.am'],detectionsStageRef, ...
                            'scales',dataIsotropy,'prop',{{'minProb',minProb},{'azimuth',sphCoordBest.azimuth},{'elevation',sphCoordBest.elevation},{'poleId',cellfun(@(x,y) x.*y,inliers,poleId,'unif',0)}});
                    end
                end
            else
                if (exist([outputDirDetect filesep 'detectionSpindleRef.mat'], 'file') == 2)
                    load([outputDirDetect filesep 'detectionSpindleRef.mat']);
                    load([outputDirDetect filesep 'sphericalCoord.mat']);
                    load([outputDirDetect filesep 'dist.mat']);
                else
                    disp('No pole detection previously executed. Disabling spindle referential.');
                    trackAfterRegistration=0;
                    spindleReferential=0;
                end

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
                

                plusTipTrackerLegacyRoot=[outputDirTrack filesep 'plustiptrackerio' filesep];
                mkdir([plusTipTrackerLegacyRoot filesep  'feat' filesep]); 
                save([plusTipTrackerLegacyRoot filesep  'feat' filesep 'movieInfo.mat'],'movieInfo');
                
                projDir=struct();
                projDir.anDir=plusTipTrackerLegacyRoot;
                projDir.imDir='';                       % tested but useless
                
                % Additional, fixed parameter
                timeRange=[1 MD.nFrames_];
                
                plusTipCometTracker3DQD(projDir,timeWindow,...
                    minTrackLength,minRadius,maxRadius,...
                    maxFAngle,maxBAngle,maxShrinkFactor,...
                    fluctRad,timeRange,[],breakNonLinearTracks, ...
                    searchRadiusMult,searchRadiusFirstIteration);
                
                %% Convert tracks final in a user-friendlier format in three versions
                % - lab frame of reference
                % - stage frame of reference 
                % - spindle frame of reference 
                trackFile=load([plusTipTrackerLegacyRoot filesep 'track' filesep 'trackResults.mat']);
                %%
                tracksFinalStripped=rmfield(trackFile.tracksFinal,'tracksCoordAmpCG');
                tracksLabRef=TracksHandle(tracksFinalStripped,detectionsLabRef);
                save([outputDirTrack filesep  'tracksLabRef.mat'],'tracksLabRef')
                tracksFinalStripped=rmfield(trackFile.tracksFinal,'tracksCoordAmpCG');
                tracksStageRef=TracksHandle(tracksFinalStripped,detectionsStageRef);
                save([outputDirTrack filesep 'tracksStageRef.mat'],'tracksStageRef')
                if(spindleReferential)
                    tracksFinalStripped=rmfield(trackFile.tracksFinal,'tracksCoordAmpCG');
                    tracksSpindleRef=TracksHandle(tracksFinalStripped,detectionsSpindleRef);
                    save([outputDirTrack filesep  'tracksSpindleRef.mat'],'tracksSpindleRef')
                end
                
            else
                if (exist([outputDirTrack  filesep 'tracksLabRef.mat'], 'file') == 2)
                    load([outputDirTrack  filesep 'tracksLabRef.mat']);
                    load([outputDirTrack  filesep 'tracksStageRef.mat']);
                    if(spindleReferential)
                        load([outputDirTrack  filesep 'tracksSpindleRef.mat']);
                    end
                end
            end
            
            
            if(printAmiraFile)
                
                %% Tracks in the lab FoR.
                amiraWriteTracks([outputDirTrack filesep 'Amira' filesep 'AmiraTrackLabRef' filesep 'tracksLabRef.am'],tracksLabRef,'scales',dataIsotropy);
                amiraWriteTracks([outputDirTrack filesep 'Amira' filesep  'AmiraTrackLabRef20plus' filesep 'trackLabRef20plus.am'],tracksLabRef([tracksLabRef.lifetime]>20),'MD',MD);

                %% Tracks in the Stage FoR.
                amiraWriteTracks([outputDirTrack filesep 'Amira' filesep  'AmiraTracksStageRef' filesep 'tracksLabRef.am'],tracksStageRef,'scales',dataIsotropy);
                amiraWriteTracks([outputDirTrack filesep 'Amira' filesep  'AmiraTrackStageRef20plus' filesep 'trackLabRef20plus.am'],tracksStageRef([tracksStageRef.lifetime]>20),'MD',MD);                
                
                %% Tracks in the spindle FoR. 
                if(spindleReferential)
                    amiraWriteTracks([outputDirTrack filesep 'Amira' filesep  'AmiraTrackSpindleRef' filesep 'trackSpindleRef.am'],tracksSpindleRef,'scales',dataIsotropy);
                    amiraWriteTracks([outputDirTrack filesep 'Amira' filesep  'AmiraTrackSpindleRef20plus' filesep 'trackSpindleRef20plus.am'],tracksSpindleRef([tracksSpindleRef.lifetime]>20),'MD',MD);
                end

            end
            
            if(runPostProcessing)
                [sphericalAzimuth,sphericalElevation,time,trackId,poleId]=sphericalDistribution(tracksStageRef,sphCoordBest.azimuth,sphCoordBest.elevation,sphCoordBest.rho,poleId,sphericalProjectionRadius);
                mkdir([MD.outputDirectory_ filesep 'EB3' filesep 'sphericalProjection' filesep 'radius-' num2str(sphericalProjectionRadius)]);
                save([MD.outputDirectory_ filesep 'EB3' filesep 'sphericalProjection' filesep 'radius-' num2str(sphericalProjectionRadius) filesep 'sphericalProjection.mat'],'sphericalAzimuth','sphericalElevation','time','trackId','poleId');                
            end;
            
        end
    end
end