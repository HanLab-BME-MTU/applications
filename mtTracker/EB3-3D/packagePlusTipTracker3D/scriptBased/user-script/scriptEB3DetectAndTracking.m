%%% Description
% This script detect and track EB3s in list of movies. The resulting tracks
% of each movies are saved in their respective "analysis" in the original
% u-track format, a more use-friendly structured format and an
% Amira-readable format for visualization. If required, stage shift can be
% registered automatically to correct detection coordinates.

%% %% USER INPUT 

%% LOADING MOVIES INFORMATION (TWO OPTIONS)
%% Loading MovieList file
% MovieList paths: 
%MLPath='/work/gdanuser/proudot/project/EB3-3D-track/packaging/alpha/plusTipTracker3D-alpha-1/data/A1_HeLa_Cells_EB1/prometaphase/analysis/'
MLPath='/project/cellbiology/gdanuser/shared/proudot/project/lattice-track/data-analysis/chemical-inhibition/2015_06_17(EB1_GFP_drugs)/analysis/';

% MovieList FileName (the combination of condition you want to compare). 
%movieListFileNames={'movieList.mat'};
movieListFileNames={'controlAnaphase.mat'};
%movieListFileNames={'controlMetaphase.mat','10nMMetaphase.mat'};
%movieListFileNames={'controlMetaphase.mat','10nMMetaphase.mat','33nMMetaphase.mat','100nMMetaphase.mat'};
%movieListFileNames={'controlInterphase.mat','10nMInterphase.mat','33nMInterphase.mat','100nMInterphase.mat',};
%movieListFileNames={'controlAnaphase.mat','10nMAnaphase.mat','33nMAnaphase.mat'};

% Build the array of MovieList (automatic)
aMovieListArray=cellfun(@(x) MovieList.loadMatFile([MLPath x]), movieListFileNames,'unif',0);
aMovieListArray=aMovieListArray{:};

%% Using movieList object that are already loaded (see loadMoviesManagementFile.m).
%  aMovieListArray=[ inter10nML inter33nML inter100nML];

% Frame to be processed ( [] for whole movie)
processFrame=[];

%% PROCESS MANAGEMENT 
%  - 1 to (re-)run the algorithm
%  - 0 to load previously computed results (if any)
runEB3Detection=0;     
runPoleDetection=0; 
runRegistration=0;  % Register sudden stage shift during acquisition.
runEB3Tracking=1; runAmiraWrite=1;      

%% DEBUGGING PARAMETER (Could void your waranty)
trackAfterRegistration=1;

%% ALGORITHM PARAMETERS

%% EB3 Detection parameters
% The detection method must be chosen depending on the data at hand
% - 'pointSourceAutoSigmaLM' is parameter-free and suitable for raw 
%   (but deskewed) data
% - 'bandPassWatershed' is more suitable for deconvolve data, requires 
%   input threshold parameter <waterThresh>.
detectionMethod='pointSourceAutoSigmaLM';
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

%% %% Process each movie
for aPoleScale=poleScale
    for k=1:length(aMovieListArray)
        ML=aMovieListArray(k);
        
        %% loop over each different cell in each condition
        for i=1:1%length(ML.movieDataFile_)            
            MD=MovieData.loadMatFile(ML.movieDataFile_{i});
            disp(['Processing movie: ' MD.getFullPath]);
            dataAnisotropy=[MD.pixelSize_ MD.pixelSize_ MD.pixelSizeZ_];
            
            %% EB3 Detection
            outputDirDetect=[MD.outputDirectory_ filesep 'EB3' filesep detectionMethod];mkdir(outputDirDetect);
            if(runEB3Detection)
                [movieInfo,lab]=detectEB3(MD,'type',detectionMethod,'waterThresh',waterThresh,'showAll',false)
                save([outputDirDetect filesep 'detection.mat'],'movieInfo');
                mkdir([outputDirDetect filesep 'amiraVertex' ]);
                amiraWriteMovieInfo([outputDirDetect filesep 'amiraVertex' filesep detectionMethod '.am'], movieInfo,'scales',dataAnisotropy);
                mkdir([outputDirDetect filesep 'mask']);
                for tidx=1:length(movieInfo)
                    stackWrite(lab{tidx},[outputDirDetect filesep 'mask' filesep 'detect_T_' num2str(tidx,'%05d') '.tif']);
                end
                mkdir([outputDirDetect filesep 'amiraVertex']);
            else
                load([outputDirDetect filesep 'detection.mat']);
            end

            
            %% Pole detection
            poleDetectionMethod=['simplex_scale_' num2str(aPoleScale,'%03d')];
            outputDirPoleDetect=[MD.outputDirectory_ filesep 'poles' filesep poleDetectionMethod];mkdir(outputDirPoleDetect);
            outputDirDist=[MD.outputDirectory_ filesep 'EB3PoleRef' filesep poleDetectionMethod filesep detectionMethod];mkdir(outputDirDist);
            if(runPoleDetection)
                % Detect Poles and write results
                [poleMovieInfo]=detectPoles(MD,'showAll',false,'processFrames',processFrame,'scales',(aPoleScale*(dataAnisotropy/dataAnisotropy(1)).^(-1)));
                save([outputDirPoleDetect filesep 'poleDetection.mat'],'poleMovieInfo');
                mkdir([outputDirPoleDetect filesep 'amiraVertex']);
                amiraWriteMovieInfo([outputDirPoleDetect filesep 'amiraVertex' filesep poleDetectionMethod '.am'],poleMovieInfo,'scales',dataAnisotropy);

            
                % Set detection in the spindle referential
                [dist,poleId,inliers,originProb,minProb,azimuth,elevation,rho,movieInfoSpindle]=poleDist(poleMovieInfo,movieInfo,'anisotropy',dataAnisotropy,'angleRef','poles');
                save([outputDirDist filesep 'dist.mat'],'dist','minProb','poleId','inliers');
                save([outputDirDist filesep 'sphericalCoord.mat'],'azimuth','elevation','rho');
                save([outputDirDist filesep 'cartesianCoord.mat'],'movieInfoSpindle');
                mkdir([outputDirDist filesep 'amiraVertex']);
                amiraWriteMovieInfo([outputDirDist filesep 'amiraVertex' filesep poleDetectionMethod '_' detectionMethod '.am'],movieInfo, ...
                    'scales',dataAnisotropy,'prop',{{'minProb',minProb},{'azimuth',azimuth},{'elevation',elevation},{'poleId',cellfun(@(x,y) x.*y,inliers,poleId,'unif',0)}});
                amiraWriteMovieInfo([outputDirDist filesep 'amiraVertexSpindle' filesep poleDetectionMethod '_' detectionMethod '.am'],movieInfoSpindle, ...
                    'scales',dataAnisotropy,'prop',{{'minProb',minProb},{'azimuth',azimuth},{'elevation',elevation},{'poleId',cellfun(@(x,y) x.*y,inliers,poleId,'unif',0)}});
            else
                load([outputDirPoleDetect filesep 'detection.mat']);
                load([outputDirDist filesep 'cartesianCoord.mat']);
            end
            

            
            %% Registration
            if(runRegistration)
                [~,jumpIdx,displacements]=registerTranslation3D(MD,'show',true,'warp',true,'computeShift',true);
                parfor i=1:MD.nFrames_
                    jIdx=find((jumpIdx<i));
                    for j=jIdx
                        movieInfo(i).xCoord(:,1)=movieInfo(i).xCoord(:,1)+displacements{j}.T(4,1);
                        movieInfo(i).yCoord(:,1)=movieInfo(i).yCoord(:,1)+displacements{j}.T(4,2);
                        movieInfo(i).zCoord(:,1)=movieInfo(i).zCoord(:,1)+displacements{j}.T(4,3);
                    end
                end
                save([outputDirDetect filesep 'detectionReg.mat'],'movieInfo');
            else
                if (exist([outputDirDetect filesep 'detectionReg.mat'], 'file') == 2)
                    load([outputDirDetect filesep 'detectionReg.mat']);
                end
            end
            
            
            %% Tracking
            outputDirTrack=[outputDirDetect filesep 'plustipTrackerio'];
            if(trackAfterRegistration);
                            outputDirTrack=[outputDirDist filesep 'plustipTrackerio'];
            end;
            if runEB3Tracking
                mkdir(outputDirTrack);
                mkdir([outputDirTrack '/feat']);
                
                if(trackAfterRegistration);
                    movieInfo=movieInfoSpindle;
                end
                
                save([outputDirTrack '/feat/movieInfo.mat'],'movieInfo');
                
      
                projDir=struct();
                projDir.anDir=outputDirTrack;
                projDir.imDir='';                       % tested but useless
                
                % Additional, fixed parameter
                timeRange=[1 MD.nFrames_];
                timeWindow=2;
                plusTipCometTracker3DQD(projDir,timeWindow,...
                    minTrackLength,minRadius,maxRadius,...
                    maxFAngle,maxBAngle,maxShrinkFactor,...
                    fluctRad,timeRange,[],breakNonLinearTracks, ...
                    searchRadiusMult,searchRadiusFirstIteration);
                
                %% Convert tracks final in a user-friendlier format
                trackFile=load([outputDirTrack filesep 'track' filesep 'trackResults.mat']);
                tracks=TracksHandle(trackFile.tracksFinal);
                save([outputDirTrack filesep 'track' filesep 'trackNewFormat.mat'],'tracks')
                
            end
            if(runAmiraWrite)
                trackFile=load([outputDirTrack filesep 'track' filesep 'trackResults.mat']);
                               
                tmp=load([outputDirDetect filesep 'detection.mat']);
                movieInfoLab=tmp.movieInfo;
                
                %% Tracks in the lab FoR.
                tracks=TracksHandle(trackFile.tracksFinal,movieInfoLab);
                amiraWriteTracks([outputDirTrack filesep 'AmiraTrackLabRef' filesep 'trackLabRef.am'],tracks,'MD',MD);
                amiraWriteTracks([outputDirTrack filesep 'AmiraTrackLabRef20plus' filesep 'trackLabRef20plus.am'],tracks([tracks.lifetime]>20),'MD',MD);

                %% Tracks in the spindle FoR. 
                tracks=TracksHandle(trackFile.tracksFinal,movieInfoSpindle);
                amiraWriteTracks([outputDirTrack filesep 'AmiraTrackSpindleRef' filesep 'trackSpindleRef.am'],tracks,'MD',MD);
                amiraWriteTracks([outputDirTrack filesep 'AmiraTrackSpindleRef20plus' filesep 'trackSpindleRef20plus.am'],tracks([tracks.lifetime]>20),'MD',MD);

            end
            
            
        end
    end
end