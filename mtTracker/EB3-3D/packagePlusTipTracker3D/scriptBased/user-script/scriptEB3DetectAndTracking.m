%%% Description
% This script detect and track EB3s in list of movies. The resulting tracks
% of each movies are saved in their respective "analysis" in the original
% u-track format, a more use-friendly structured format and an
% Amira-readable format for visualization. If required, stage shift can be
% registered automatically to correct detection coordinates.

%% %% USER INPUT 

%% INPUT MOVIES
% After movie loading (see loadMoviesManagementFile.m). provide an array of
% MovieList. 
aMovieListArray=[ inter10nML inter33nML inter100nML];%interControlNoco

%% pixelSizes (while movieData does not work properly)
pixelSizeLat=100;
pixelSizeAxial=235;

% Frame to be processed ( [] for whole movie)
processFrame=[];


%% PROCESS MANAGEMENT 
%  - 1 to (re-)run the algorithm
%  - 0 to load previously computed results (if any)
runEB3Detection=1;     
runPoleDetection=1; 
runRegistration=0;  % Register sudden stage shift during acquisition.
runEB3Tracking=1;      

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
        for i=1:length(ML.movieDataFile_)            
            MD=ML.getMovie(i);
            disp(['Processing movie: ' MD.getFullPath]);
            MD.sanityCheck;
            disp(['Processing movie: ' MD.getFullPath]);
            MD.pixelSize_=pixelSizeLat;
            MD.pixelSizeZ_=pixelSizeAxial;
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
            if(runPoleDetection)
                [poleMovieInfo]=detectPoles(MD,'showAll',false,'processFrames',processFrame,'scales',(aPoleScale*(dataAnisotropy/dataAnisotropy(1)).^(-1)));
                save([outputDirPoleDetect filesep 'detection.mat'],'poleMovieInfo');
                mkdir([outputDirPoleDetect filesep 'amiraVertex']);
                amiraWriteMovieInfo([outputDirPoleDetect filesep 'amiraVertex' filesep poleDetectionMethod '.am'],poleMovieInfo,'scales',dataAnisotropy);

            
                
                outputDirDist=[MD.outputDirectory_ filesep 'EB3PoleRef' filesep poleDetectionMethod filesep detectionMethod];mkdir(outputDirDist);
                [dist,poleId,inliers,originProb,minProb,azimuth,elevation,rho]=poleDist(poleMovieInfo,movieInfo,'anisotropy',dataAnisotropy,'angleRef','poles');
                save([outputDirDist filesep 'dist.mat'],'dist','minProb','poleId','inliers');
                save([outputDirDist filesep 'sphericalCoord.mat'],'azimuth','elevation','rho');
                mkdir([outputDirDist filesep 'amiraVertex']);
                amiraWriteMovieInfo([outputDirDist filesep 'amiraVertex' filesep poleDetectionMethod '_' detectionMethod '.am'],movieInfo, ...
                    'scales',dataAnisotropy,'prop',{{'minProb',minProb},{'azimuth',azimuth},{'elevation',elevation},{'poleId',cellfun(@(x,y) x.*y,inliers,poleId,'unif',0)}});
            else
                load([outputDirPoleDetect filesep 'detection.mat']);
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
            if runEB3Tracking
                mkdir(outputDirTrack);
                mkdir([outputDirTrack '/feat']);
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
                tmp=load([outputDirTrack filesep 'track' filesep 'trackResults.mat']);
                tracks=TracksHandle(tmp.tracksFinal);
                save([outputDirTrack filesep 'track' filesep 'trackNewFormat.mat'],'tracks')
                
                %% load track results and save them to Amira
                mkdir([outputDirTrack filesep 'AmiraTrack']);
                amiraWriteTracks([outputDirTrack filesep 'AmiraTrack' filesep 'test.am'],tracks,'scales',dataAnisotropy);
                
            end
            
            
        end
    end
end