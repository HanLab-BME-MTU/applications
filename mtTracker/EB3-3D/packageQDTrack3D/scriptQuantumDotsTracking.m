%%% Description
% This script detect and track detections in list of movies. The resulting tracks
% of each movies are saved in their respective "analysis" in the original
% u-track format, a more use-friendly structured format and an
% Amira-readable format for visualization.

%% %% USER INPUT 

%% Load the movie list files (without checking channels)
MLPath='/project/bioinformatics/Danuser_lab/shared/proudot/3d-vis/utrackPackage/scriptBased/UTrack-QD-v1/data/analysis'
movieListFileNames={'singleCellList.mat'};

aMovieListArray=cellfun(@(x) MovieList.loadMatFile([MLPath filesep x]), movieListFileNames,'unif',0);
aMovieListArray=aMovieListArray{:};

%% Spatiotemporal ROI
% Frame to be processed ( [] for whole movie)
processFrames=[];
% ROI mask name
maskFileName='';

%% PROCESS MANAGEMENT 
%  - 1 to (re-)run the algorithm
%  - 0 to load previously computed results (if any)
runDetection=0;     
runTracking=1;      
runAmiraRendering=1;

%% ALGORITHM PARAMETERS

%% Detection parameters
% The detection method must be chosen depending on the data at hand
% - 'pointSourceAutoSigmaLM' is parameter-free and suitable for raw 
%   (but deskewed) data
% - 'bandPassWatershed' is more suitable for deconvolve data, requires 
%   input threshold parameter <waterThresh>.
detectionMethod='pointSourceAutoSigmaFit';
alpha=0.01;

%% Tracking parameters


for k=1:length(aMovieListArray)
    ML=aMovieListArray(k);
    
    %% loop over each different cell in each condition
    for i=1:length(ML.movieDataFile_)
        MD=MovieData.loadMatFile(ML.movieDataFile_{i});
        disp(['Processing movie: ' MD.getFullPath]);
        dataAnisotropy=[MD.pixelSize_ MD.pixelSize_ MD.pixelSizeZ_];
        
        % Read mask
        mask=[];
        if(exist([MD.outputDirectory_ filesep 'mask' filesep maskFileName],'file'))
            mask=readtiff([MD.outputDirectory_ filesep 'mask' filesep maskFileName]);
        end
        
        %% Detection     
        outputDirDetect=[MD.outputDirectory_ filesep 'detection' filesep detectionMethod];mkdir(outputDirDetect);
        if(runDetection)
            [movieInfo,lab,scales]=detectEB3(MD,'type',detectionMethod,'Alpha',alpha,'processFrames',processFrames, ... 
                'printAll',true,'showAll',false,'subDirectory','detection', 'mask',mask,'scales',[1.3141 0.8950]);
        else
            load([outputDirDetect filesep 'detection.mat']);
        end
        
        
        %% Tracking
        outputDir=[outputDirDetect filesep 'tracks'];mkdir(outputDir);
        if runTracking
            [gapCloseParam,costMatrices,kalmanFunctions,...
                probDim,verbose]=QDTrackerParam();
            watch_KF_iter=0;
            
            saveResults.dir =  outputDir; %directory where to save input and output
            saveResults.filename = 'trackResults.mat'; %name of file where input and output are saved
        
            
            % Correct for zCoord scale while keeping pixel values (to be
            % reflected in later scale. 
            zRatio=MD.pixelSizeZ_/MD.pixelSize_;
            for i=1:length(movieInfo) movieInfo(i).zCoord=movieInfo(i).zCoord*zRatio; end;
            
            [tracksFinal,kalmanInfoLink,errFlag] = ...
                trackCloseGapsKalmanSparse(movieInfo, ...
                costMatrices,gapCloseParam,kalmanFunctions,...
                probDim,saveResults,verbose);
      
            tracks=TracksHandle(tracksFinal);
            save([outputDir filesep 'tracksHandle.mat'],'tracks');
        end
        
        if(runAmiraRendering)
            load([outputDir filesep 'tracksHandle.mat']);
            s=[MD.pixelSize_ MD.pixelSize_ MD.pixelSize_ MD.timeInterval_];
            
            amiraWriteTracks([outputDir filesep 'AmiraTracks' filesep 'tracking_' detectionMethod '.am'],tracks,'scales',s);
        end
        
    end
end