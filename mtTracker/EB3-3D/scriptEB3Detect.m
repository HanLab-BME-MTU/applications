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
aMovieListArray=movieListArray;%interControlNoco

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
waterThresh=80; % Only for the bandPassWatershed method


%% %% Process each movie
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
    end
end
