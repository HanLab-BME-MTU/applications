function batchEB3analysis
% BATCHEB3ANALYSIS allows user to detect, track, and analyze EB comets for all movies in user-selected projList

% Before running this function,
% 1) Set the parameters below.
% 2) Make sure you know where the appropriate projList.mat variable is
%    saved. projList is created by the function getProj, which saves projList
%    in the user-selected top-level directory.  You will be asked to select
%    this variable when you run the batch.  It contains the list of
%    projects you want to analyze.

% USER-SET PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Which parts of the analysis do you want to run?

% use 0 to run detection, 1 to skip
doDetect=0;

% use 0 to run tracking, 1 to skip
doTrack=0; 

% use 0 to run post-processing function metaEB3analysis, 1 to skip
doMeta=0;


% DETECTION parameters (will only matter if doDetect=1) %%%%%%%%%%%%%%%%%%%

timeRangeDetect = []; % frame range [start end], use [] to use all frames
bitDepth  = 14; % change according to your camera (should be 12, 14, or 16)
savePlots = 1;  % 1 to save overlay plots of detection results in subfolder; 0 if not (may run faster)

% TRACKING parameters
% open file using "edit scriptTrack_EB3" and check tracking parameters!!!

% META parameters (will only matter if doMeta=1) %%%%%%%%%%%%%%%%%%%%%%%%%%
% (i.e. k: 2.0s,105nm;  y: 2.0s,84nm;  c: 0.8s,110nm)
secPerFrame=2; % frame rate
pixSizeNm=84; % real-space pixel size (nanometers)

doPopHist=1; % 1 to make population histograms, 0 to skip

doFeatVelMovie=0; % 1 to make velocity movie, 0 to skip
timeRangeMovie = [1 30]; % frame range [start end], use [] to use all frames
velLimit  = [30]; % max speed (microns/minute) to use for color min/max, use [] for full range
% (i.e. all tracks faster than velLimit will be the same shade of red; 
% all shrinkage events faster than -velLimit will be the same shade of blue) 

% END OF USER-SET PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



close all

[fileName,pathName] = uigetfile('*.mat','Select projList.mat file containing projects to run');
load(formatPath([pathName filesep fileName]));
clear 'fileName' 'pathName'

for i=1:size(projList,1)
    
    if ~isempty(strfind(projList(i).anDir,'sub'))
        continue
    end
    
    %try
        
        if doDetect==1
            disp(['detecting project ' num2str(i) ': ' projList(i).anDir])
            [movieInfo]=eb1SpotDetector(projList(i),timeRangeDetect,bitDepth,savePlots);
        end
        
%         if changeMovieInfo==1
%             homeDir=pwd;
%             [tempDir] = formatPath(projList(i).anDir);
%             cd([tempDir filesep 'feat'])
%             load(['movieInfo.mat'])
% 
%             for j=1:length(movieInfo)
%                 movieInfo(j,1).xCoord(:,2)=.5;
%                 movieInfo(j,1).yCoord(:,2)=.5;
%             end
%             save('movieInfo.mat','movieInfo')
%             cd(homeDir)
% 
%         end
        
        if doTrack==1
            homeDir=pwd;
            [tempDir] = formatPath(projList(i).anDir);
            cd(tempDir)
            disp(['tracking project ' num2str(i) ': ' projList(i).anDir])
            scriptTrack_EB3 
            cd(homeDir)
        end
        
        if doMeta==1
            [temp.anDir] = formatPath(projList(i).anDir);
            [temp.imDir] = formatPath(projList(i).imDir);

            [projData]=metaEB3analysis(temp,secPerFrame,pixSizeNm);
        
            if doPopHist==1
                popHist(projData);
            end
            if doFeatVelMovie==1
                roiYX=imread([projData.anDir filesep 'roiMask.tif']);
                featVelMovie(projData,timeRangeMovie,velLimit,roiYX);
            end
        end
       
%     catch
%         disp('problem with')
%         projList(i).anDir
%     
%     end
end

close all







