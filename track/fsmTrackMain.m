function [fsmParam,status]=fsmTrackMain(fsmParam)
% fsmTrackMain is the main function of the fsmTrack module
%
% SYNOPSIS   [fsmParam,status]=fsmTrackMain(fsmParam)
%
% INPUT      fsmParam:   general parameter structure
%
% OUTPUT     fsmParam:   modified (when needed) parameter structure
%            status  :   it decides whether the program should continue after this module 
%                        or should stop because of errors;
%                        status is set to 0 (error) in the beginning of a module and will
%                        be set to 1 at the end if the module completed successfully.
%                        status = 1 - if the module completed successfully
%                        status = 0 - if the module did not complete successfully
%
% DEPENDENCES   fsmTrackMain uses {}
%               fsmTrackMain is used by { fsmMain }
%
% Aaron Ponti, October 7th, 2002


% Set initial module status
status=0;

% Check input parameter
if nargin~=1
    error('Input parameter fsmParam expected');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% READ NEEDED PARAMETERS FROM fsmParam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

threshold=fsmParam.track.threshold;
influence=fsmParam.track.influence;
userPath=fsmParam.main.path;
TRACKER=fsmParam.track.tracker;
enhanced=fsmParam.track.enhanced;
grid=fsmParam.track.grid;
strg=fsmParam.specific.formString;
n=fsmParam.specific.imageNumber;
imgSize=fsmParam.specific.imgSize;
firstIndex=fsmParam.specific.firstIndex;
d0=fsmParam.track.corrLength;
init=fsmParam.track.init;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CHANGE TO WORK PATH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(userPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% TRACK SPECKLES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate grid size for enhanced tracking if needed

gridSize=0;
if enhanced==1
    nP=min(imgSize/11);
    if nP<5
        step=2*round(5/nP)+1;
        gridSize=round([step step]);
    else
        gridSize=[11 11];
    end
    
    % If no grid is needed, gridSize has to be set to 0
    if grid==0
        gridSize=0;
    end
    
end

% Initialize void coupling matrix
emptyM=zeros(1,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                             %
% POSSIBLE COMBINATIONS TRACK-INITIALIZER                                                                                     %
%                                                                                                                             %
%    init = 0: no tracker initialization                                                                                      %
%    init = 1: initialization through imKymoAnalysis - the user must have run imKymoAnalysis before SpeckTackle               %
%        TRACKER = 1: fsmTrackTrackerBMTNNMain takes care of checking that the files exists                                   %
%        TRACKER = 2: LAP does not support correlation - this combination should not be possible                              %
%    init = 2: initialization through TFT - this can be run here                                                              %
%       TRACKER = 1: tft will be run independently of whether initialization files already exist or not                       %
%       TRACKER = 2: LAP will call TFT internally                                                                             %
%                                                                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%  %                                                                                                                        % %
%  %                            THE SELECTED TRACKER IS NEURAL NETWORK TRACKER OR 'NO TRACKER'                              % % 
%  %                                                                                                                        % % 
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if TRACKER==1 | TRACKER ==3

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                                                                                        %
    %                                   THE SELECTED INITIALIZER IS 'Correlation'                                            % 
    %                                                                                                                        % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if init==1
    
        % Nothing to do here. 
        % The neural network tracker will make sure that the init files will be loaded correctly.
        % If 'no tracker', one can quit here.
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                                                                                        %
    %                                   THE SELECTED INITIALIZER IS 'TFT'                                                    % 
    %                                                                                                                        % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if init==2 % TFT has to be run BEFORE the neural network tracker
        
        % Number of images to be used for TFT
        lastImage=n-2;
        
        % Initializing progress bar
        h = waitbar(0,'Running TFT initializer...');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                               %
        %     PREPARE SPECKLE LISTS     %
        %                               %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Prepare three images and then track them with TFT
        for counter1=1:lastImage

            % Current index
            currentIndex=counter1+firstIndex-1;
            
            if counter1==1   % Only in the first round we need to process three files
            
                % Load speckle information for the FIRST image from cands###.mat
                indxStr=sprintf(strg,currentIndex);
                eval(['load ',userPath,filesep,'cands',filesep,'cands',indxStr,'.mat;']);  % Load cands for the current image
%                 % load from cands###_spa.mat if subpixel accuracy is
%                 % determined
%                 if fsmParam.prep.subpixel==1
%                     eval(['load ',userPath,filesep,'cands',filesep,'cands',indxStr,'_spa.mat;']);
%                     cands=candsSP;
%                 end
                
                
                % Create speckle map
                [I,img]=fsmTrackFillSpeckleList(cands,imgSize);
                if isempty(I)
                    error('Aborting tracking: no particles found in current image.');
                end
                clear cands;
                
                % Load speckle information for the SECOND image from cands###.mat
                indxStr=sprintf(strg,currentIndex+1);
                eval(['load ',userPath,filesep,'cands',filesep,'cands',indxStr,'.mat;']);  % Load cands for the current image
%                 % load from cands###_spa.mat if subpixel accuracy is
%                 % determined
%                 if fsmParam.prep.subpixel==1
%                     eval(['load ',userPath,filesep,'cands',filesep,'cands',indxStr,'_spa.mat;']);
%                     cands=candsSP;
%                 end
                
                % Create speckle map
                [J,img2]=fsmTrackFillSpeckleList(cands,imgSize);
                if isempty(J)
                    error('Aborting tracking: no particles found in current image.');
                end
                clear cands;
                
                % Load speckle information for the THIRD image from cands###.mat
                indxStr=sprintf(strg,currentIndex+2);
                eval(['load ',userPath,filesep,'cands',filesep,'cands',indxStr,'.mat;']);  % Load cands for the current image
%                 % load from cands###_spa.mat if subpixel accuracy is
%                 % determined
%                 if fsmParam.prep.subpixel==1
%                     eval(['load ',userPath,filesep,'cands',filesep,'cands',indxStr,'_spa.mat;']);
%                     cands=candsSP;
%                 end
                
                % Create speckle map
                [K,img3]=fsmTrackFillSpeckleList(cands,imgSize);
                if isempty(K)
                    error('Aborting tracking: no particles found in current image.');
                end
                clear cands;
                
            else % Only the third image needs to be processed
                
                img=img2; clear img2; % The first image is the second from the previous run
                I=J;      clear J;
                
                img2=img3; clear img3; % The second image is the third from the previous run    
                J=K;      clear K;
                    
                % Load only speckle information for the THIRD image from cands###.mat
                indxStr=sprintf(strg,currentIndex+2);
                eval(['load ',userPath,filesep,'cands',filesep,'cands',indxStr,'.mat;']);  % Load cands for the current image
                % load from cands###_spa.mat if subpixel accuracy is
                % determined
                if fsmParam.prep.subpixel==1
                    eval(['load ',userPath,filesep,'cands',filesep,'cands',indxStr,'_spa.mat;']);
                    cands=candsSP;
                end
                
                % Create speckle map
                [K,img3]=fsmTrackFillSpeckleList(cands,imgSize);
                clear cands;
            end
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                               %
            %           RUN TFT             %
            %                               %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            links = tft(I,J,K,threshold); 
            aux1_indxStr = sprintf(strg,currentIndex);
            savefile = [userPath,filesep,'links',filesep,'links',aux1_indxStr,'.mat'];
            save(savefile,'links')
            
            aux1_links = links(:,1:4);
            clear links
            if counter1 > 1
                aux2_indxStr=sprintf(strg,currentIndex-1);
                aux_path = [userPath,filesep,'links',filesep,'links',aux2_indxStr,'.mat'];
                load(aux_path)
                aux2_links = links(:,3:6);
                clear links
                aux_links = [aux1_links;aux2_links];
            else
                aux_links = aux1_links;
            end          
            flow = vectorFieldInterp(aux_links,aux_links(:,1:2),d0,[]);
            savefile = [userPath,filesep,'flow',filesep,'flow',aux1_indxStr,'.mat'];
            save(savefile,'flow')
            
            % Update wait bar
            waitbar(counter1/lastImage,h);
    
        end
        
        % Close waitbar
        close(h);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                               %
    %           NOW TRACK           %
    %                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    if TRACKER==1 % If needed, run the Neural Network tracker
        
        % Number of images to be used for the Nural Network tracker
        lastImage=n-1;
        
        % Initializing progress bar
        h = waitbar(0,'Tracking speckles...');
        
        for counter1=1:lastImage
            
            % Current index
            currentIndex=counter1+firstIndex-1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                                     % 
            % PREPARE THE IMAGES FOR TRACKING (FROM cands###.mat) %
            %                                                     % 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if counter1==1   % Only in the first round we need to process two files
                
                % Load speckle information for the FIRST image from cands###.mat
                indxStr=sprintf(strg,currentIndex);
                eval(['load ',userPath,filesep,'cands',filesep,'cands',indxStr,'.mat;']);  % Load cands for the current image
                % load from cands###_spa.mat if subpixel accuracy is
                % determined
                if fsmParam.prep.subpixel==1
                    eval(['load ',userPath,filesep,'cands',filesep,'cands',indxStr,'_spa.mat;']);
                    cands=candsSP;
                end
                
                % Create speckle map
                [I,img]=fsmTrackFillSpeckleList(cands,imgSize);
                if isempty(I)
                    error('Aborting tracking: no particles found in current image.');
                end
                clear cands;
                
                % Load speckle information for the SECOND image from cands###.mat
                indxStr=sprintf(strg,currentIndex+1);
                eval(['load ',userPath,filesep,'cands',filesep,'cands',indxStr,'.mat;']);  % Load cands for the current image
                % load from cands###_spa.mat if subpixel accuracy is
                % determined
                if fsmParam.prep.subpixel==1
                    eval(['load ',userPath,filesep,'cands',filesep,'cands',indxStr,'_spa.mat;']);
                    cands=candsSP;
                end
                
                % Create speckle map
                [J,img2]=fsmTrackFillSpeckleList(cands,imgSize);
                if isempty(J)
                    error('Aborting tracking: no particles found in current image.');
                end
                clear cands;
                
            else % Only the SECOND image needs to be processed
                
                img=img2; clear img2; % The first image is the second from the previous run
                I=J;      clear J;
                
                % Load only speckle information for the third image from cands###.mat
                indxStr=sprintf(strg,currentIndex+1);
                eval(['load ',userPath,filesep,'cands',filesep,'cands',indxStr,'.mat;']);  % Load cands for the current image
                % load from cands###_spa.mat if subpixel accuracy is
                % determined
                if fsmParam.prep.subpixel==1
                    eval(['load ',userPath,filesep,'cands',filesep,'cands',indxStr,'_spa.mat;']);
                    cands=candsSP;
                end
                
                % Create speckle map
                [J,img2]=fsmTrackFillSpeckleList(cands,imgSize);
                if isempty(J)
                    error('Aborting tracking: no particles found in current image.');
                end
                clear cands;
            end
            
            % Track with the Neural Network tracker
            tmp=fsmTrackTrackerBMTNNMain(I,J,threshold,influence,fsmParam,counter1,gridSize);
            
            % Check for unsupported error
            if isempty(tmp)
                error('Something went wrong during tracking...');
            end
            
            % Store speckle position information
            M(1:size(emptyM,1),1:size(emptyM,2),counter1,1)=emptyM;
            M(1:size(tmp,1),1:size(tmp,2),counter1,1)=tmp;
            
            % Update wait bar
            waitbar(counter1/lastImage,h);
            
        end
        
        %Delete any saved temporary 'tackFlow'.
        initPath=fsmParam.track.initPath;
        tackAvgFlowField = [initPath filesep 'tackFlow.mat'];
        if exist(tackAvgFlowField,'file') == 2
            delete(tackAvgFlowField);
        end
        
        % Close waitbar
        close(h);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                                                                                        %
    %                               THE SELECTED INITIALIZER IS NEITHER 'CORRELATION' NOR 'TFT'                              % 
    %                                                                                                                        % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if init~=0 & init~=1 & init~=2
        
        % Inform the user that this is a bug and quit
        error('The selected initializer does not exist. This is a bug. Please report it.');
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%  %                                                                                                                        % %
%  %                                   THE SELECTED TRACKER IS THE LINEAR ASSIGNMENT TRACKER                                % % 
%  %                                                                                                                        % % 
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if TRACKER==2 % Linear assignment tracker

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                                                                                        %
    %                                   THE SELECTED INITIALIZER IS 'Correlation'                                            % 
    %                                                                                                                        % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if init==1
    
        error('The LAP tracker cannot be initialized by ''Correlation''. This is a bug. Please report it.');
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                                                                                        %
    %                                   THE SELECTED INITIALIZER IS 'TFT'                                                    % 
    %                                                                                                                        % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if init==2
        
        % LAP can call TFT internally
        tempFileName = 'pointData.txt';
        tempFileName = [userPath, filesep, 'pointFiles', filesep, tempFileName];
        
        indxStr=sprintf(strg,fsmParam.specific.firstIndex);
        candsFileList=getFileStackNames([userPath,filesep,'cands',filesep,'cands',indxStr,'.mat']);
        % load from cands###_spa.mat if subpixel accuracy is
        % determined
        if fsmParam.prep.subpixel==1
            candsFileList=getFileStackNames([userPath,filesep,'cands',filesep,'cands',indxStr,'_spa.mat']);
        end
        LAP_pointDataFileGeneration(tempFileName, candsFileList)
        
        tempFileName = 'intensityProfile.txt';
        tempFileName = [userPath, filesep, 'pointFiles', filesep, tempFileName];
        LAP_intensityProfileGeneration(tempFileName, candsFileList);
        
        tempFileName = 'LAPconfig.txt';
        tempFileName = [userPath, filesep, 'pointFiles', filesep, tempFileName];
        
        % By passing smParam.track.init as input parameter, tft will be turned on or off depending on its value
        LAP_configurationFileGeneration(tempFileName, userPath, fsmParam.specific.imgSize(1), fsmParam.specific.imgSize(2), lastImage+1, threshold, fsmParam.track.init);
        
        LAPpath = which('LAPTrack65');
        indxSep=findstr(LAPpath,filesep);
        LAPpath=LAPpath(1:indxSep(end));
        
        LAP_batchFileGeneration(tempFileName, userPath, LAPpath)
        
        msg=['Please manually start the batch file ''lapbatch.bat'' in ',userPath,filesep,'pointFiles, and click on OK when the tracking is done.'];
        uiwait(msgbox(msg));
        
        tempFileName = 'LAPconfig_track.txt';
        tempFileName = [userPath, filesep, 'pointFiles', filesep, tempFileName];
        track = trackReadAll(tempFileName)
        
        % convert to M
        M = LAPtrack2M(track, lastImage+1);
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%  %                                                                                                                        % %
%  %                                          IF 'NO TRACKER' WAS SELECTED, RETURN HERE                                     % % 
%  %                                                                                                                        % % 
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if TRACKER==3
    status=1;
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%  %                                                                                                                        % %
%  %                                    THE SELECTED TRACKER IS NEITHER OF THE EXISTING OPTION                              % % 
%  %                                                                                                                        % % 
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if TRACKER~=1 & TRACKER~=2 & TRACKER~=3
    
    % Inform the user that this is a bug and quit
    error('The selected tracker does not exist. This is a bug. Please report it.');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%  %                                                                                                                        % %
%  %                                                  GAP CLOSER AND LINKER                                                 % % 
%  %                                                                                                                        % % 
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               %
%           GAP CLOSER          %
%                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Correct dimensions
cM=M>0; [i j k]=find(cM); M=M(1:max(i),:,:,:);

if size(M,3)>1 % In case only two frames have been tracked, no gaps can be closed
    
    if TRACKER~=2 % No gap closer for the Linear Assignment tracker
        
        % Not running the gapCloser, no gapList files will be saved. This will allow next time to create them.
        if enhanced==1
            % Call external function to perform enhanced gap closing
            [M,closedGap]=fsmTrackEnhancedGapCloser(M,strg,threshold,userPath,firstIndex);
        else
            % Close gaps in M
            closedGap=0; [M,closedGap]=fsmTrackGapCloser(M,threshold,strg,userPath,firstIndex);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                               %
    %            LINKER             %
    %                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Build speckle sequence
    [MPM,M]=fsmTrackLinker(M);
    
else
    
    % For two frames, MPM is equal to M
    MPM=M;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   SAVE MAGIC POSITION MATRIX
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

infoH=fsmGuiInfo; ch=get(infoH,'Children');
textString=['Saving Magic Postion Matrix to ',userPath,'\mpm.mat'];
set(ch(3),'String','TRACKER MODULE');
set(ch(2),'String',textString);
set(ch(1),'String','Prease wait...');
save mpm.mat MPM M;
set(ch(1),'String','Done!');
pause(1);
close(infoH);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SETTING MODULE STATUS TO 1 AND RETURNING
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the status to 1 to mean that the module successfully finished
status=1;
