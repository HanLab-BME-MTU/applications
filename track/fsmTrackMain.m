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

if enhanced==1
    nP=min(imgSize/11);
    if nP<5
        step=2*round(5/nP)+1;
        gridSize=round([step step]);
    else
        gridSize=[11 11];
    end

    % Calculate d0 for interpolation
    d0=3*max(gridSize);
    
    % If no grid is needed, gridSize has to be set to 0
    if grid==0
        gridSize=0;
    end
    
end

% Initialize void coupling matrix
emptyM=zeros(1,4);

% Initializing progress bar
h = waitbar(0,'Tracking speckles...');

% Setting last image to be evaluated (depends on the tracker)
if TRACKER==3
    lastImage=n-2;
else
    lastImage=n-1;
end

for counter1=1:lastImage

    % Current index
    currentIndex=counter1+firstIndex-1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % PREPARE THE IMAGES FOR TRACKING (FROM cands###.mat)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if counter1==1   % Only in the first round we need to process three files
                
        % Load speckle information for the first image from cands###.mat
        indxStr=sprintf(strg,currentIndex);
        %eval(strcat('load cands',filesep,'cands',indxStr,'.mat;')); % Load cands for the current image
        eval(['load ',userPath,filesep,'cands',filesep,'cands',indxStr,'.mat;']);  % Load cands for the current image
        
        % Create speckle map
        [I,img]=fsmTrackFillSpeckleList(cands,imgSize);
        if isempty(I)
            error('Empty list of coordinates from image 1.');
        end
        clear cands;
        
        % Load speckle information for the second image from cands###.mat
        indxStr=sprintf(strg,currentIndex+1);
        %eval(strcat('load cands',filesep,'cands',indxStr,'.mat;')); % Load cands for the current image
        eval(['load ',userPath,filesep,'cands',filesep,'cands',indxStr,'.mat;']);  % Load cands for the current image
       
        % Create speckle map
        [J,img2]=fsmTrackFillSpeckleList(cands,imgSize);
        if isempty(J)
            error('Empty list of coordinates from image 2.');
        end
        clear cands;
        
        if TRACKER==3
            % Load speckle information for the third image from cands###.mat
            indxStr=sprintf(strg,currentIndex+2);
            %eval(strcat('load cands',filesep,'cands',indxStr,'.mat;')); % Load cands for the current image
            eval(['load ',userPath,filesep,'cands',filesep,'cands',indxStr,'.mat;']);  % Load cands for the current image

            % Create speckle map
            [K,img3]=fsmTrackFillSpeckleList(cands,imgSize);
            if isempty(K)
                error('Empty list of coordinates from image 3.');
            end
            clear cands;
        end
           
    else % Only the third image needs to be processed
        
        img=img2; clear img2; % The first image is the second from the previous run
        I=J;      clear J;
        
        if TRACKER==3    
            img2=img3; clear img3; % The second image is the third from the previous run    
            J=K;      clear K;
            
            % Load only speckle information for the third image from cands###.mat
            indxStr=sprintf(strg,currentIndex+2);
            %eval(strcat('load cands',filesep,'cands',indxStr,'.mat;')); % Load cands for the current image
            eval(['load ',userPath,filesep,'cands',filesep,'cands',indxStr,'.mat;']);  % Load cands for the current image
            
            % Create speckle map
            [K,img3]=fsmTrackFillSpeckleList(cands,imgSize);
            clear cands;
        else
        
            % Load only speckle information for the third image from cands###.mat
            indxStr=sprintf(strg,currentIndex+1);
            %eval(strcat('load cands',filesep,'cands',indxStr,'.mat;')); % Load cands for the current image
            eval(['load ',userPath,filesep,'cands',filesep,'cands',indxStr,'.mat;']);  % Load cands for the current image

            % Create speckle map
            [J,img2]=fsmTrackFillSpeckleList(cands,imgSize);
            if isempty(J)
                error('Empty list of coordinates from image 3.');
            end
            clear cands;
        end
        
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % TRACK
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    switch TRACKER
    case 1
%         [tmp,n1,n2,nc1,nc2]=fsmTrackTrackerA(img,img2,img,img2);
        tmp=fsmTrackTrackerBMTNN(I,J,threshold,influence);
    case 2
        tmp=fsmTrackTrackerP(img,img2,threshold);    
    case 3
        tmp=fsmTrackTrackerPP(img,img2,img3,threshold);    
    otherwise
        error('Please select an EXISTING tracker.');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % ENHANCED TRACKING
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if enhanced==1
        
        % Call external function to perform enhanced tracking
        tmp=fsmTrackEnhancedTracker(tmp,I,J,strg,currentIndex,gridSize,d0,userPath,threshold,influence,size(img));
        
    end
    
    % Check for unsupported error
    if isempty(tmp)
        error('Something went wrong during tracking...');
    end
    
    % Store speckle position information
    M(1:size(emptyM,1),1:size(emptyM,2),counter1,1)=emptyM;
    M(1:size(tmp,1),1:size(tmp,2),counter1,1)=tmp;
    
    % Update wait bar
    waitbar(counter1/(n-1),h);

    % Force matlab to update the waitbar
    drawnow;
 
end

% Close waitbar
close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   GAP CLOSER
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Correct dimensions
cM=M>0; [i j k]=find(cM); M=M(1:max(i),:,:,:);

if size(M,3)>1 % In case only two frames have been tracked, no gaps can be closed
    % Not running the gapCloser, no gapList files will be saved. This will allow next time to create them.
    if enhanced==1
        % Call external function to perform enhanced gap closing
        [M,closedGap]=fsmTrackEnhancedGapCloser(M,strg,threshold,userPath,firstIndex);
    else
        % Close gaps in M
        closedGap=0; [M,closedGap]=fsmTrackGapCloser(M,threshold,strg,userPath,firstIndex);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   LINKER
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
