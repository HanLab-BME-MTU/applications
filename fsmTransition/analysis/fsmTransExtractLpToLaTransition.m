function [fPositionsScore,fPositionsVel,fPositionsTrans,dpRatio,scoreProfiles,speedProfiles]=fsmTransExtractLpToLaTransition(tackProjDir,toggleKinetic,toggleSpeed,toggleDpRatio,profiles,dist,n,DEBUG,medianFiltFrames,minDist)
% fsmTransExtractLpToLaTransition extracts the lamellipodium-to-lamellum transition by analyzing kinetic and speed profiles at the leading edge
%
% SYNOPSIS   [fPositionsScore,fPositionsVel,fPositionsTrans,dpRatio,scoreProfiles,speedProfiles]=
%                  fsmTransExtractLpToLaTransition(tackProjDir,toggleKinetic,toggleSpeed,toggleDpRatio,profiles,dist,n,DEBUG,medianFiltFrames,minDist)
%
% INPUT      tackProjDir        : tack project subdirectory
%            toggleKinetic      : [ 0 | 1 ] turns on|off the calculation of kinetic profiles as a criterion for transition extraction
%            toggleSpeed        : [ 0 | 1 ] turns on|off the calculation of kinematic profiles as a criterion for transition extraction
%            toggleDpRatio      : calculates the ratio of the integrals of depoly and poly rates from the leading edge to the transition
%            profiles           : profiles structure as returned by fsmTransReadPrPanelDatFiles
%            dist               : defines the maximum distance (support) in pixels from a given position along the profile for kinetic
%                                 and kinematic measurements (scores and vectors, respectively) to be taken into account for averaging 
%                                 onto the profile
%            n                  : numer of frames for temporal integration of kinetic score / profiles
%                                 set n=[] for default (all frames)
%                                 set n=-1 for all frames ( = size(M,3))
%            DEBUG              : draws (and saves) debug figures (CAREFUL: this can create A LOT of figures)
%            medianFiltFrames   : order of the media filter (to smoothen the obtained transition line)
%                                 (optional, default = 5 - if you want to pass a value for minDist, you MUST pass a value for medianFiltFrames)
%            minDist            : Minimum distance from the leading edge (in pixels) for a transition to be valid
%                                 (optional, default = 9)
%
% OUTPUT     positionsScore     : transition coordinates for all frames calculated from kinetic criterion only

%            fPositionsScore    : transition coordinates for all frames calculated from kinetic criterion only
%                                 and median filtered with order = medianFiltFrames (see INPUT)
%                                 size = [ number of profiles x 2 x frames]
%            fPositionsVel      : transition coordinates for all frames calculated from kinematic criterion only
%                                 and median filtered with order = medianFiltFrames (see INPUT)
%                                 size = [ number of profiles x 2 x frames]
%            fPositionsTrans    : transition coordinates for all frames calculated from averaged kinetic and kinematic criteria
%                                 and median filtered with order = medianFiltFrames (see INPUT)
%                                 size = [ number of profiles x 2 x frames]
%            dpRatio            : depoly/poly ratio along the transition; size = [number of profiles x frames ]
%            scoreProfile       : all calculated kinetic profiles; size = [number of profiles x length of profile (pixels) x frames ]
%            speedProfile       : all calculated kinematic profiles; size = [number of profiles x length of profile (pixels) x frames ]
%
% REMARK     THIS FUNCTION ASSUMES THAT THE DIFFERENCE IN FRAMES DUE TO THE TIME-AVERAGING OF THE TRANSITION HAS ALREADY BEEN TAKEN INTO
%            ACCOUNT WHEN THE STRUCTURE OF SPLINES HAS BEEN GENERATED (this is the case if the function is called through the user 
%            interface FSMTRANSITION).
%
% DEPENDENCES   fsmTransExtractLpToLaTransition uses { }
%               fsmTransExtractLpToLaTransition is used by { fsmTransition }
%
% Aaron Ponti, September 6th, 2004

% Debug info flag
if DEBUG~=0 & DEBUG~=1
    DEBUG=0;
end

% Global variable
global uFirst uLast

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS, OUTPUTS AND PARAMETERS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default value for the median filter order medianFiltFrames and minDist
if nargin==8
    medianFiltFrames=5; % Order of median filter
    minDist=9;
end

if nargin==9
    minDist=9;
end

% Default value for the number of time points n if n is empty
if isempty(n)
    n=0;
end

% Check that at least one of the criteria is turned on
if toggleKinetic==0 & toggleSpeed==0
    error('At least one criterion should be on');
end

% Initialize outputs
positionsScore=[];positionsVel=[];positionsTrans=[];fPositionsTrans=[];dpRatio=[];scoreProfiles=[];speedProfiles=[];

% Initialize vectors for plots
Pvel=[];Pkin=[];velTransition=[];scoreTransition=[];velProfile=[];scoreProfile=[];
gScore=[];g2Score=[];gVel=[];g2Vel=[]; 

% Some hardwired parameters
d0_score=5;
d0_vel=5; % Turned off low-pass filtering of kinetic maps
sampling=1;


% Turn off the 'Plot empty' warning
warning off 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LOOK FOR AND CHECK fsmParam.mat
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fsmParam,status]=fsmPostGetFsmParam(tackProjDir);
if status==0
    return
end
bitDepth=log2(fsmParam.main.normMax+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CHECK THAT THE TRACKING MODULE IS UP-TO-DATE AND LOAD mpm.mat
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(fsmParam.track,'uptodate')

    % Check that the tracking module is up-to-date
    if fsmParam.track.uptodate==0
        
        % The tracking module is not up-to-date. Inform the user and return    
        errordlg('The file ''mpm.mat'' is not up-to-date. Please re-run the Tracking module in SpeckTackle.','Error','modal');
        return

    else
        uptodate=1;    
    end
    
else
    
    % Old version of fsmParam.mat. Inform the user that he/she has to make sure that everything is up-to-date
    uiwait(msgbox('Since ''fsmParam.mat'' has been created by an old version of SpeckTackle, I cannot make sure that ''mpm.mat'' is up-to-date. Continue at your own risk.','Warning'));
    uptodate=-1;
    
end

% Load mpm.mat
if exist([tackProjDir,filesep,'mpm.mat'],'file');
    load([tackProjDir,filesep,'mpm.mat']);
    clear MPM;
else
    errordlg('Could not find ''mpm.mat''.','Error','modal');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CHECK THAT THE KINETIC MODULE IS UP-TO-DATE AND GET THE LIST OF kinScore###.mat FILES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(fsmParam.kin,'uptodate')

    % Check that the tracking module is up-to-date
    if fsmParam.kin.uptodate==0
        
        % The tracking module is not up-to-date. Inform the user and return    
        errordlg('The ''kinScore###.mat'' files are not up-to-date. Please re-run the Kinetic module in SpeckTackle.','Error','modal');
        return

    else
        uptodate=1;    
    end
    
else
    
    % Old version of fsmParam.mat. Inform the user that he/she has to make sure that everything is up-to-date
    uiwait(msgbox('Since ''fsmParam.mat'' has been created by an old version of SpeckTackle, I cannot make sure that the ''kinScore###.mat'' files are up-to-date. Continue at your own risk.','Warning'));
    uptodate=-1;
    
end

% Get list of kinScore files
[kinScoreFileList,success]=fsmPostGetSubProjFileList([tackProjDir,filesep,'kinScore'],'kinScore');
len=length(kinScoreFileList);

% First and last kinScore indices
[tmpPath,tmpBody,firstKinScoreIndex,tmpExt]=getFilenameBody(char(kinScoreFileList(1)));
[tmpPath,tmpBody,lastKinScoreIndex,tmpExt]=getFilenameBody(char(kinScoreFileList(end)));

% Format string for numerical index
strg=fsmParam.specific.formString;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GET THE LIST OF IMAGE FILE NAMES AND FIRST/LAST INDICES FROM fsmParam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[imageFileList,firstImageIndex,lastImageIndex,status]=fsmPostGetImageListFromFsmParam(fsmParam);
if status==0
    return
end
imgSize=fsmParam.specific.imgSize;

% Make sure that images, kinScores, profiles, and M all match
if str2num(firstKinScoreIndex)~=firstImageIndex | str2num(lastKinScoreIndex)>lastImageIndex | length(kinScoreFileList)~=(size(M,3)+1) | length(kinScoreFileList)~=size(profiles,3)
    error('Some of images, kinScores, profiles, and M do not match. This will be extended in the future to be more flexible.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% START
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First and last frames
first=str2num(firstKinScoreIndex);
last=str2num(lastKinScoreIndex);

% Select range of frames
guiH=fsmTrackSelectFramesGUI;
set(findall(0,'Tag','pushOkay'),'UserData',n-1); % Minimum range 
title='Select images to be processed:';
set(findall(0,'Tag','editFirstFrame'),'String',num2str(first));
set(findall(0,'Tag','editLastFrame'),'String',num2str(last));
set(findall(0,'Tag','SelectFramesGUI'),'Name',title);
sSteps=[1/(last-first) 1/(last-first)];
set(findall(0,'Tag','sliderFirstFrame'),'SliderStep',sSteps,'Max',last,'Min',first,'Value',first);
set(findall(0,'Tag','sliderLastFrame'),'SliderStep',sSteps,'Max',last,'Min',first,'Value',last);
waitfor(guiH); % The function waits for the dialog to close (and to return values for uFirst and uLast)

if uFirst==-1
    return % The user closed the dialog
end

% Keep only the file names in the user-selected range
kinScoreFileList=kinScoreFileList(uFirst:uLast);
imageFileList=imageFileList(uFirst:uLast,:);

% Keep only the profiles for the frames in the user-selected range
profiles=profiles(:,:,uFirst:uLast);

% Keep only tracks in the user-selected range
M=M(:,:,uFirst:uLast-1);

% Update first and last
first=1;
last=length(kinScoreFileList);

if n==0
    n=length(kinScoreFileList);
    last=first;
else
    last=length(kinScoreFileList)-n+1;
end

% Initialize vector of figure handles
h=zeros(1,last-first+1);

% Allocate memory for some of the outputs
positionsScore=zeros(size(profiles,1),2,last-first+1);
positionsVel=positionsScore;
positionsTrans=positionsScore;
fPositionsTrans=positionsScore;

% Initialize waitbar
hWait=waitbar(0,'Please wait...');

for j=first:last % Go through all images
    
    % Current index
    current=first+j-1; % With current implementation, current corresponds to j

    % Load current image
    eval(['img=imreadnd2(''',(char(imageFileList(current,:))),''',1,2^bitDepth-1);']);

    if DEBUG==1
        h(current)=figure;imshow(img,[]);
        hold on
    end
    
    for i=1:size(profiles,1) % Go through all profiles
       
        % Take one profile out of profiles
        lineDescr=[profiles(i,1,current) profiles(i,3,current);profiles(i,2,current) profiles(i,4,current)];
        
        if toggleKinetic==1
            
            % Kinetic profile
            [scoreProfile,polyProfile,depolyProfile,Pkin,posScores,whatever]=fsmKineticProfile(char(kinScoreFileList(current)),imgSize,n,lineDescr,dist,d0_score,sampling);

            % Calculate first and second derivative of the kinetic map
            gScore=gradient(scoreProfile);
            g2Score=gradient(gScore);
            
            % Find a zero in the first derivative of the kinetic profile
            %    The abs of the diff of the sign of the gradient is 2 where the 
            %    sign goes from positive to negative
            gZerosScore=[[abs(diff(sign(gScore)))==2]; 0]; % To maintain vector length
            
            % Look for a maximum in the kinetic activity (not the first one)
            %   (second derivative < 0)
            scoreTransition=find((gZerosScore.*g2Score)<0); % Maximum 
            scoreTransition=scoreTransition(find(scoreTransition>minDist));
            if isempty(scoreTransition)
                scoreTransition=minDist;
                disp('No minima found in kinetic activity profile');
            end
            scoreTransition=scoreTransition(1);
            
            % Calculate depoly/poly ratio from leading edge to transition (using the kinetic transition)
            if toggleDpRatio==1
                currentDpRatio=sum(abs(depolyProfile(1:scoreTransition)))/(sum(polyProfile(1:scoreTransition))+eps); % To prevent division by zero
            end
            
        end
        
        if toggleSpeed==1
            
            % Speed profile (use n-1 since each matrix in M represents a frame PAIR)
            [velProfile,Pvel,Mi,list]=fsmVelocityProfile(M,n-1,lineDescr,dist,d0_vel,sampling);
            
            
            % Calculate first and second derivative of the velocity
            gVel=gradient(velProfile);
            g2Vel=gradient(gVel);
            
            % Find a zero in the first derivative of the speed profile
            %    The abs of the diff of the sign of the gradient is 2 where the 
            %    sign goes from positive to negative
            gZerosVel=[[abs(diff(sign(gVel)))==2]; 0]; % To maintain vector length
            
            % Look for the first minimum (second derivative > 0)
            velTransition=find((gZerosVel.*g2Vel)>0); 
            velTransition=velTransition(find(velTransition>minDist));
            if isempty(velTransition)
                velTransition=minDist;
                disp('No minima found in speed profile');
            end
            velTransition=velTransition(1);
            
        end
        
        % Plot and store found transitions
        if DEBUG==1
            figure(h(current));
            plot([lineDescr(2) lineDescr(4)],[lineDescr(1) lineDescr(3)],'r-');
        end
        
        if toggleKinetic==1
            
            if DEBUG==1
                plot(Pkin(scoreTransition,2),Pkin(scoreTransition,1),'y.');
            end
            
            % Store value
            positionsScore(i,1:2,current)=Pkin(scoreTransition,1:2);
            
            % Store score profile - allocate memory if first time
            if i==1 & j==first % first profile of first image
                scoreProfiles=zeros(size(profiles,1),length(scoreProfile),last-first+1);
            end
            scoreProfiles(i,1:length(scoreProfile),current)=scoreProfile';
            
            if toggleDpRatio==1
                dpRatio(i,current)=currentDpRatio;
            end
        end
        
        if toggleSpeed==1
            
            if DEBUG==1
                plot(Pvel(velTransition,2),Pvel(velTransition,1),'g.');
            end
            
            % Store value
            positionsVel(i,1:2,current)=Pvel(velTransition,1:2);
            
            % Store speed profile - allocate memory if first time
            if i==1 & j==first % first profile of first image
                speedProfiles=zeros(size(profiles,1),length(velProfile),last-first+1);
            end
            speedProfiles(i,1:length(velProfile),current)=velProfile';

        end
        
        % The transition is set to be at half-way between the two transitions (if both exist)
        if toggleKinetic==1 & toggleSpeed==1
            positionsTrans(i,1:2,current)=mean([Pkin(scoreTransition,1:2);Pvel(velTransition,1:2)],1);
        elseif toggleKinetic==1 & toggleSpeed==0
            positionsTrans(i,1:2,current)=Pkin(scoreTransition,1:2);
        elseif toggleKinetic==0 & toggleSpeed==1
            positionsTrans(i,1:2,current)=Pvel(velTransition,1:2);
        else
            error('At least one criterion should be on');
        end

        
        if DEBUG==1
            
            % Display profiles and derivatives
            figure,
            % Score
            subplot(6,1,1), plot([1:length(Pkin)],scoreProfile,'r-');legend('k');
            hold on; plot(scoreTransition,scoreProfile(scoreTransition),'r*');
            % Velocity
            subplot(6,1,2), plot([1:length(Pvel)],velProfile,'k-');legend('v');
            hold on; plot(velTransition,velProfile(velTransition),'k*');
            % Score - first derivative
            subplot(6,1,3), plot([1:length(Pkin)],gScore,'b-');legend('dk');
            hold on; plot(scoreTransition,gScore(scoreTransition),'b*');
            plot([1:length(Pkin)],max(gScore).*sign(gScore),'b-')
            % Score - second derivative
            subplot(6,1,4), plot([1:length(Pkin)],g2Score,'c-');legend('d2k');
            hold on; plot(scoreTransition,g2Score(scoreTransition),'c*');
            plot([1:length(Pkin)],max(g2Score).*sign(g2Score),'c-')
            % Velocity - first derivative
            subplot(6,1,5), plot([1:length(Pvel)],gVel,'g-');legend('dv');
            hold on; plot(velTransition,gVel(velTransition),'g*');
            plot([1:length(Pvel)],max(gVel).*sign(gVel),'g-')
            % Velocity - second derivative
            subplot(6,1,6), plot([1:length(Pvel)],g2Vel,'m-');legend('d2v');
            hold on; plot(velTransition,g2Vel(velTransition),'m*');
            plot([1:length(Pvel)],max(g2Vel).*sign(g2Vel),'m-')
            
        end
            
        % Update waitbar
        waitbar(((j-1)*size(profiles,1)+i)/((last-first+1)*size(profiles,1)),hWait)
        
    end
    
    % Plot transition if needed
    if DEBUG==1
        figure(h(current));
        plot(positionsTrans(:,2,current),positionsTrans(:,1,current),'y-.');
    end
    
    % Median filter the positions calculated with the combined critera ...
    fPositionsTrans(:,1,current)=medfilt1(positionsTrans(:,1,current),medianFiltFrames);
    fPositionsTrans(:,2,current)=medfilt1(positionsTrans(:,2,current),medianFiltFrames);
    
    % ... with the kinetic criterion
    fPositionsScore(:,1,current)=medfilt1(positionsScore(:,1,current),medianFiltFrames);
    fPositionsScore(:,2,current)=medfilt1(positionsScore(:,2,current),medianFiltFrames);
    
    % ... and with the speed criterion
    fPositionsVel(:,1,current)=medfilt1(positionsVel(:,1,current),medianFiltFrames);
    fPositionsVel(:,2,current)=medfilt1(positionsVel(:,2,current),medianFiltFrames);
    
    % Plot filtered transition if needed
    if DEBUG==1
        plot(fPositionsTrans(:,2,current),fPositionsTrans(:,1,current),'r-','LineWidth',1);
    end

    % Save figure if needed
    if DEBUG==1
        
        % Save current figure
        indxStr=sprintf(strg,current);
        name=['transition',indxStr,'.fig'];
        saveas(h(current),name,'fig');
    end
                

end

% Close waitbar
close(hWait);

% Turn the warnings back on
warning on

