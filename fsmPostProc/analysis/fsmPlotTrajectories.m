function fsmPlotTrajectories(projDir,method,colorCode)
% fsmPlotTrajectories plots speckle trajectories
%
% SYNOPSIS      fsmPlotTrajectories(projDir,method,colorCode)
%
% REMARK        Run fsmPlotTrajectories through the graphical user interface fsmCenter 
%               ----------------------------------------------------------------------
%
% INPUT         projDir  : string pointing to the project directory (where the fsmParam.mat file is 
%                             located). Pass projDir=[] to manually select fsmParam.mat via a dialog.
%               method   : 1 - All the matches in a selected frame range are displayed in 
%                              one plot. Trajectories longer than the range will be cut.
%                          2 - All the trajectories originitating within the selected frame 
%                              range are displayed in one plot. Trajectories longer than the 
%                              range will be displayed in their entirety.
%               colorCode: [ 0 | 1 ] Turns off | on the color-coding of the trajectories. Optional, default = 0;
%                          Colors represent time. For method:
%                                  1 - Trajectories are multicolor. Color of a segment along the
%                                      trajectory corresponds to the frame.
%                                  2 - All trajectories starting at a given frame will be colored the same.
%
% The function will require the user to set the range of interest via a dialog.
%
% OUTPUT        None  
%
% DEPENDENCES   fsmPlotTrajectories uses { }
%               fsmPlotTrajectories is used by { fsmCenter }
%
% Aaron Ponti, January 22th, 2004

if nargin==1
    colorCode=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LOOK FOR AND CHECK fsmParam.mat
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fsmParam,status]=fsmPostGetFsmParam(projDir);
if status==0
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GET THE LIST OF IMAGE FILE NAMES AND FIRST/LAST INDICES FROM fsmParam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[imageFileList,imageFirstIndex,imageLastIndex,status]=fsmPostGetImageListFromFsmParam(fsmParam);
if status==0
    return
end

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
if exist([projDir,filesep,'mpm.mat'],'file');
    load([projDir,filesep,'mpm.mat']);
    %clear MPM;
else
    errordlg('Could not find ''mpm.mat''.','Error','modal');
    return
end

% Select first and last frame
last=size(M,3)+1;
if method==1
    first=1;
    minDistFrame=1; 
    clear MPM; % Not needed
else
    first=2;
    minDistFrame=0;
    clear M; % Not needed
end

% Select frames
[uFirst,uLast]=fsmTrackSelectFramesGUI(first,last,minDistFrame,'Plot trajectories within frames:');

% make sure that values are integers!!!
uFirst = int16(uFirst);
uLast  = int16(uLast);

% Check whether the dialog was closed (_CloseRequestFcn)
if uFirst==-1
    return % This will return an error (status=0)
end

if colorCode==1
    % Initialize colormap
    nc = uint16(uLast-uFirst+1);
    cmap=colormap(jet(122));
    close(gcf);
    %cmap=colormap(jet(156));close(gcf);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GET SOME IMAGE INFO
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load first image
img=double(imread(char(imageFileList(int16(uFirst),:))));

% Store image size
imgSize=size(img);

% Open figure and axes
figPlotH=figure;
img=(img-max(img(:)))/(min(img(:))-max(img(:)));
imshow(img,[]); 
hold on
set(figPlotH,'NumberTitle','off');
switch method
    case 1, strg=['Trajectories contained in frame range ',num2str(uFirst),' - ',num2str(uLast)];
    case 2, strg=['Trajectories originating in frame range ',num2str(uFirst),' - ',num2str(uLast)];
    otherwise, error('Wrong method selected');
end
set(figPlotH,'Name',strg);

% Initialize index for colormap
cmapIndex=0;

% Draw trajectories
if method==1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % METHOD 1: Plot all matches per frames over the range uFirst:uLast
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Go through frames
    for i=uFirst:uLast-1
        
        % Update color map index
        cmapIndex=cmapIndex+1;
        
        Mo=M(:,:,i);
        Mv=Mo(find(Mo(:,1)~=0 & Mo(:,3)~=0),:);
        
        % Check for MATLAB version - quiver 7.0 is no longer compatible
        v=ver('MATLAB');
        pointPos=findstr(v.Version,'.');
        if ~isempty(pointPos)
            v.Version=v.Version(1:pointPos(1)+1);
        end
        
        %Plot arrows
        if str2num(v.Version)<7
            h=quiver(Mv(:,2),Mv(:,1),Mv(:,4)-Mv(:,2),Mv(:,3)-Mv(:,1),0);
        else
            h=quiver('v6',Mv(:,2),Mv(:,1),Mv(:,4)-Mv(:,2),Mv(:,3)-Mv(:,1),0);
        end
        if colorCode==1
            set(h(1),'Color',cmap(cmapIndex,:)); set(h(1),'LineWidth',1);
            set(h(2),'Color',cmap(cmapIndex,:)); set(h(1),'LineWidth',1);
        else
            set(h(1),'Color','red');
            set(h(2),'Color','red');
        end
        
    end 
    
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % METHOD 2: Plot entire trajectories ORIGINATING between frames uFirst and uLast (SLOW!)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i=uFirst:uLast
        
        % Update color map index
        cmapIndex=cmapIndex+1;
        
        % Extract all trajectories starting in the current frame
        MPM0=MPM(:,2*i-3:2*i);
        [y,x]=find(MPM0(:,1)==0 & MPM0(:,3)~=0); % Check [0 0] in the previous frame
        
        % Initialize enough space to store trajectories
        traj=zeros(1,length(y)*fix(0.5*size(MPM,2)));
        
        % Find end of trajectories
        currentTraj=MPM(y,2*i-1:end);


        maxLength=size(currentTraj,2);
        
        % Store all trajectories for current frame in a long array
        %    Trajectories are separated by NaN
        currentPos=1;
        for j=1:size(currentTraj,1)
            indx=find(currentTraj(j,:)==0);
            if isempty(indx)
                indx=maxLength+1;
            end
            if length(indx)>1
                indx=indx(1);
            end
            indx=indx-1;
            
            % Discard ghost speckles
            if indx>2
                traj(currentPos:currentPos+indx-1)=currentTraj(j,1:indx);
                currentPos=currentPos+indx;
                traj(currentPos:currentPos+1)=[NaN NaN]; % Division with the next trajectory
                currentPos=currentPos+2;
            else
                % 
            end
        end
        
        % Crop not-used entries
        traj=traj(1:currentPos-3);

        % Plot
        if colorCode==1
            plot(traj(2:2:end),traj(1:2:end-1),'Color',cmap(cmapIndex,:));
            % Mark beginning of trajectories
            bg=find(isnan(traj)); bg=bg(2:2:end)+1;
            plot(traj(bg+1),traj(bg),'.','Color',cmap(cmapIndex,:));
            
        else
            plot(traj(2:2:end),traj(1:2:end-1),'r');
            % Mark beginning of trajectories
            bg=find(isnan(traj)); bg=bg(2:2:end)+1;
            plot(traj(bg+1),traj(bg),'r.');
        end
          
    end
    
end

% Reformat axes
axis ij
axis equal
axis([1 imgSize(2) 1 imgSize(1)]);

% if uLast-uFirst>1 % We don't need a colorbar if only one frame pair is displayed
if (colorCode==1 & method==1 & uLast-uFirst>1) | (colorCode==1 & method==2 & uLast-uFirst>0)

    % Add colorbar
    figH=figure;
    cbar=fix(repmat(1:0.01:size(cmap,1),100,1));
    imshow(cbar,[]);
    colormap(cmap);
    set(figH,'NumberTitle','off');
    strg=['Selected interval: frame ',num2str(uFirst),' (blue) -> frame ',num2str(uLast),' (green)'];
    set(figH,'Name',strg);
end

% Add common tools menu
fsmCenterCB_addToolsMenu(figPlotH);

