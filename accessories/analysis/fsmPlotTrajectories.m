function fsmPlotTrajectories(method)
% fsmPlotTrajectories plots speckle trajectories
%
% REMARK        Run fsmPlotTrajectories through the graphical user interface fsmCenter 
%               ------------------------------------------------------------------
%
% INPUT         method   : [1 | 2] Method 2 is SLOW!
%                          1 - All the matches in a selected frame range are displayed in 
%                              one plot. Trajectories longer than the range will be cut.
%                          2 - All the trajectories originitating within the selected frame 
%                              range are displayed in one plot. Trajectories longer than the 
%                              range will be displayed in their entirety.
%
% The function will require the user to set the range of interest via a dialog.
%
% OUTPUT        None  
%
% DEPENDENCES   fsmPlotTrajectories uses { }
%               fsmPlotTrajectories is used by { fsmCenter }
%
% Aaron Ponti, January 22th, 2004

global uFirst uLast

% Load image - this is needed to get image size
[fName,dirName] = uigetfile(...
    {'*.tif;*.tiff;*.jpg;*.jpeg','Image Files (*.tif,*.tiff,*.jpg,*.jpeg)';
    '*.tif','TIF files (*.tif)'
    '*.tiff','TIFF files (*.tiff)'
    '*.jpg;','JPG files (*.jpg)'
    '*.jpeg;','JPEG files (*.jpeg)'
    '*.*','All Files (*.*)'},...
    'Select starting image');
if(isa(fName,'char') & isa(dirName,'char'))
    img=nrm(imread([dirName,fName]),1);
    % Store image size
    imgSize=size(img);
else
    return 
end

% Load MPM.mat
[fName,dirName] = uigetfile(...
    {'MPM.mat;','MPM.mat';
    '*.*','All Files (*.*)'},...
    'Select MPM.mat');
if ~(isa(fName,'char') & isa(dirName,'char'))
    return 
end
load([dirName,fName]);
clear MPM;

% Select first and last frame
if method==1
    first=1;
else
    first=2;
end
last=size(M,3)-1;
guiH=fsmTrackSelectFramesGUI; ch=get(guiH,'Children');
set(findobj('Tag','pushOkay'),'UserData',1); % Sets the min distance allowed between the sliders
title='Plot trajectories originating between frames:';
set(findobj('Tag','editFirstFrame'),'String',num2str(first));
set(findobj('Tag','editLastFrame'),'String',num2str(last));
set(findobj('Tag','SelectFramesGUI'),'Name',title);
sSteps=[1/(last-first) 1/(last-first)];
set(findobj('Tag','sliderFirstFrame'),'SliderStep',sSteps,'Max',last,'Min',first,'Value',first);
set(findobj('Tag','sliderLastFrame'),'SliderStep',sSteps,'Max',last,'Min',first,'Value',last);
waitfor(guiH); % The function waits for the dialog to close (and to return values for uFirst and uLast)

% Check whether the dialog was closed (_CloseRequestFcn)
if uFirst==-1
    return % This will return an error (status=0)
end

% Open figure and axes
figure;
img=(img-max(img(:)))/(min(img(:))-max(img(:)));
imshow(img,[]); 
%h=axes;
hold on

% Draw trajectories
if method==1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % METHOD 1: Plot all matches per frames over the range uFirst:uLast
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Go through frames
    for i=uFirst:uLast
        Mo=M(:,:,i);
        Mv=Mo(find(Mo(:,1)~=0 & Mo(:,3)~=0),:);
        
        % Plot arrows
        quiver(Mv(:,2),Mv(:,1),Mv(:,4)-Mv(:,2),Mv(:,3)-Mv(:,1),0,'r')

    end 
    
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % METHOD 2: Plot entire trajectories ORIGINATING between frames uFirst and uLast (SLOW!)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i=uFirst:uLast
        
        % Extract all trajectories staring in the current frame
        Mo=M(:,:,i-1); % Check [0 0] in the previous frame
        [y,x]=find(Mo(:,1)==0 & Mo(:,3)~=0);
    
        % Plot them
        for j=1:length(y)
        
            % Extract all trajectories starting from here
            trj=M(y(j),:,i:end);
            
            ln=1;
            while ln<size(trj,3) & ~any(find(trj(1,:,ln)==0))
                plot([trj(1,2,ln) trj(1,4,ln)],[trj(1,1,ln) trj(1,3,ln)],'r');
                ln=ln+1;
            end
            
        end
        
    end
end

% Reformat axes
axis ij
axis equal
axis([1 imgSize(2) 1 imgSize(1)]);
