function fsmCentercb_loadNewImage
% fsmCentercb_loadNewImage
%
% This function is a callback for the "More Tools" menu in the figure
% opened by the pushImShow_Callback function of fsmCenter
%
% INPUT   None
%
% OUTPUT  None
%
% REMARK  A variable img is returned into Matlab base workspace
%
% Aaron Ponti, 11/18/2003

% Current directory
oldDir=cd;

% Read current image path from fsmCenter
hFsm=findall(0,'Tag','fsmCenter','Name','fsmCenter');
if ~isempty(hFsm)
    handles=guidata(hFsm);
    imagePath=get(handles.textCurrentImage,'String');
    if ~isempty(imagePath)
        if isdir(imagePath)
            cd(imagePath);
        end
    end
end

% Select image
[fileName,dirName] = uigetfile(...
    {'*.tif;*.tiff;*.jpg;*.jpeg','Image Files (*.tif,*.tiff,*.jpg,*.jpeg)';
    '*.tif','TIF files (*.tif)'
    '*.tiff','TIFF files (*.tiff)'
    '*.jpg;','JPG files (*.jpg)'
    '*.jpeg;','JPEG files (*.jpeg)'
    '*.*','All Files (*.*)'},...
    'Select image to load');
if(isa(fileName,'char') & isa(dirName,'char'))

    % Load image
    img=double(imread([dirName,fileName]));
    h=gcf;

    % Clear current window's children
    clf;
    
    % Display new image
    imshow(img,[]);
    
    % Make sure the axis are correct
    axis([1 size(img,2) 1 size(img,1)]);
    
    % Add title
    titleStr=[dirName,fileName,' (',num2str(size(img,1)),'x',num2str(size(img,2)),')'];
    set(h,'Name',titleStr,'NumberTitle','off');

    % Add univocal tag
    set(h,'Tag','ViewPanel');
    
    % Return image to MATLAB base workspace
    assignin('base','img',img);
else
    close(gcf);
    return
end

% Add common menus
fsmCenterCB_addToolsMenu(gcf);

% Create menu with additional tools
figure(gcf);
hMenu=uimenu('Label','More Tools');
uimenu(hMenu,'Label','Load image','Callback','fsmCenterCB_loadNewImage;','Accelerator','L');
uimenu(hMenu,'Label','Find edges','Callback','fsmCenterCB_findEdges;','Accelerator','E','Separator','On');

% Go back to old directory
cd(oldDir);
