function fsmCentercb_loadNewImage(img)
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

    % Return image to MATLAB base workspace
    assignin('base','img',img);
else
    close(gcf);
    return
end

% Create menu with additional tools
figure(gcf);
hMenu=uimenu('Label','More Tools');
uimenu(hMenu,'Label','Load image','Callback','fsmCenterCB_loadNewImage;','Accelerator','L');
uimenu(hMenu,'Label','Draw polygon','Callback','fsmCenterCB_drawPolygon;','Accelerator','D','Separator','On');
uimenu(hMenu,'Label','Find edges','Callback','fsmCenterCB_findEdges;','Accelerator','E');
uimenu(hMenu,'Label','Filter image','Callback','fsmCenterCB_filterImage;','Accelerator','F');
% uimenu(hMenu,'Label','Auto contrast','Callback','fsmCenterCB_autoContrast;','Accelerator','T');
uimenu(hMenu,'Label','Show speckle selection','Callback','fsmCenterCB_showCands;','Accelerator','K','Separator','On');
uimenu(hMenu,'Label','Invert image','Callback','fsmCenterCB_invertImage;','Accelerator','I');

