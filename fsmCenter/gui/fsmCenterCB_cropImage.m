function fsmCenterCB_cropImage
% fsmCenterCB_cropImage
%
% This function is a callback for the "More Tools" menu in the figure
% opened by the pushImShow_Callback function of fsmCenter
%
% This function allows the user to crop the image displayed in the ViewPanel 
%
% INPUT   None
%
% OUTPUT  None
%
% Aaron Ponti, 08/13/2004

% Crop
figure(gcf);
[img,rect]=imcrop;

% Update figure
clf;
imshow(img,[]);

% Update title
figureTitle=get(gcf,'Name');
pos=findstr(figureTitle,' (');
imSize=size(img);
newTitle=[figureTitle(1:pos-1),' (',num2str(imSize(1)),'x',num2str(imSize(2)),')'];
set(gcf,'Name',newTitle);

% Print out the coordinates
fprintf(1,'Cropped region: y=%4.0f:%4.0f, x=%4.0f:%4.0f\n',rect(2),rect(1),rect(2)+rect(4),rect(1)+rect(3));

% Export to base workspace
assignin('base','img',img);

% Add common menus
fsmCenterCB_addToolsMenu(gcf);

% Create menu with additional tools
figure(gcf);
hMenu=uimenu('Label','More Tools');
uimenu(hMenu,'Label','Load image','Callback','fsmCenterCB_loadNewImage;','Accelerator','L');
uimenu(hMenu,'Label','Find edges','Callback','fsmCenterCB_findEdges;','Accelerator','E','Separator','On');


