function fsmCenterCB_addToolsMenu(handle)
% fsmCenterCB_addToolsMenu
%
% This function adds the menu "More Tools" to the figure with handle 'handle'
% allowing the functions of the fsmCenter "View Panel" to be called.
%
% INPUT   handle: handle of the figure to which the menu has to be attached
%
% OUTPUT  None
%
% Aaron Ponti, 04/07/2004

% Create menu with additional tools
figure(handle);
hMenu=uimenu('Label','Common Tools');
uimenu(hMenu,'Label','Draw polygon','Callback','fsmCenterCB_drawPolygon;','Accelerator','D');
uimenu(hMenu,'Label','Load polygon','Callback','fsmCenterCB_loadPolygon;');
uimenu(hMenu,'Label','Crop image','Callback','fsmCenterCB_cropImage;');
uimenu(hMenu,'Label','Filter image','Callback','fsmCenterCB_filterImage;','Accelerator','F');
uimenu(hMenu,'Label','Invert image','Callback','fsmCenterCB_invertImage;','Accelerator','I');
uimenu(hMenu,'Label','Show speckle selection','Callback','fsmCenterCB_showCands;','Accelerator','K');

