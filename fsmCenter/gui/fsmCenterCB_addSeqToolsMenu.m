function fsmCenterCB_addSeqToolsMenu(handle)
% fsmCenterCB_addSeqToolsMenu
%
% This function adds the menu "More Tools" to the figure with handle 'handle'
% allowing the functions of the fsmCenter "View Panel" to be called.
%
% INPUT   handle: handle of the figure to which the menu has to be attached
%
% OUTPUT  None
%
% Andre Kerstens, 12/15/2004

% Create menu with additional tools
figure(handle);
hMenu=uimenu('Label','Common Tools');
uimenu(hMenu,'Label','Filter image','Callback','fsmCenterCB_filterSeqImage;','Accelerator','F');
uimenu(hMenu,'Label','Invert image','Callback','fsmCenterCB_invertSeqImage;','Accelerator','I');
uimenu(hMenu,'Label','Show speckle selection','Callback','fsmCenterCB_showSeqCands;','Accelerator','K');

