function fsmCenterCB_invertImage
% fsmCenterCB_invertImage
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
% Aaron Ponti, 01/22/2004

% Retrieve image from figure
hImg=findall(gca,'Type','Image');
img=get(hImg,'CData');

% Invert
iImg=(img-max(img(:)))/(min(img(:))-max(img(:)));

% Refresh figure
set(hImg,'CData',iImg);
set(gca,'CLimMode','manual');
set(gca,'CLim',[0 1]);
refresh;

% Export to base workspace
assignin('base','iImg',iImg);
