function fsmCenterCB_invertSeqImage
% fsmCenterCB_invertSeqImage
%
% This function is a callback for the "More Tools" menu in the figure
% opened by the pushImShow_Callback function of fsmCenter. Used for the
% sequence of images
%
% INPUT   None
%
% OUTPUT  None
%
% REMARK  A variable img is returned into Matlab base workspace
%
% Andre Kerstens, 12/15/2004

% Get handle to the invert menu
hMenu = findobj('Label','Invert image');

if strcmp(get(hMenu, 'Checked'),'on')
    set(hMenu, 'Checked', 'off');
else 
    set(hMenu, 'Checked', 'on');
end

% Retrieve image from figure
hImg=findall(gca,'Type','Image');
img=get(hImg,'CData');

% Invert
iImg=(img-max(img(:)))/(min(img(:))-max(img(:)));

% Refresh figure
set(hImg,'CData',iImg);
%set(gca,'CLimMode','manual');
%set(gca,'CLim',[0 1]);
refresh;

% Export to base workspace
assignin('base','iImg',iImg);