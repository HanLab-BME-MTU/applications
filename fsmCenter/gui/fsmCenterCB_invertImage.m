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
children=get(gca,'Children');
if length(children)>1
    children=children(end); % We take the last one, which is the image (first in last out)
end
img=get(children,'CData');

% Invert
iImg=(img-max(img(:)))/(min(img(:))-max(img(:)));

% Refresh figure
set(children,'CData',iImg);
set(gca,'CLimMode','manual');
set(gca,'CLim',[0 1]);
refresh;

% Export to base workspace
assignin('base','iImg',iImg);
