function fsmCentercb_drawPolygon
% fsmCentercb_drawPolygon
%
% This function is a callback for the "More Tools" menu in the figure
% opened by the pushImShow_Callback function of fsmCenter
%
% INPUT   None
%
% OUTPUT  None
%
% REMARK  A black and white mask (bwMask) and the drawn polygon are
%         returned into Matlab base workspace
%
% Aaron Ponti, 11/18/2003

% The user should draw a polygon
[bwMask,x,y]=roipoly;

% Plot the polygon
hold on;
plot(x,y,'r-');

% Return B&W mask and polygon to MATLAB base workspace
assignin('base','bwMask',bwMask);
assignin('base','polygon',[x y]);