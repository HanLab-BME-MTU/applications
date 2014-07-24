function LG_changeTagColor
%LG_CHANGECOLOR changes the color for a tag
%
% SYNOPSIS: LG_changeColor
%
% INPUT none
%
% OUTPUT none
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 15-May-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get handles
[naviHandles, movieWindowHandles] = LG_getNaviHandles;

% get colormap
colorMap = movieWindowHandles.colorMap;


% find current tag:
% read tag index in userData of current object
tagIdx = get(gco,'UserData');

% launch colorselection dialog
oldColor = colorMap(tagIdx,:);
newColor = uisetcolor(oldColor);

if all(oldColor == newColor)
    % do nothing if user didn't select anything new
    return
end

% write colorMap
movieWindowHandles.colorMap(tagIdx,:) = newColor;

% save colorMap
guidata(movieWindowHandles.LG_movieWindow, movieWindowHandles);

% plot
LG_gotoFrame(LG_getCurrentTime);