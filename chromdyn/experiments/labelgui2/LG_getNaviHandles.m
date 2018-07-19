function [naviHandles, movieWindowHandles] = LG_getNaviHandles(noWarn)
%LG_getNaviHandles finds the full handle structure of the navigatior window
% movieWindowHandles.otherWindows contains all the handles to other
% labelgui2-windows
% naviHandles: navigator window
% movieWindowHandles: movie window

if nargin == 0
    noWarn = 0;
end


% find navigator window handle
naviH = findall(0,'Tag','LG_navigator');

% check that handle exists
if isempty(naviH)
    if ~noWarn
        h = errordlg('Lablegui Navigator window has been closed',...
            'Handle not found!');
        uiwait(h);
    end
    naviHandles = [];
    movieWindowHandles = [];
    return
end

% return handles
naviHandles = guidata(naviH);

% if there are two output arguments: Give also handles to current movie
% window
if nargout > 1
    % change here in case we ever allow two movies open at the same time
    if ~isempty(naviHandles.movieWindowH)
        movieWindowHandles = guidata(naviHandles.movieWindowH);
    else
        movieWindowHandles = [];
    end
end