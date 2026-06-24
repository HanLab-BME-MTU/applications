function fixLccbGuiFonts(varargin)
%FIXLCCBGUIFONTS  Shrink oversized fonts in lccb GUIDE windows on Linux.
%
% The lccb GUIDE GUIs (movieViewer, movieSelectorGUI, packageGUI, the u-track
% setting dialogs, ...) hard-code FontSize values that render too large under
% Linux, so package names and labels get clipped. This rescales the fonts of
% the relevant figures' UI controls without touching any lccb source file.
%
% USAGE
%   fixLccbGuiFonts                      % rescale all currently open figures
%   fixLccbGuiFonts(h)                   % rescale a specific figure handle
%   fixLccbGuiFonts('scale', 0.8)        % custom scale factor (default 0.78)
%   fixLccbGuiFonts('force', true)       % apply even when not on Linux
%
% TIP: to apply automatically, call this right after opening a GUI, e.g.
%   movieSelectorGUI; fixLccbGuiFonts;
% or add 'fixLccbGuiFonts' to your startup.m so it can be run on demand.
% Sangyoon J. Han / 2026

ip = inputParser;
ip.addOptional('h', [], @(x) isempty(x) || all(ishghandle(x)));
ip.addParameter('scale', 0.78, @(x) isnumeric(x) && isscalar(x) && x>0);
ip.addParameter('minSize', 7, @isnumeric);     % never go below this (points)
ip.addParameter('force', false, @islogical);
ip.addParameter('auto', false, @islogical);    % install listener for new GUIs
ip.parse(varargin{:});
h = ip.Results.h; scale = ip.Results.scale;
minSize = ip.Results.minSize; force = ip.Results.force;

% only act on Linux unless forced (Windows/Mac render these fine)
if ~force && (~isunix || ismac)
    return;
end

% Auto mode: every figure created from now on is rescaled shortly after it
% finishes building. Call 'fixLccbGuiFonts('auto',true)' once (e.g. in
% startup.m) and you never have to call it manually again.
if ip.Results.auto
    persistent autoTimer %#ok<TLEV>
    if isempty(autoTimer) || ~isvalid(autoTimer)
        autoTimer = timer('ExecutionMode','fixedSpacing','Period',1.0, ...
            'TimerFcn', @(~,~) rescaleNew(scale, minSize));
        start(autoTimer);
        fprintf('fixLccbGuiFonts: auto-rescale enabled (Linux).\n');
    end
    return;
end

if isempty(h)
    h = findall(0, 'Type', 'figure');
end

for f = reshape(h, 1, [])
    rescaleFigure(f, scale, minSize);
end
end

% =====================================================================
function rescaleFigure(f, scale, minSize)
% Shrink every FontSize-bearing control in figure f by 'scale'. Marks the
% figure as done (appdata flag) so auto mode does not rescale it twice.
if ~ishghandle(f), return; end
if isappdata(f, 'lccbFontFixed') && getappdata(f, 'lccbFontFixed')
    return;   % already rescaled
end
kids = findall(f, '-property', 'FontSize');
for c = reshape(kids, 1, [])
    try
        fs = get(c, 'FontSize');
        if isempty(fs) || ~isnumeric(fs), continue; end
        newfs = max(minSize, round(fs * scale));
        if newfs < fs
            set(c, 'FontSize', newfs);
        end
    catch
        % some objects expose FontSize but reject set; ignore
    end
end
try, setappdata(f, 'lccbFontFixed', true); catch, end %#ok<CTCH>
end

% =====================================================================
function rescaleNew(scale, minSize)
% Timer callback for auto mode: rescale any figure not yet handled. The
% appdata flag set by rescaleFigure prevents repeated work, so this is cheap
% even when run once per second.
figs = findall(0, 'Type', 'figure');
for f = reshape(figs, 1, [])
    if isappdata(f, 'lccbFontFixed') && getappdata(f, 'lccbFontFixed')
        continue;
    end
    rescaleFigure(f, scale, minSize);
end
end
