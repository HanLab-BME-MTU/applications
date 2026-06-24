function fixLccbGuiFonts(varargin)
%FIXLCCBGUIFONTS  Shrink oversized fonts in lccb GUIDE windows on Linux.
%
% The lccb GUIDE GUIs (movieViewer, movieSelectorGUI, packageGUI, the u-track
% setting dialogs, ...) hard-code large FontSize values that render too big
% under Linux, so package names and labels get clipped. This CAPS only the
% oversized fonts at a maximum size and leaves already-small fonts (menus,
% list boxes) untouched, without editing any lccb source file.
%
% USAGE
%   fixLccbGuiFonts                      % cap fonts on all open figures
%   fixLccbGuiFonts(h)                   % cap a specific figure handle
%   fixLccbGuiFonts('maxSize', 9)        % cap large fonts at 9 pt (default 10)
%   fixLccbGuiFonts('minTouch', 11)      % only shrink fonts above this (default 11)
%   fixLccbGuiFonts('widen', true)       % also widen clipped text controls
%   fixLccbGuiFonts('auto', true)        % auto-apply to new GUIs (Linux)
%   fixLccbGuiFonts('force', true)       % apply even when not on Linux
%
% Only fonts strictly larger than minTouch are reduced, and never below
% maxSize, so small controls are left alone while big clipped labels shrink.
% Sangyoon J. Han / 2026

ip = inputParser;
ip.addOptional('h', [], @(x) isempty(x) || all(ishghandle(x)));
ip.addParameter('maxSize', 10, @(x) isnumeric(x) && isscalar(x) && x>0);
ip.addParameter('minTouch', 11, @(x) isnumeric(x) && isscalar(x) && x>0);
ip.addParameter('widen', false, @islogical);
ip.addParameter('force', false, @islogical);
ip.addParameter('auto', false, @islogical);
ip.parse(varargin{:});
h = ip.Results.h;
maxSize = ip.Results.maxSize; minTouch = ip.Results.minTouch;
widen = ip.Results.widen; force = ip.Results.force;

% only act on Linux unless forced (Windows/Mac render these fine)
if ~force && (~isunix || ismac)
    return;
end

% Auto mode: cap fonts on every new figure shortly after it is built.
if ip.Results.auto
    persistent autoTimer %#ok<TLEV>
    if isempty(autoTimer) || ~isvalid(autoTimer)
        autoTimer = timer('ExecutionMode','fixedSpacing','Period',1.0, ...
            'TimerFcn', @(~,~) rescaleNew(maxSize, minTouch, widen));
        start(autoTimer);
        fprintf('fixLccbGuiFonts: auto font-cap enabled (Linux).\n');
    end
    return;
end

if isempty(h)
    h = findall(0, 'Type', 'figure');
end
for f = reshape(h, 1, [])
    rescaleFigure(f, maxSize, minTouch, widen);
end
end

% =====================================================================
function rescaleFigure(f, maxSize, minTouch, widen)
% Cap oversized fonts in figure f. Fonts at or below minTouch are left as-is;
% fonts above minTouch are reduced to maxSize. An appdata flag prevents
% double processing in auto mode.
if ~ishghandle(f), return; end
if isappdata(f, 'lccbFontFixed') && getappdata(f, 'lccbFontFixed')
    return;
end
kids = findall(f, '-property', 'FontSize');
for c = reshape(kids, 1, [])
    try
        fs = get(c, 'FontSize');
        if isempty(fs) || ~isnumeric(fs), continue; end
        if fs > minTouch                 % only the oversized ones
            set(c, 'FontSize', maxSize);
            if widen, widenControl(c); end
        end
    catch
        % some objects expose FontSize but reject set; ignore
    end
end
try, setappdata(f, 'lccbFontFixed', true); catch, end %#ok<CTCH>
end

% =====================================================================
function widenControl(c)
% Give a clipped text control a little more width so the (now smaller) text
% fits. Only touches static text / labels with pixel-like units.
try
    if ~isprop(c, 'Position'), return; end
    u = '';
    if isprop(c,'Units'), u = get(c,'Units'); end
    if ~any(strcmpi(u, {'pixels','points','characters'})), return; end
    ty = '';
    if isprop(c,'Style'), ty = get(c,'Style'); end
    if ~any(strcmpi(ty, {'text','pushbutton','radiobutton','checkbox',''})), return; end
    pos = get(c,'Position');
    pos(3) = pos(3) * 1.25;          % 25% wider
    set(c,'Position',pos);
catch
end
end

% =====================================================================
function rescaleNew(maxSize, minTouch, widen)
figs = findall(0, 'Type', 'figure');
for f = reshape(figs, 1, [])
    if isappdata(f, 'lccbFontFixed') && getappdata(f, 'lccbFontFixed')
        continue;
    end
    rescaleFigure(f, maxSize, minTouch, widen);
end
end
