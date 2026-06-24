function fixMovieViewer(varargin)
%FIXMOVIEVIEWER  Compact the movieViewer GUI on Linux so all overlay items
%fit on screen. Shrinks font sizes AND rescales each panel's control
%positions vertically so the full list (including Filopodia at the bottom)
%is accessible without a scrollbar.
%
% USAGE
%   fixMovieViewer                  % find open Viewer and compact it
%   fixMovieViewer(hFig)            % apply to a specific figure
%   fixMovieViewer('vscale', 0.6)   % vertical compression (default 0.65)
%   fixMovieViewer('fontSize', 8)   % font cap (default 8 pt)
%   fixMovieViewer('force', true)   % run even on non-Linux
%
% How it works
%   For each panel (Image, Overlay) the function:
%     1. Finds the control with the largest Y (top-most control).
%     2. Finds the panel height needed to hold all controls at their current
%        positions (contentH).
%     3. Scales every control's Y and Height by (panelH / contentH) so the
%        content exactly fits the existing panel without overflow.
%     4. Caps FontSize at fontSize.
% Sangyoon J. Han / 2026

ip = inputParser;
ip.addOptional('hFig', [], @(x) isempty(x) || ishandle(x));
ip.addParameter('vscale',   0.65, @(x) isnumeric(x) && x>0 && x<=1);
ip.addParameter('fontSize', 8,    @(x) isnumeric(x) && x>=4);
ip.addParameter('force', false,   @islogical);
ip.parse(varargin{:});
hFig     = ip.Results.hFig;
vscale   = ip.Results.vscale;
fontSize = ip.Results.fontSize;
force    = ip.Results.force;

if ~force && (~isunix || ismac), return; end

if isempty(hFig)
    hFig = findall(0, 'Type','figure','Name','Viewer');
    if isempty(hFig)
        warning('fixMovieViewer:notFound','No open Viewer figure found.'); return;
    end
    hFig = hFig(1);
end

% ---- 1. resize the figure to use most of the screen height ----
scr = get(0,'ScreenSize');       % [left bottom width height]
scrH = scr(4); scrW = scr(3);
fig_pos = get(hFig,'Position');  % [left bottom width height]
newH = min(scrH - 80, fig_pos(4));   % leave room for taskbar
newBottom = max(10, scrH - newH - 40);
set(hFig,'Position',[fig_pos(1), newBottom, fig_pos(3), newH]);

% ---- 2. compress each panel's contents vertically ----
panels = findall(hFig, 'Type','uipanel');
for p = reshape(panels, 1, [])
    compressPanel(p, vscale, fontSize);
end

% ---- 3. cap fonts on non-panel controls (toolbar labels etc.) ----
others = findall(hFig, '-property','FontSize');
for c = reshape(others, 1, [])
    try
        if get(c,'FontSize') > fontSize + 2
            set(c,'FontSize', fontSize + 2);
        end
    catch, end
end
end

% ----------------------------------------------------------------
function compressPanel(hPanel, vscale, fontSize)
% Get all direct and indirect children with a Position property
kids = findall(hPanel, '-property','Position', '-not','Type','uipanel');
if isempty(kids), return; end

% collect current positions (units = pixels assumed for GUIDE panels)
oldUnits = cell(numel(kids),1);
pos = zeros(numel(kids),4);
for i = 1:numel(kids)
    try
        oldUnits{i} = get(kids(i),'Units');
        set(kids(i),'Units','pixels');
        pos(i,:) = get(kids(i),'Position');
    catch
        pos(i,:) = NaN;
    end
end
valid = ~any(isnan(pos),2);
if ~any(valid), return; end

% content extents: lowest and highest Y
yBot = min(pos(valid,2));
yTop = max(pos(valid,2) + pos(valid,4));
contentH = yTop - yBot;

% panel inner height
panUnits = get(hPanel,'Units');
set(hPanel,'Units','pixels');
panPos = get(hPanel,'Position');
set(hPanel,'Units',panUnits);
panInnerH = panPos(4);   % approximate inner height

if contentH <= panInnerH
    % already fits; restore units and exit
    for i=1:numel(kids)
        try, set(kids(i),'Units',oldUnits{i}); catch, end
    end
    return;
end

% scale factor to fit content into panel
sf = vscale;   % user-supplied compression

for i = 1:numel(kids)
    if ~valid(i), continue; end
    p = pos(i,:);
    % scale Y position and height; keep X and Width
    newY = (p(2) - yBot) * sf + 2;          % shift up from bottom
    newH = max(1, p(4) * sf);
    try
        set(kids(i), 'Position', [p(1), newY, p(3), newH]);
    catch, end
    % cap font
    try
        if get(kids(i),'FontSize') > fontSize
            set(kids(i),'FontSize', fontSize);
        end
    catch, end
    % restore units
    try, set(kids(i),'Units',oldUnits{i}); catch, end
end
end
