function patch_strainEnergyGUI_addCalciumControls(figPath)
% patch_strainEnergyGUI_addCalciumControls Add calcium-related UI controls to
% strainEnergyCalculationProcessGUI.fig (GUIDE figure).
%
% Usage:
%   patch_strainEnergyGUI_addCalciumControls('strainEnergyCalculationProcessGUI.fig')
%
% This makes a timestamped backup next to the fig and then overwrites fig.

if nargin<1 || isempty(figPath)
    figPath = 'strainEnergyCalculationProcessGUI.fig';
end
assert(exist(figPath,'file')==2, 'Fig not found: %s', figPath);

% Backup
[folder,name,ext] = fileparts(figPath);
ts = datestr(now,'yyyymmddTHHMMSS');
bak = fullfile(folder, sprintf('%s_backup_%s%s', name, ts, ext));
copyfile(figPath, bak, 'f');
fprintf('[Backup] %s\n', bak);

% Open fig invisibly
hFig = openfig(figPath,'invisible','reuse');
set(hFig,'Units','normalized');

% Helper: create if missing
makeIfMissing = @(tag, maker) createIfMissing(hFig, tag, maker);

% A small block in the lower-left area. Adjust if you want.
x0 = 0.05; y0 = 0.02; w = 0.40; h = 0.05;
dy = 0.055;

makeIfMissing('checkbox_readCalcium', @() uicontrol(hFig,'Style','checkbox', ...
    'Units','normalized','Position',[x0 y0+4*dy w h], ...
    'String','Read calcium intensity (optional)','Value',0, ...
    'Tag','checkbox_readCalcium'));

makeIfMissing('text_calciumChan', @() uicontrol(hFig,'Style','text', ...
    'Units','normalized','Position',[x0 y0+3*dy 0.18 h], ...
    'String','Ca channel','HorizontalAlignment','left', ...
    'Tag','text_calciumChan'));

makeIfMissing('edit_calciumChan', @() uicontrol(hFig,'Style','edit', ...
    'Units','normalized','Position',[x0+0.20 y0+3*dy 0.08 h], ...
    'String','1','BackgroundColor','white', ...
    'Tag','edit_calciumChan'));

makeIfMissing('checkbox_calciumDFF0', @() uicontrol(hFig,'Style','checkbox', ...
    'Units','normalized','Position',[x0 y0+2*dy w h], ...
    'String','Compute dF/F0','Value',1, ...
    'Tag','checkbox_calciumDFF0'));

makeIfMissing('text_calciumBaselineFrames', @() uicontrol(hFig,'Style','text', ...
    'Units','normalized','Position',[x0 y0+dy 0.18 h], ...
    'String','F0 frames','HorizontalAlignment','left', ...
    'Tag','text_calciumBaselineFrames'));

makeIfMissing('edit_calciumBaselineFrames', @() uicontrol(hFig,'Style','edit', ...
    'Units','normalized','Position',[x0+0.20 y0+dy 0.18 h], ...
    'String','1:10','BackgroundColor','white', ...
    'Tag','edit_calciumBaselineFrames'));

makeIfMissing('text_calciumMeasure', @() uicontrol(hFig,'Style','text', ...
    'Units','normalized','Position',[x0 y0 0.18 h], ...
    'String','Ca metric','HorizontalAlignment','left', ...
    'Tag','text_calciumMeasure'));

makeIfMissing('popupmenu_calciumMeasure', @() uicontrol(hFig,'Style','popupmenu', ...
    'Units','normalized','Position',[x0+0.20 y0 0.18 h], ...
    'String',{'mean','median','sum'}, 'Value',1, ...
    'Tag','popupmenu_calciumMeasure'));

% Save & close
savefig(hFig, figPath);
close(hFig);
fprintf('[Saved] Updated fig: %s\n', figPath);

end

function h = createIfMissing(hFig, tag, maker)
h = findall(hFig,'Tag',tag);
if isempty(h)
    h = maker();
else
    h = h(1);
end
end
