function unhideFits(figH)
%UNHIDEFITs shows all fits from trajectoryAnalysis in a figure
%
% SYNOPSIS unhideFits(figH)
%
% INPUT    figH handle of the figure you want to turn the fits on
%           complement to hideFits). Fits have to have been
%           drawn with trajectoryAnalysis.m or they have to be tagged 'TAfit'
%
% c: 05/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find lineHandles in figure
lineHandles = findall(figH,'Type','line');

for lh = lineHandles'
    if findstr(get(lh,'Tag'),'TAfit')
        set(lh,'Visible','on')
    end
end