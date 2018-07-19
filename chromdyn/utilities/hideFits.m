function hideFits(figH)
%HIDEFIT hides all fits from trajectoryAnalysis in a figure
%
% SYNOPSIS hideFits(figH)
%
% INPUT    figH handle of the figure you want to turn the fits off
%           complement to unhideFits). Fits have to have been
%           drawn with trajectoryAnalysis.m or they have to be tagged 'TAfit'
%
% c: 05/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% is there a figure
if nargin == 0 | isempty(figH)
    figH = get(0,'CurrentFigure');
    if isempty(figH)
        error('no figure open')
    end
end

% find lineHandles in figure
lineHandles = findall(figH,'Type','line');

for lh = lineHandles'
    if findstr(get(lh,'Tag'),'TAfit')
        set(lh,'Visible','off')
    end
end