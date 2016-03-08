function [ idxWindInFinal ] = GCASubRegionalGetWindows(windows,roiMask )
%GCASubRegionalGetWindows

% bummer but need to test if empty - easiest way I could think to quickly
% do this.
xyCoordsWindEdge = cell(numel(windows),1);
linIdxWindEdge = cell(numel(windows),1);
for iWind = 1:numel(windows)
    if ~isempty(windows{iWind});
        xyCoordsWindEdge{iWind,1} = windows{iWind}{1}{end}; % get the edge coords
        % change to linear indices ( just easier to do quick compare)
        linIdxWindEdge{iWind,1} = sub2ind(size(roiMask), ...
            round(xyCoordsWindEdge{iWind,1}(2,:)),round(xyCoordsWindEdge{iWind,1}(1,:)));
    end
end
% get the coordinates of the outer windows - couldn't do it hte below way because erroring for empty winds.
%xyCoordsWindEdge = arrayfun(@(x) windows{x}{1}{end},1:numel(windows)-1,'uniformoutput',0);

% convert to linearIdx
%linIdxWindEdge = cellfun(@(x) sub2ind(size(roiMask),round(x(2,:)),round(x(1,:))),xyCoordsWindEdge,'uniformoutput',0);

linIdxMask = find(roiMask);

% check for overlap with roiMaks
test = cellfun(@(x) intersect(x,linIdxMask),linIdxWindEdge,'uniformoutput',0);
idxWindInFinal = cellfun(@(x) ~isempty(x),test);
% idxWindIn = arrayfun(@(x) length(test{x})==length(unique(linIdxWindEdge{x})),1:numel(test));
% idxWindIn2 = arrayfun(@(x) length(test{x}) == length(unique(linIdxWindEdge{x}))-1,1:numel(test));
% idxEmpty = arrayfun(@(x) isempty(linIdxWindEdge{x}), 1:numel(test));
%idxWindInFinal = (idxWindIn | idxWindIn2 ) & ~idxEmpty ;
end

