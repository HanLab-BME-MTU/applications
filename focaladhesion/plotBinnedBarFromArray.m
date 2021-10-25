function [edges,Ydiscrete]=plotBinnedBarFromArray(X,Y,bin)
%function []=plotBinnedBarFromArray(X,Y,bin) plot a bar graph by binning X
%with bin.
%   Usage: plotBinnedBarFromArray(LF,PT,20)
edges = 0:bin:ceil(max(X)/50)*50;
Ydiscrete = discretize(X,edges);

YPerBin = arrayfun(@(x) nanmean(Y(Ydiscrete==x)),1:numel(edges));
pB = bar(edges+bin/2,YPerBin);
pB.FaceAlpha=0;
pB.EdgeColor='r';
pB.LineWidth = 2;
stdYPerBin = arrayfun(@(x) stdErrMean(Y(Ydiscrete==x)),1:numel(edges));
eH = errorbar(edges+bin/2,YPerBin, stdYPerBin,'Color', 'r');
eH.YNegativeDelta = [];
eH.LineStyle='none';
eH.LineWidth = 2;
end

