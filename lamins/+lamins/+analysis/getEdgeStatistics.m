function [ stat ] = getEdgeStatistics( skeletons, property, scale )
%plotFaceStatistics
%
% skeletons - an array of lamins.classes.Skeletons
% scale - normalization factor to scale statistic by

% no parameters, present gui
if(nargin < 1 || isempty(skeletons))
    [filename, pathname] = uigetfile('*_list.mat');
    skeletons = [pathname filesep filename];
end

% default property is Area
if(nargin < 2)
    property = 'Area';
end

% default scale is 1000 nm.^2
if(nargin < 3)
    scale = 1000;
end

if(ischar(skeletons))
    skeletons = MovieList.load(skeletons);
end
if(isa(skeletons,'MovieList'))
    data = lamins.functions.aggregateMovieListData(skeletons);
    if(isfield(data,'Sf2'))
        stat(1) = lamins.analysis.getEdgeStatistics([data.Sf], property, scale);
        stat(2) = lamins.analysis.getEdgeStatistics([data.Sf2], property, scale);
        return;
    else
        skeletons = [data.Sf];
    end
end

switch(property)
    case 'EdgesPerVertex'
        out = arrayfun(@(S) S.getEdgesPerVertex,skeletons,'UniformOutput',false);
        rps(numel(skeletons)) = struct(property,[]);
        [rps.(property)] = out{:};
        rps = num2cell(rps);
    otherwise
        rps = arrayfun(@(S) regionprops(S.edges,property),skeletons,'UniformOutput',false);
end

stat.all.rp = vertcat(rps{:});
[stat.all.data] = [stat.all.rp.(property)]/scale;
stat.all.median = median(stat.all.data(:));
stat.all.std = std(stat.all.data(:));
stat.all.mean = mean(stat.all.data(:));
stat.all.top = prctile(stat.all.data(:),95);
stat.all.cutoff = stat.all.median + 3*stat.all.std;
% stat.all.cutoff = 1;
stat.all.outliers = stat.all.data > stat.all.cutoff;
[stat.all.N,stat.all.edges] = histcounts(stat.all.data(~stat.all.outliers));
stat.all.cutoff = max(stat.all.cutoff,stat.all.edges(end));
stat.all.histmax = max(stat.all.N);
stat.all.centers = ( stat.all.edges(1:end-1) + stat.all.edges(2:end) )/2;

stat.movie(length(skeletons)) = struct();

for i=1:length(skeletons)
    [stat.movie(i).data] = [rps{i}.(property)]/scale;
    stat.movie(i).median = median(stat.movie(i).data(:));
    stat.movie(i).std = std(stat.movie(i).data(:));
    stat.movie(i).mean = mean(stat.movie(i).data(:));
    stat.movie(i).top = prctile(stat.movie(i).data(:),95);
    stat.movie(i).outliers = stat.movie(i).data > stat.all.cutoff;
    [stat.movie(i).N,stat.movie(i).edges] = histcounts(stat.movie(i).data(~stat.movie(i).outliers),stat.all.edges);
    stat.movie(i).histmax = max(stat.movie(i).N);
    stat.movie(i).centers = ( stat.movie(i).edges(1:end-1) + stat.movie(i).edges(2:end) )/2;

    
end


end

