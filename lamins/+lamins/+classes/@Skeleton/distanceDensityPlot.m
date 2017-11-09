function [ h ] = distanceDensityPlot( S, I , varargin)
%distanceDensityPlot Analyze intensity density versus density from lamin
%centerlines
% S lamins.classes.Skeleton
% I image

ip = inputParser;
ip.addParameter('metric','euclidean');
ip.addParameter('do3D',false);
ip.addParameter('showDensity',false);
ip.parse(varargin{:});

in = ip.Results;
if(strcmp(ip.Results.metric,'roundEuclidean'))
    in.metric = 'euclidean';
end

% keyboard;
dist = bwdist(S.bw,in.metric);
if(strcmp(ip.Results.metric,'roundEuclidean'))
    dist = round(dist);
end
mask = S.getMask;
[C,~,ic] = unique(dist(mask));
% Inte
dataByDist = accumarray(ic,mat2gray(I(mask)),[],@(x) {x});
h = figure;
if(ip.Results.showDensity || ip.Results.do3D)
    [f,xi] = cellfun(@(x) ksdensity(double(x)),dataByDist,'Unif',false);
end
hold on;
if(ip.Results.do3D)
    plotstats = @plot3d;
else
    plotstats = @(x,y) plot(x,y,'.-');
end
plotstats(C,cellfun(@(x) max(x),dataByDist));
plotstats(C,cellfun(@(x) prctile(x,99.5),dataByDist));
plotstats(C,cellfun(@(x) prctile(x,99),dataByDist));
plotstats(C,cellfun(@(x) prctile(x,95),dataByDist));
plotstats(C,cellfun(@(x) prctile(x,75),dataByDist));
plotstats(C,cellfun(@median,dataByDist));
plotstats(C,cellfun(@mean,dataByDist));
plotstats(C,cellfun(@(x) prctile(x,25),dataByDist));
plotstats(C,cellfun(@(x) prctile(x,5),dataByDist));
plotstats(C,cellfun(@(x) prctile(x,1),dataByDist));
plotstats(C,cellfun(@(x) prctile(x,0.5),dataByDist));
plotstats(C,cellfun(@(x) min(x),dataByDist));
if(ip.Results.showDensity)
    cellfun(@(x,y,z,l) plot3(repmat(x,size(y)),y,z,'-','Color',[l l l]),num2cell(C),xi,f,num2cell((0:length(C)-1)'/length(C)),'Unif',false);
end

xlim([0 max(C)]);
ylim([0 1]);
if(ip.Results.showDensity)
    zlim([0 prctile(cellfun(@max,f),90)])
end

legend({ ...
    'Max', ...
    '99.5th Percentile', ...
    '99th Percentile', ...
    '95th Percentile', ...
    '75th Percentile', ...
    'Median', ...
    'Mean', ...
    '25th Percentile', ...
    '5th Percentile', ...
    '1st Percentile', ...
    '0.5th Percentile', ...
    'Min' ...
    }, ...
'Location','northeastoutside' ...
);

xlabel([upper(ip.Results.metric(1)) ip.Results.metric(2:end) ' Distance from Lamins (pixels)'])
ylabel('Fluorescence Intensity (AU)');
zlabel('Fluorescence KSDensity');

    function plot3d(x,y)
        z = cellfun(@(x,y) ksdensity(double(x),y),dataByDist,num2cell(y));
        plot3(x,y,z,'o:');
    end

end

