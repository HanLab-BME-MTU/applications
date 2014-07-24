function figureHandle = possibilisticClusteringPlot(varargin)
%POSSIBILISTICCLUSTERINGPLOT plots the membership function overlaid over data points
%
% SYNOPSIS: possibilisticClusteringPlot(data,membership,centers,imSize,membershipSelect, classification, initialGuess, trueCenters)
%           possibilisticClusteringPlot(ax,...)
%
% INPUT data: Data input for clustering. Has to be 2-d "coordinates", but
%               can be MDS coordinates
%		membership (opt): membership array from possibilisticClustering.m
%       centers (opt): structure array with field .mean for every center
%       imSize (opt): [nX, nY] number of bins in x,y. Default [50,50]
%       membershipSelect (opt)
%              1: plot for every pixel the maximum membership (def)
%              2: plot for every pixel the minimum membership
%       classification (opt): array of the size of membership (or with one
%               additional column) with 1 at ij if data point i belongs to
%               cluster j. The additional column is for noise points.
%       initialGuess (opt): nCenters-by-2 array with initial center guesses
%       trueCenters (opt): nCenters-by-2 array with true centers
%       ax: handles to plot axes
%
% OUTPUT figureHandle:  figure handle
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 24-Oct-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def_membershipSelect = 1;
cLength = 256; %lenght of colormap

% check input
if nargin < 2
    error('not enough input arguments');
end

potentialHandle = varargin{1};
if isempty(potentialHandle)
    error('can''t have first argument empty');
end
if isscalar(potentialHandle) && ishandle(potentialHandle)
    ah = potentialHandle;
    varargin(1) = [];
else
    ah = [];
end


data = returnRightVector(varargin{1},2);
nData = size(data,1);
if length(varargin)<2 || isempty(varargin{2})
    doMap = false;
else
    doMap = true;
    membership = returnRightVector(varargin{2},nData,'r');
end

if length(varargin)<3 || isempty(varargin{3}) ||...
        ~isstruct(varargin{3}) || ~isfield(varargin{3},'mean')
    centers = [];
else
    centers = varargin{3};
    centers = cat(1,centers.mean);
end


if length(varargin) < 5 || isempty(varargin{5})
    membershipSelect = def_membershipSelect;
else
    membershipSelect = varargin{5};
end

if length(varargin)< 6 || isempty(varargin{6})
    classification = true(nData,1);
    if ~doMap
        error('need either membership or classification')
    end
else
    classification = logical(varargin{6});
end
if length(varargin)< 7 || isempty(varargin{7})
    initialGuess = [];
else
    initialGuess = (varargin{7});
end
if length(varargin)< 8 || isempty(varargin{8})
    trueCenters = [];
else
    trueCenters = (varargin{8});
end

if doMap
    switch membershipSelect
        case 2
            mapMembership = min(membership,[],2);
        case 1
            mapMembership = max(membership,[],2);
        otherwise
            error('membershipSelect %i has not been implemented',...
                membershipSelect);
    end

    % count nClusters so that there is no confusion with centers
    nClusters = size(membership,2);
else
    nClusters = size(classification,2)-1;
end

% % norm data to [0,1]
% oldData = data;
% minData = min(data,[],1);
% data = data - repmat(minData,nData,1);
% maxData = max(data,[],1);
% data = data./repmat(maxData,nData,1);

% norm data to [0,1]
oldData = data;
minData = min(data(:),[],1) * [1,1];
data = data - repmat(minData,nData,1);
maxData = max(data(:),[],1) * [1,1];
data = data./repmat(maxData,nData,1);

% check bins
if length(varargin) < 4 || isempty(varargin{4})
    imSize = [50,50];
else
    imSize = varargin{4};
end

% multiply by bins
data = data.*repmat(imSize-1,nData,1)+1;

if doMap

    % create map
    grid = [(1:imSize(1))',(1:imSize(2))'];
    map = imDataMap(imSize,data,mapMembership,'infLen',min(imSize)/5,...
        'gridSmoothing',min(imSize)/15,'grid',grid);
    % normalize map to cLength
    minMap = min(map(:));
    maxMap = max(map(:))-minMap;
    map = round((map-minMap)/(maxMap)*cLength-1)+1;

    % if nClusters < 4, we can do individual contours
    cmap = flipud(gray(cLength));
    if nClusters < 4
        individualMaps = zeros([imSize,nClusters]);
        colors = {'rw','gw','bw'};
        % loop up to 3 times
        for c = 1:nClusters
            tmp = ...
                imDataMap(imSize,data,membership(:,c),...
                'infLen',min(imSize)/5,'gridSmoothing',...
                min(imSize)/15,'grid',grid);
            % normalize individual maps
            individualMaps(:,:,c) = ...
                round((tmp-min(tmp(:)))/(max(tmp(:))-min(tmp(:)))*cLength-1)+1;
            cmap = [cmap;isomorphicColormap(colors{c},cLength)];

        end

    end % if less than 4 clusters
else
    cmap = jet(64);
end % if doMap

% Create figure
if isempty(ah)
    fh = figure('Colormap',cmap);
else
    fh = get(ah,'Parent');
    set(fh,'Colormap',cmap);
    axes(ah);
end

if doMap
    % Create contour
    contour(linspace(minData(1),maxData(1)+minData(1),imSize(1)),...
        linspace(minData(2),maxData(2)+minData(2),imSize(2)),...
        map',...
        'Fill','on',...
        'LevelStep',1,...
        'LineStyle','none');
    % add contourlines if possible
    if nClusters < 4
        hold on
        for c=1:nClusters
            contour(linspace(minData(1),maxData(1)+minData(1),imSize(1)),...
                linspace(minData(2),maxData(2)+minData(2),imSize(2)),...
                individualMaps(:,:,c)'+c*cLength,'LevelStep',cLength/10)
        end
    end
end

% plot data points
hold on
% loop according to classification
nClasses = size(classification,2);
if nClusters < 4
    colors = 'rgb';
else
    colors = repmat('k',1,nClusters);
end
if nClasses == 1
    plot(oldData(:,1),oldData(:,2),'o','MarkerFaceColor','k',...
        'MarkerEdgeColor','w','MarkerSize',4)
else
    % loop to plot
    markers = 's^dv<>ph';
    markerSize = [4,5,5,5,5,5,5,6];
    nMarkers = length(markers);


    for c=1:nClusters
        markerIdx = mod(c-1,nMarkers)+1;

        plot(oldData(classification(:,c),1),...
            oldData(classification(:,c),2),markers(markerIdx),...
            'MarkerFaceColor',colors(c),...
            'MarkerEdgeColor','w','MarkerSize',markerSize(markerIdx))
    end

    % add noise-data
    if nClasses>nClusters
        plot(oldData(classification(:,end),1),...
            oldData(classification(:,end),2),'o',...
            'MarkerFaceColor','k',...
            'MarkerEdgeColor','w','MarkerSize',4)
    end


end
if ~isempty(centers)
    if nClusters < 4
        for c = 1:nClusters
            plot(centers(c,1),centers(c,2),['o',colors(c)])
            plot(centers(c,1),centers(c,2),['+',colors(c)])
        end
    else
        plot(centers(:,1),centers(:,2),'og')
    end
end
if ~isempty(initialGuess)
    if nClusters < 4
        for c = 1:nClusters
            plot(initialGuess(c,1),initialGuess(c,2),['o',colors(c)],...
                'MarkerSize',12)
        end
    else
        plot(initialGuess(:,1),initialGuess(:,2),'og',...
            'MarkerSize',12)
    end
end
if ~isempty(trueCenters)
    if nClusters < 4
        for c = 1:nClusters
            plot(trueCenters(c,1),trueCenters(c,2),['+',colors(c)],...
                'MarkerSize',16)
        end
    else
        plot(trueCenters(:,1),trueCenters(:,2),'+g',...
            'MarkerSize',16)
    end
end
% make axis square
axis square

if nargout > 0
    figureHandle = fh;
end
