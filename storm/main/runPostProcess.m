%% LAUNCH IMARIS

im = Imaris();
im.setupScene();
dis = Show(data,im);

disp('------------------------------------------------');

% DISPLAY MODELS

dis.points
% dis.pointsRaw
% dis.pointsWithWidth
dis.models
% dis.projections
% dis.maxGaps
% dis.nullCluster
% dis.clusters
% dis.modelGroundTruth

disp('------------------------------------------------');

%% POST-PROCESS
pro = Processor(data);
thresholdLength = 150;
pro.dissolveModelsShorterThan(thresholdLength);
densityThreshold = 0.05;
pro.dissolveModelsLessDenseThan(densityThreshold);
% pro.dissolveClustersSmallerThan(5);

disp('------------------------------------------------');

%% ANALYZER MAIN
ana = Analysis(data);

% % % % Maximum gap
% % % figure(1);
% % % ana.maxGap();

% Residuals
% % % figure(2);
% % % rms = ana.residuals();
% % % [sigX,sigY,sigZ] = ana.residualsXYZ()

% % % % Bundles
% % % edgeRadius = 100;
% % % distanceThreshold = 100;
% % % % minOverlap = 0.50;
% % % % ana.bundles(im,edgeRadius,distanceThreshold,minOverlap);
% % % minOverlapDistance = 30;
% % % ana.bundles(im,edgeRadius,distanceThreshold,minOverlapDistance);

% % % % Branches
% % % edgeRadius = 50;
% % % % minAngle = 65;
% % % % maxAngle = 75;
% % % minAngle = 65;
% % % maxAngle = 75;
% % % maxDistance = 20;
% % % figure(3);
% % % ana.branches(im,edgeRadius,minAngle,maxAngle,maxDistance);

% % % % Maximum curvature
% % % figure(4);
% % % ana.maxCurvature();

% % % % Spacings
% % % figure(5);
% % % medianSpacing = ana.spacings();

% Point density
figure(6);
ana.density();

% % % % Intensity
% % % figure(7);
% % % ana.intensity();

% % % % Model length
% % % figure(8);
% % % averageLength = ana.length()

% % % % Cluster size
% % % figure(9);
% % % ana.clusterSize();

% % % % Statistics
% % % ana.statistics();

% % % % Reconstruction score
% % % sampleInterval = 5; edgeRadius = 100;
% % % ana.reconstructionScore(sampleInterval,edgeRadius);

% % % % Width histogram
% % % figure(10);
% % % ana.width();

% % % % Orientation
% % % scalingFactor = 2000;
% % % smoothing = 5;
% % % displayText = false;
% % % maxHistValue = 4000;
% % % zoom = 1.0;
% % % figure(11);
% % % orientation = ana.orientation(smoothing,displayText,maxHistValue,zoom);


disp('------------------------------------------------');

%% CREATE MAP

% Get a list with all the .dat-files in _map
path = 'Y:\fsm\harvard\data\Zhuang\_map\';
list = dir(path);
itemIdx = 0;
nDataSets = numel(list)-2;
items = cell(nDataSets,2);
for i=3:numel(list)
    if list(i).isdir
        sublist = dir([path list(i).name '\']);
        for k=3:numel(sublist)
            if strcmp(sublist(k).name(end-3:end),'.cfg')
                % Config found - look for corresponding data file
                for j=3:numel(sublist)
                    if strcmp(sublist(j).name(end-5:end),'.p.dat')
                        itemIdx = itemIdx + 1;
                        items{itemIdx,1} = [path list(i).name '\' sublist(k).name];
                        items{itemIdx,2} = [path list(i).name '\' sublist(j).name];
                        break;
                    end
                end
                break;
            end
        end
    end
end
items = items(1:itemIdx,:);

% Load preview file
cfg = Config();
cfg = cfg.loadWithFullPath(items{1,1});
roiSel = ROISelector();
roiSel = roiSel.load([cfg.path,cfg.fileName]);
minX = roiSel.minX;
minY = roiSel.minY;
scaleFac = roiSel.previewScaleFactor;

% Create figure 
hFig = figure(150);
pos = get(hFig,'Position');
set(hFig,'Position',[100 100 size(roiSel.preview)]);

% Plot background
hAxes = axes();
set(hAxes,'Position',[0 0 1 1]);
% imagesc(roiSel.preview(:,end:-1:1));
imagesc(roiSel.preview);
% imagesc(ones([size(roiSel.preview),3]).*repmat(reshape([1,0,1],1,1,3),[size(roiSel.preview),1]));

% Loop through all the files
for i=1:size(items,1)
    % Read the .dat-file and .cfg-file
    cfg = Config();
    data = Data();
    cfg = cfg.loadWithFullPath(items{i,1});
    data = data.loadWithFullPath(items{i,2});

    % Read the plot position and convert the units
    siz = data.roiSize(1:2)*scaleFac./[size(roiSel.preview,2) size(roiSel.preview,1)];
    pos = ((data.roiPosition(1:2)-[minX,minY])*scaleFac)./[size(roiSel.preview,2) size(roiSel.preview,1)];   
    
    pos(2) = 1-pos(2);
    pos = pos-[0 1].*siz;
    
    % Create subplot
    h = axes();
    scale = 1;
    set(h,'Position',[pos-scale/4*siz scale*siz]);
    
    % Compute the statistic   
    ana = Analysis(data);
    smoothing = 5;
    displayText = false;
    maxHistValue = 10000;
    zoom = 1.0;
    ana.predominantOrientationXY(smoothing,displayText,maxHistValue,zoom);
%     ana.length();
    
    delete(findall(h,'type','text'));
    set(h,'xticklabel',[]);
    set(h,'yticklabel',[]);
    drawnow;
    
    if i==1
        xLim = get(h,'XLim');
        yLim = get(h,'YLim');
    else
        set(h,'XLim',xLim);
        set(h,'YLim',yLim);
    end
end

disp('------------------------------------------------');

%% LINK CLUSTERS

clc;

maxDist = 80; % 80
closeGapSize = 200; % 200
angleTol = 30; % 30

% Get all the model edges up to maxDist
pro = Processor(data);
pro.initEdges(2*closeGapSize+maxDist);
pro.updateEdges();

% Clear edges leading to a model-less cluster
idx = ismember(data.edges(:,1),0);
edges = data.edges(~idx,:);
idx = ismember(edges(:,2),0);
edges = edges(~idx,:);

nEdges = size(edges,1);

linksStart = [0 0 0];
linksEnd = [0 0 0];
linksEdges = [0 0];
linksWeights = 0;

% Loop through all the edges
for e=1:nEdges
    clusterIdx1 = edges(e,1);
    clusterIdx2 = edges(e,2);
    
    % Bezier curve control points
    cP1 = data.modelBezCP{clusterIdx1};
    cP2 = data.modelBezCP{clusterIdx2};

    % Computing the tangent vectors at the endpoints pointing away from the
    % model
    t = [0,1];
    [~,normalT1] = tangentBezier(cP1,t');
    normalT1(1,:) = -normalT1(1,:);
    [~,normalT2] = tangentBezier(cP2,t');
    normalT2(1,:) = -normalT2(1,:);
    
    % Compute the distance
    cP1StartA = cP1(1,:)+closeGapSize*normalT1(1,:);
    cP1EndA = cP1(end,:)+closeGapSize*normalT1(end,:);
    cP2StartA = cP2(1,:)+closeGapSize*normalT2(1,:);
    cP2EndA = cP2(end,:)+closeGapSize*normalT2(end,:);
    %     cP1StartB = cP1(1,:)-closeGapSize*normalT1(1,:);
    %     cP1EndB = cP1(end,:)-closeGapSize*normalT1(end,:);
    %     cP2StartB = cP2(1,:)-closeGapSize*normalT2(1,:);
    %     cP2EndB = cP2(end,:)-closeGapSize*normalT2(end,:);
    backGap = 0.25;
    cP1StartB = cP1(1,:)-closeGapSize*backGap*normalT1(1,:);
    cP1EndB = cP1(end,:)-closeGapSize*backGap*normalT1(end,:);
    cP2StartB = cP2(1,:)-closeGapSize*backGap*normalT2(1,:);
    cP2EndB = cP2(end,:)-closeGapSize*backGap*normalT2(end,:);
    %     cP1StartB = cP1(1,:);
    %     cP1EndB = cP1(end,:);
    %     cP2StartB = cP2(1,:);
    %     cP2EndB = cP2(end,:);
    
    dMat(1,1) = segments_dist_3d(cP1StartA',cP1StartB',cP2StartA',cP2StartB');
    dMat(1,2) = segments_dist_3d(cP1StartA',cP1StartB',cP2EndA',cP2EndB');
    dMat(2,1) = segments_dist_3d(cP1EndA',cP1EndB',cP2StartA',cP2StartB');
    dMat(2,2) = segments_dist_3d(cP1EndA',cP1EndB',cP2EndA',cP2EndB');
    
    % Compute the angle between the tangents
    angle = abs(acosd(normalT1*normalT2')-180);
    angleMat(1,1) = abs(acosd(normalT1(1,:)*normalT2(1,:)')-180);
    angleMat(1,2) = abs(acosd(normalT1(1,:)*normalT2(2,:)')-180);
    angleMat(2,1) = abs(acosd(normalT1(2,:)*normalT2(1,:)')-180);
    angleMat(2,2) = abs(acosd(normalT1(2,:)*normalT2(2,:)')-180);
    
    % Compute the weights
    weightMat = 1 - (angleMat/angleTol + dMat/maxDist);
    
    weightMat(angleMat>angleTol) = -1;
    weightMat(dMat>maxDist) = -1;
    
    % Find the max weight of all the model start/end combinations
    [colVal rowIdx] = max(weightMat);
    [weightMatMax colIdx] = max(colVal);
    rowIdx = rowIdx(colIdx);
   
%     dis.imaris.displaySegments(cP1StartA,cP1StartB,'Model: Linear Model');
%     dis.imaris.displaySegments(cP1EndA,cP1EndB,'Model: Linear Model');
%     dis.imaris.displaySegments(cP2StartA,cP2StartB,'Model: Linear Model');
%     dis.imaris.displaySegments(cP2EndA,cP2EndB,'Model: Linear Model');
%     dis.imaris.displayPoints(cP1StartA,15,[1.0 0.0 0.0 0.0]);
%     dis.imaris.displayPoints(cP1StartB,30,[0.0 0.0 1.0 0.0]);
%     dis.imaris.displayPoints(cP1EndA,15,[1.0 0.0 0.0 0.0]);
%     dis.imaris.displayPoints(cP1EndB,30,[0.0 0.0 1.0 0.0]);
%     dis.imaris.displayPoints(cP2StartA,15,[1.0 0.0 0.0 0.0]);
%     dis.imaris.displayPoints(cP2StartB,30,[0.0 0.0 1.0 0.0]);
%     dis.imaris.displayPoints(cP2EndA,15,[1.0 0.0 0.0 0.0]);
%     dis.imaris.displayPoints(cP2EndB,30,[0.0 0.0 1.0 0.0]);
    
    if weightMatMax > 0
        cP1End = cP1([1,end],:);
        cP2End = cP2([1,end],:);
        cP1EndLinked = cP1End(rowIdx,:);
        cP2EndLinked = cP2End(colIdx,:);
        
        % Copy color
        % data.clusterColor(clusterIdx2,:) = data.clusterColor(clusterIdx1,:);
 
        linksStart = [linksStart;cP1EndLinked];
        linksEnd = [linksEnd;cP2EndLinked];
        linksEdges = [linksEdges;clusterIdx1,clusterIdx2];
        linksWeights = [linksWeights;weightMatMax]; % Unused so far
        
    end
end

linksStart = linksStart(2:end,:);
linksEnd = linksEnd(2:end,:);
linksEdges = linksEdges(2:end,:);
linksWeights = linksWeights(2:end,:);

%%%% SUBSET SELECTION ###############################################
% Nodes
pnts = [linksStart;linksEnd];
nPnts = size(pnts,1);
nEdges = size(linksStart,1);
[unqPnts,~,n] = unique(pnts,'rows'); % Unique points
nUnqPnts = size(unqPnts,1);
nodes = 1:nUnqPnts;
lut = nodes(n); % Look-up table
edgs = [lut(1:nPnts/2);lut(nPnts/2+1:nPnts)]';
wghts = linksWeights;

% Find the maximum weight matching 
matching = maxWeightedMatching(nPnts,edgs,wghts);

linksStart = linksStart(matching,:);
linksEnd = linksEnd(matching,:);
linksEdges = linksEdges(matching,:);
linksWeights = linksWeights(matching,:);
%%%% SUBSET SELECTION ###############################################

dis.imaris.displaySegments(linksStart(1:end,:),linksEnd(1:end,:),'Post Process: Cluster links');

% Backup cluster color
clusterColorBak = data.clusterColor;

% Set new color
map = linksEdges;
for e=1:size(linksEdges,1)
    parent = map(e,1);
    child = map(e,2);
    map(map==child) = parent; % Parent inherits the color to its child and the children's children
end

% Concatenate map
map = map(:);

% Find unique elements
[linksEdges,idx] = unique(linksEdges(:)); % List of unique node indices
map = map(idx); % The indx of the node color

% Copy color
data.clusterColor(linksEdges,:) = data.clusterColor(map,:);

% Display clusters
% % dis.clusters

% Revert cluster color
% % data.clusterColor = clusterColorBak;

disp('------------------------------------------------');

%%
%%%% REMOVE CLUSTERS SHORTER THAN XXX IF THEY ARE NOT ALIGNED ###############################################
% Find long models
thresholdLength = 700;
isLongCluster = data.modelLength>=thresholdLength;
isShortAndUnlinkedIdx = setdiff(find(~isLongCluster),linksEdges);
isShortAndUnlinked = false(size(isLongCluster));
isShortAndUnlinked(isShortAndUnlinkedIdx) = true;

% Update clusters
shortClusters = data.clusters(isShortAndUnlinked);
data.clusters = data.clusters(~isShortAndUnlinked);

% Update model lists
pro = Processor(data);
pro.removeModels(isShortAndUnlinked);

% Update null cluster
unclustered = horzcat(shortClusters{:})';
data.nullCluster = [data.nullCluster;unclustered];
%%%% REMOVE CLUSTERS SHORTER THAN XXX IF THEY ARE NOT ALIGNED ###############################################
disp('The End')

%%


%% COUNT MODELS
clc;
n = nnz(data.modelLength)
n = numel(data.simModelBezCP)
% l2=data.modelLength;
% numel(data.modelLength)
% data.nPoints
% data.clusters;
% size(data.points,1)
% size(data.clusters)
% s = cellfun(@numel,data.clusters)
% s = s(s>1);
% nelem = sum(s)
% % [l l2]



