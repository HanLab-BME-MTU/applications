function measureSpindles(idlist,dicMaxProj)
%MEASURESPINDLES groups tags from several cells and measures spindle length
%
% SYNOPSIS: measureSpindles
%
% INPUT 
%
% OUTPUT 
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 28-Sep-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0 || isempty(idlist)
% load idlist
[idlistName, idlistDir] = uigetfile;
id = load(fullfile(idlistDir,idlistName));
fn = fieldnames(id);
idlist = id.(fn{1});
else
    idlistDir = pwd;
end

if nargin < 2 || isempty(dicMaxProj)
    % load DIC image. For the moment, just don't plot
    doPlot = 0;
else
    doPlot = 1;
end
    
% calculate distances
coords = idlist(1).linklist(:,9:11);
distList = pdist(coords); % default is euclidean

% cluster with complete linkage
links = linkage(distList,'complete');
% make clusters if the longest distance between two points is less than 2
% microns
figure,dendrogram(links,0,'Colorthreshold',2);
labels = cluster(links,'cutoff',2,'criterion','distance');

% for every "cell": find longest distance, number of points, angle opposite
% the longest side, flag whether the intensities of the two tags associated
% with the longest side are more similar than the other tag(s)
nCells = max(labels);
cellData = zeros(nCells, 4);

for c = 1:nCells
    cIdx = find(labels == c);
    nTags = length(cIdx);
    dist = pdist(coords(cIdx,:));
    [maxDist,maxDistIdx] = max(dist);
    
    
    if nTags > 2
        % find which tags the maxDist belongs to
        [r,c]=find(tril(ones(4),-1));
        maxDistCidx = cIdx([r(maxDistIdx);c(maxDistIdx)]);
        minDistCidx = setdiff(cIdx,maxDistCidx);
        
        maxAmp = idlist(1).linklist(maxDistCidx,8);
        deltaMaxAmp = diff(maxAmp);
        minAmp = idlist(1).linklist(minDistCidx,8);
        
        % for each minDistCidx: Find angle, similarAmplitude
        for i=length(minDistCidx)-1:1
            % angle
            v1 = coords(maxDistCidx(1),:)-coords(minDistCidx(i),:);
            v2 = coords(maxDistCidx(2),:)-coords(minDistCidx(i),:);
            angle(i) = 180/pi * acos(dot(v1,v2)/(norm(v1)*norm(v2)));
            
            % amplitude
            similarAmplitude =  deltaMaxAmp < max(maxAmp-minAmp(i));
        end
        
    else
        angle = 180;
        similarAmplitude = 1;        
    end
    
    cellData(c,:) = [nTags,max(dist),min(angle),min(similarAmplitude)];
    
    angle = [];
    similarAmplitude = [];
    
end

% display statistics
metaIdxL = cellData(:,1) == 4 | (cellData(:,1) == 3 & cellData(:,4) == 1);
proMetaIdxL = cellData(:,1) == 3 & ~metaIdxL;
g1IdxL = cellData(:,1) <= 2;
nM3 = sum(cellData(metaIdxL,1) == 3);
nM4 = sum(cellData(metaIdxL,1) == 4);
nPro = sum(proMetaIdxL);
nG1 = sum(g1IdxL);
nM = sum(metaIdxL);

figure,boxplot(cellData(metaIdxL,2),cellData(metaIdxL,1),'notch','on')

disp(sprintf('%1.2f%% G1, %1.2f%% proM, %1.2f%% M, %1.2f%% 3spM, %1.2f%% 4spM, Ntot %i, NM %i',...
    nG1/nCells*100,nPro/nCells*100,nM/nCells*100,nM3/nM*100,nM4/nM*100,nCells,nM))

if doPlot
    figure,imshow(dicMaxProj',[]);
    hold on
    for c = 1:nCells
        cIdx = find(labels == c);
        xy = coords(cIdx,[2,1]);
        plot(xy(:,1)/0.0663,xy(:,2)/0.066,'o','MarkerSize',4,'Color',extendedColors(c));
    end
end