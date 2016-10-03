function [outlierIndex, r] = detectVectorFieldOutliersTFM(data,varargin)
% detectVectorFieldOutliersTFM detect and return the outliers in a vector field
%
% Synopsis:        outlierIdx = detectVectorFieldOutliersTFM(data)
%                  [outlierIdx,r] = detectVectorFieldOutliersTFM(data,2)
%
% This function detects outliers within a vectorial field using an extended
% version of the 'median test' of Westerweel et al. 1994, adapted for PTV.
% After finding neighbors within an average bead distance using KDTree,
% the algorithm calculates the directional fluctuation and norm of vector fluctuation
% with respect to the neighborhood median residual for each vertex. 
% A threshold is then applied to this quantity to extract
% the outliers. Directional fluctuation weights more than vector norm
% fluctuation in usual TFM experiment.
%
% Input:
%      data - a vector field, i.e. a matrix of size nx4 where the two first
%      columns give the positions and the two last the displacement
%
%      threshold (optional) - a threshold for the detection criterion.
%      Usually values are between 2-4 depending on the stringency.
%
%      weighted (optional) - a boolean. If true, neighbors influence is
%      weighted using their relative distance to the central point.
%
% Output
%      outlierIndx - the index of the outlier along the first dimension of
%      data
%
%      r - the values of the normalized fluctuation for each element of the
%      vector field
%
% For more information, see:
% J. Westerweel & F. Scarano, Exp. Fluids, 39 1096-1100, 2005.
% J. Duncan et al., Meas. Sci. Technol., 21 057002, 2010.

% Sebastien Besson, Aug 2011
% Sangyoon Han, Mar 2015

% Input check
ip=inputParser;
ip.addRequired('data',@(x) size(x,2)==4);
ip.addOptional('threshold',2,@isscalar);
ip.addOptional('weighted' ,1,@isscalar);
ip.addParamValue('epsilon',.1,@isscalar);
ip.parse(data,varargin{:})
threshold=ip.Results.threshold;
weighted=ip.Results.weighted;
epsilon=ip.Results.epsilon;

% Filter out NaN from the initial data (but keep the index for the
% outliers)
ind=find(~isnan(data(:,3)));
data=data(ind,:);
    
% Take out duplicate points (Sangyoon)
[data,idata,~] = unique(data,'rows'); %data2 = data(idata,:),data = data2(iudata,:)
    
% calculate maximum closest distance
distance=zeros(length(data),1);
% distanceAll=cell(length(data),1); % the second closest distance
distance2=zeros(length(data),1); % the second closest distance
neiBeadsWhole = data(:,1:2);
parfor i=1:length(data)
    neiBeads = neiBeadsWhole;
    neiBeads(i,:)=[];
    [~,distance(i)] = KDTreeClosestPoint(neiBeads,data(i,1:2));
    [~,curDistanceAll] = KDTreeBallQuery(neiBeads,data(i,1:2),5*distance(i));
    if length(curDistanceAll{1})>1
        distance2(i) = curDistanceAll{1}(2);
    else
        distance2(i) = NaN;
    end
end

% Discard vectors that are in sparse location
opts = statset('maxIter', 200);
objDist = cell(3,1);
n = 0;
lastwarn('');
[~, msgidlast] = lastwarn;
warning('off','stats:gmdistribution:FailedToConverge')
warning('off','stats:gmdistribution:MissingData')
fitError = false;
while ~strcmp(msgidlast,'stats:gmdistribution:FailedToConverge') && n<4 && ~fitError
    n = n+1;
    try
        objDist{n} = gmdistribution.fit(distance2, n, 'Options', opts);
        [~, msgidlast] = lastwarn;
    catch
        fitError=true;
    end
end
if n>1
    objDist = objDist(1:n-1);
    [~,idx] = min(cellfun(@(i) i.BIC, objDist));
    objDist = objDist{idx};
    [mu,idx] = sort(objDist.mu);
    svec = squeeze(objDist.Sigma(:,:,idx));
    if length(mu)>1
        threshDist= max(mu(1)+4*svec(1),mu(2)+4*svec(2));
    else
        threshDist= mu+4*svec;
    end
elseif n==1 && fitError
    threshDist = 2*mean(distance2);
end
idxCloseVectors = distance2<threshDist & ~isnan(distance2);
dataFiltered = data(idxCloseVectors,:);
idCloseVectors = find(idxCloseVectors);
idAwayVectors = find(~idxCloseVectors);
neighborhood_distance = 5*max(distance(idxCloseVectors));%quantile(distance,0.95);%mean(distance);%size(refFrame,1)*size(refFrame,2)/length(beads);

% Find neighbors and distances
% [idx,neiDist] = KDTreeBallQuery(data(:,1:2), data(:,1:2), neighborhood_distance);
[idx,neiDist] = KDTreeBallQuery(dataFiltered(:,1:2), dataFiltered(:,1:2), neighborhood_distance);

% Get the triangulation edges and calculate all distances
% edges= tri.edges;
% dp=(data(edges(:,2),1:2)-data(edges(:,1),1:2));
% D= cellfun(@(x) sqrt(sum(x.^2)),neiDist);
% if ~weighted
%     D=ones(size(D));
% end

% Convert the edge list into an adjacency list 
% nodes = unique([edges(:,1)' edges(:,2)']);
% N=cell(numel(nodes),1);
% E=cell(numel(nodes),1);
% for e=1:size(edges,1); 
%     N{edges(e,1)}=[N{edges(e,1)},edges(e,2)]; 
%     E{edges(e,1)}=[E{edges(e,1)},e]; 
%     N{edges(e,2)}=[N{edges(e,2)},edges(e,1)]; 
%     E{edges(e,2)}=[E{edges(e,2)},e];
% end
N = idx;

% Measure weighted local and neighborhood velocities
options = {1:size(dataFiltered,1),'Unif',false};
% d =arrayfun(@(x)D(E{x}),options{:});
d = neiDist;
localVel=arrayfun(@(x)dataFiltered(x,3:4)/(median(d{x})+epsilon),options{:});
neighVel=arrayfun(@(x)dataFiltered(N{x},3:4)./repmat(d{x}+epsilon,1,2),options{:});

% Get median weighted neighborhood velocity
medianVel=cellfun(@median,neighVel,'Unif',false);

% Calculate normalized fluctuation using neighborhood residuals
medianRes=arrayfun(@(x) median(abs(neighVel{x}-repmat(medianVel{x},size(neighVel{x},1),1))),options{:});
normFluct = arrayfun(@(x) abs(localVel{x}-medianVel{x})./(medianRes{x}+epsilon),options{:});
r=cellfun(@norm, normFluct);

% Filter outliers using threshold
outlierIndex = ind(idata([idCloseVectors(r>threshold); idAwayVectors]));