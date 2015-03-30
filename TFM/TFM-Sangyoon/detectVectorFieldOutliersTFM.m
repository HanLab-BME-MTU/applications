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
[data,idata,iudata] = unique(data,'rows'); %data2 = data(idata,:),data = data2(iudata,:)
    
% calculate maximum closest distance
distance=zeros(length(data),1);
for i=1:length(data)
    neiBeads = data(:,1:2);
    neiBeads(i,:)=[];
    [~,distance(i)] = KDTreeClosestPoint(neiBeads,data(i,1:2));
end
neighborhood_distance = 3*max(distance);%quantile(distance,0.95);%mean(distance);%size(refFrame,1)*size(refFrame,2)/length(beads);

% Find neighbors and distances
[idx,neiDist] = KDTreeBallQuery(data(:,1:2), data(:,1:2), neighborhood_distance);

% Get the triangulation edges and calculate all distances
% edges= tri.edges;
% dp=(data(edges(:,2),1:2)-data(edges(:,1),1:2));
D= cellfun(@(x) sqrt(sum(x.^2)),neiDist);
if ~weighted
    D=ones(size(D));
end

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
options = {1:size(data,1),'Unif',false};
% d =arrayfun(@(x)D(E{x}),options{:});
d = neiDist;
localVel=arrayfun(@(x)data(x,3:4)/(median(d{x})+epsilon),options{:});
neighVel=arrayfun(@(x)data(N{x},3:4)./repmat(d{x}+epsilon,1,2),options{:});

% Get median weighted neighborhood velocity
medianVel=cellfun(@median,neighVel,'Unif',false);

% Calculate normalized fluctuation using neighborhood residuals
medianRes=arrayfun(@(x) median(abs(neighVel{x}-repmat(medianVel{x},size(neighVel{x},1),1))),options{:});
normFluct = arrayfun(@(x) abs(localVel{x}-medianVel{x})./(medianRes{x}+epsilon),options{:});
r=cellfun(@norm, normFluct);

% Filter outliers using threshold
outlierIndex = ind(idata(r>threshold));