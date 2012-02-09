function updateEdgesAnisotropic(obj,dRef,alpha,samplePeriod,dMaxAlong,dMinAway)

% Maximum distance from the sample points
len = cellfun(@lengthBezier,obj.data.modelBezCP);
dAway = max(dRef - len * alpha,dMinAway);

% Extend and sample models
dAlong = min(dRef + len * alpha,dMaxAlong);
ext = dAlong./len; % Extend the model
nSamples = arrayfun(@ceil,len/samplePeriod);
s = arrayfun(@(a,b,c) linspace(a,b,c)',-ext,1+ext,nSamples,'UniformOutput',false);
t = cellfun(@(a,b) arcLengthToNativeBezierParametrization(a,b,'off'),obj.data.modelBezCP,s,'UniformOutput',false);
samplePnts = cellfun(@renderBezier,obj.data.modelBezCP,t,'UniformOutput',false);
num = cellfun(@(a) size(a,1),samplePnts);
last = cumsum(num);
first = last-num+1;
samplePntsVector = vertcat(samplePnts{:});

% Find the neighbor points of the sample points
dAway = arrayfun(@(a,b) ones(a,1)*b,nSamples,dAway,'UniformOutput',false);
dAwayVector = vertcat(dAway{:});
neighPntIdxVector = KDTreeBallQuery(obj.data.points,samplePntsVector,dAwayVector);
neighPntIdx = arrayfun(@(a,b) vertcat(neighPntIdxVector{a:b}),first,last,'UniformOutput',false);

% Remove duplicate neighbor point indices
neighPntIdx = cellfun(@unique,neighPntIdx,'UniformOutput',false);
neighPntIdxVector = vertcat(neighPntIdx{:});

% Find the parent of the neighbor points
neighModelIdxVector = obj.data.parents(neighPntIdxVector);

% Build the edge list
nNeighbors = cellfun(@numel,neighPntIdx);
modelIdx = arrayfun(@(a,b) ones(a,1)*b,nNeighbors',1:obj.data.nClusters,'UniformOutput',false);
modelIdxVector = vertcat(modelIdx{:});
obj.data.edges = [modelIdxVector neighModelIdxVector];

% Put the smaller cluster index to the left [1 2] and [2 1] => [1 2]
obj.data.edges = sort(obj.data.edges,2);

% Remove the edge duplicates
obj.data.edges = unique(obj.data.edges,'rows');

% Remove self-edges
obj.data.edges = obj.data.edges(obj.data.edges(:,1)~=obj.data.edges(:,2),:);

% Initialize the weights
obj.data.weights = -ones(size(obj.data.edges,1),1);

disp('Process: Edges updated!');
end