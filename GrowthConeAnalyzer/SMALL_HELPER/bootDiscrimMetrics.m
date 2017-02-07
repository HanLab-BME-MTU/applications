function [DMObserved,pValue,plotInfo,CI] = bootDiscrimMetrics(s1, s2, varargin)
% bootDiscrimMetrics:  Quick hack to do some bootstrapping on the inverseDB
%
% Inputs:
%         s1, s2 : 
%      r feature x c neurite matrix 
%      where r is the number of features and c is the number of
%      observations (often cells/movies) 

%
% Options: 
%          alpha : alpha value, default: 0.05.
%
%           nrep : # of permutations, default: 1900. This gives a coefficient of
%                  variation <=0.10 for alpha = 0.05. Calculated as
%                  nrep = (1-alpha)/cv^2/alpha. See [1] for details.
%
% Reference:
% [1] Efron, B. and Tibshirani, R., "An introduction to the Bootstrap," Ch. 15, 1993.
%
%  Example data from [1]
%   y = [10 27 31 40 46 50 52 104 146];
%   z = [16 23 38 94 99 141 197];

% MB 02/07/2017

ip = inputParser;
ip.addRequired('s1', @isnumeric);
ip.addRequired('s2', @isnumeric);
ip.addOptional('alpha', 0.05, @isscalar);
ip.addOptional('nrep', 1900, @isscalar);

ip.addParamValue('CmpFunction', [] );
ip.parse(s1, s2, varargin{:})
nrep = ip.Results.nrep;
if isempty(ip.Results.CmpFunction)
    fct = @inverseDB; 
else
    fct = ip.Results.CmpFunction;
end
% s1 = s1(:);
% s2 = s2(:);
% s1 = s1(~isnan(s1));
% s2 = s2(~isnan(s2)); 

% combine the data : r = features c = neurites 
sAll = [s1  s2];

n1 = size(s1,2);
n2 = size(s2,2);
N = n1+n2;
[DMObserved, plotInfo] = fct(s1,s2);
% Calculate the number of permutations. If small, run exact test
w = warning('off', 'MATLAB:nchoosek:LargeCoefficient');
%nperms = nchoosek(n1+n2, n1);
warning(w);
% if nperms<=nrep % calculate all permutations
%     P = false(N, nperms);
%     pidx = nchoosek(1:N, n1); % returns row index of class 'sample 1'
%     % convert to linear index
%     pidx = pidx + repmat(N*(0:nperms-1)', [1 n1]);
%     % category (1->sample1, 0->sample2) matrix for all permutations
%     P(pidx) = true;
%     DM = zeros(nperms,1);
%     for i = 1:nperms
%         % perform the bootstrp discm metrics -
%         DM(i) = fct(sAll(:,P(:,i)),sAll(:,~P(:,i)));
%     end
%     ns = nperms;
% else % compute 'nrep' random permutations
    DM = zeros(nrep,1);
    for i = 1:nrep
        idx = randperm(N); % calculate random permutation of the samples
        DM(i) = fct(sAll(:,idx(1:n1)), sAll(:,idx(n1+1:end)));
    end
    ns = nrep;
%end
%% calculate confidence intervals 
% label s1 
ID1 = repmat(1,[size(s1,2),1]); 
% label s2 
ID2 = repmat(2,[size(s2,2),1]); 

IDAll = [ID1; ID2]; 
forCI = zeros(nrep,1); 
for i = 1:nrep
    % select N values with replacement keeping the labels. 
    idx2 = randi(N,N,1);
    IDsC = IDAll(idx2);
    sC = sAll(:,idx2);
    forCI(i) = fct(sC(:,IDsC==1),sC(:,IDsC==2));
    
    
end
 
CI(1) = prctile(forCI,ip.Results.alpha/2*100); % lower
CI(2) = prctile(forCI,(1-ip.Results.alpha/2)*100); % upper
%  [mu,s] = fitGaussianModeToPDF(DM); 
%  CI = 3*s; 
% 
% just formulate as the probability of observing a value at least this
% extreme by chance if assume the entire dataset is one population, you are 
% randomly sampling values of N sample from it and calculating the inverse 
% DB statistic 
pValue = sum(DM>=DMObserved)/ns; 

% deltaRef = fct(s1)-fct(s2);
% 
% switch ip.Results.tail
%     case 'both'
%         pValue = sum(abs(delta)>=abs(deltaRef))/ns;
%     case 'right'
%         pValue = sum(delta>=deltaRef)/ns;
%     case 'left'
%         pValue = sum(delta<=deltaRef)/ns;
% end
% 
% H = pValue <= ip.Results.alpha;
end 
function [inverseDB,plotInfo] = inverseDB(featsS1,featsS2)
% featsSample = rxc double array
%      where r is the number of features and c is the number of
%      observations (often cells/movies) of the perturbation condition

meanS1 = nanmean(featsS1,2)'; % get the average of each feature for 
meanS2 = nanmean(featsS2,2)';

d = pdist2(meanS1,meanS2); % get the euclidean distance between these two vectors 

% intra cluster distance to center
intraClustDistS1 = pdist2(featsS1',meanS1);
intraClustDistS2 = pdist2(featsS2',meanS2);



% genesStd = std(genesDist);
% controlStd = std(controlDist);
c1 = nanmean(intraClustDistS1); % genes mean distance to cluster center

c2 = nanmean(intraClustDistS2); % control mean distance to cluster center

% (c1+c2) / d
inverseDB = d / (c1+c2);
plotInfo.c1 = c1; 
plotInfo.c2 = c2; 
plotInfo.d = d; 
plotInfo.intraClustDistS1 = intraClustDistS1; 
plotInfo.intraClustDistS2 = intraClustDistS2; 
plotInfo.meanS1 = meanS1; 
plotInfo.meanS2 = meanS2; 
end