
%% Nearest Neighbor classification with multilpe tsne (temporarily using different subset of features)
% Performed at the well level
% Limited to two labels!

% Input: features & labels
% Output: 

% Assaf Zaritsky, Oct. 2016
% addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition'));
% addpath(genpath('/home2/azaritsky/code/extern/tsne'));

%% Assaf Nov. 2016: excluding date and same cell type from classifciation

function [accracyPerFeature] = pcKnnCls(feats,labels,cellTypeLabels,dates)
close all;

addpath(genpath('/home2/azaritsky/code/extern/tsne'));

inds = ~strcmp(labels,'');
labels = labels(inds);
feats = feats(:,inds);
cellTypeLabels = cellTypeLabels(inds);
dates = dates(inds);

n = length(labels);

% for i = 1 : n
%     train = feats(:,1:i-1),feats(:,i+1:end)
% end

%% todo: remove to use a standard set of features
featInds = {1:80,1:20,1:40,1:60,[1:20, 41:60],[1:20,41:80],[1:20,61:80],21:40,21:60,21:80,[21:40,61:80],41:60,41:80};
% featInds = {1:80,1:20,1:40,1:60};
nf = length(featInds);

% clsRight = nan(1,nf);
% clsWrong = nan(1,nf);
% for i = 1 : nf
%     [clsRight(i),clsWrong(i)] = doTsneNnCls(feats(featInds{i},:),labels);
% end
% 
% accracyPerFeature = clsRight./(clsRight+clsWrong);
% 
% mean(accracyPerFeature)
% 
% retrun;

%

params.init_dims = 15;
params.perplexity = 5;%30;
params.k = 1;

clsRight = nan(1,nf);
clsWrong = nan(1,nf);
for i = 1 : nf
    [clsRight(i),clsWrong(i)] = doTsneNnClsByCellTypeAndDate(feats(featInds{i},:),labels,cellTypeLabels,dates,params);
end
accracyPerFeature = clsRight./(clsRight+clsWrong);
mean(accracyPerFeature)

end

%%
function [clsRight,clsWrong] = doTsneNnClsByCellTypeAndDate(feats,labels,cellTypes,dates,params)

clsSuccess = 0;

n = length(labels);
for i = 1 : n
    relevantInds = excludeCellTypeAndDateInds(i,cellTypes,dates);    
    relevantFeats = [feats(:,i),feats(:,relevantInds)];
    relevantLabels = [labels(i),labels(relevantInds)];
    relevantI = 1;
    
    clsSuccess = clsSuccess + doClassify(relevantI,relevantFeats,relevantLabels,params.init_dims,params.perplexity,params.k);    
end

clsRight = clsSuccess;
clsWrong = n - clsSuccess;

end

%% Return all indices not from same cell time or day
function relevantInds = excludeCellTypeAndDateInds(i,cellTypes,dates)
currType = cellTypes{i};
curDate = dates{i};

relevantInds = find(~(strcmp(cellTypes,currType) | strcmp(dates,curDate)));
end


%% how many samples were classified right/wrong 
function [clsRight,clsWrong] = doTsneNnCls(feats,labels)
init_dims = 15;
perplexity = 5;%30;
k = 1;

clsSuccess = 0;

n = length(labels);
for i = 1 : n
    clsSuccess = clsSuccess + doClassify(i,feats,labels,init_dims,perplexity,k);    
end

clsRight = clsSuccess;
clsWrong = n - clsSuccess;
end

%% classifies wells is, based on the other wells (select randomly to make even number)
function success = doClassify(is,feats,labels,init_dims,perplexity,k)

% success0 = 0;
% success1 = 0;

% nrep = 11;
% 
% for i = 1 : nrep
%     [trainFeats,trainLabels,newI] = selectEqual(feats,labels,is); % _old
%     [trainFeats1,trainLabels1,newI1] = selectEqual_old(feats,labels,is); % _old
%     
%     mappedTrainFeats = tsne(trainFeats', [], 2, init_dims, perplexity);
%     mappedTrainFeats1 = tsne(trainFeats1', [], 2, init_dims, perplexity);
%     
%     k = 1;
%     success0 = success0 + knnTSne(mappedTrainFeats,trainLabels,newI,k);
%     success1 = success1 + knnTSne(mappedTrainFeats1,trainLabels1,newI1,k);
% end
% 
% assert((success0 > nrep/2) == (success1 > nrep/2));
% 
% if success0 > nrep/2 && success1 > nrep/2
%     success = 1;
% else
%     if success0 < nrep/2 && success1 < nrep/2
%         success = 0;
%     end
% end

% [trainFeats,trainLabels,newI] = selectEqual(feats,labels,is); % _old
[trainFeats1,trainLabels1,newI1] = selectEqual_old(feats,labels,is); % _old

% mappedTrainFeats = tsne(trainFeats', [], 2, init_dims, perplexity);
mappedTrainFeats1 = tsne(trainFeats1', [], 2, init_dims, perplexity);

% k = 1;
% success = knnTSne(mappedTrainFeats,trainLabels,newI,k); % _old
success = knnTSne(mappedTrainFeats1,trainLabels1,newI1,k); % _old
end

%% Select equal size groups for classification that include curI
function [trainFeats,trainLabels,newI] = selectEqual(feats,labels,curI)

indsCurI = zeros(1,length(labels));
indsCurI(curI) = 1;

uniqueLabels = unique(labels);
nUniqueLabels = length(uniqueLabels);

assert(nUniqueLabels == 2);

labels1 = strcmp(labels,uniqueLabels{1});
labels2 = strcmp(labels,uniqueLabels{2});

indsLabels1 = find(labels1 & ~indsCurI);
indsLabels2 = find(labels2 & ~indsCurI);
n = min(length(indsLabels1),length(indsLabels2));

subInds1 = randsample(indsLabels1,n);
subInds2 = randsample(indsLabels2,n);
trainFeats = [feats(:,subInds1),feats(:,subInds2),feats(:,curI)];
trainLabels = [labels(subInds1),labels(subInds2),labels(curI)];

newI = 2*n + 1 : length(trainLabels);
end



%% Allow testing multiple indices to test (not included in the training)
% To allow cross-validation without including the same cell type twice
function success = knnTSne(mappedTrainFeats,trainLabels,is,k)
inds4NN = ones(1,length(trainLabels));
inds4NN(is) = 0;

success = 0;

for i = 1 : length(is)
    curI = is(i);
    inds = find(inds4NN);
    distances = pdist2(mappedTrainFeats(curI,:),mappedTrainFeats(inds,:));
    [sortedDist,sortedInds] = sort(distances); % the list here does not include the data at test (at the end of the list)

    knnInds = inds(sortedInds);
    
    kInds = knnInds(1:k);
    kLabels = trainLabels(kInds);
    
    trueLabel = trainLabels{curI};
    
    if sum(strcmp(kLabels,trueLabel)) > k/2
        success = success + 1;
    end
end
end

%%

% Select equal size groups for classification that include curI
function [trainFeats,trainLabels,newI] = selectEqual_old(feats,labels,curI)

uniqueLabels = unique(labels);
nUniqueLabels = length(uniqueLabels);

assert(nUniqueLabels == 2);

% minimal number of measurements
n = min([sum(strcmp(labels,uniqueLabels{1})),sum(strcmp(labels,uniqueLabels{2})),sum(strcmp(labels,labels{curI}))-1]);

trainFeats = [];
trainLabels = {};
for i = 1 : 2
    curLabel = uniqueLabels{i};
    curInds = find(strcmp(labels,curLabel));       
    
    ind = find(curInds==curI); %#ok<EFIND>
    if isempty(ind)
        % not the one under consideration
        subInds = randsample(curInds,n);
        trainFeats = [trainFeats,feats(:,subInds)];
        trainLabels = [trainLabels,labels(subInds)];
    else
        % include i!
        featsNoI = feats(:,curInds([1:ind-1,ind+1:length(curInds)]));
        labelsNoI = labels(:,curInds([1:ind-1,ind+1:length(curInds)]));
        curIndsNoI = curInds([1:ind-1,ind+1:length(curInds)]);
        subInds = randsample(curIndsNoI,n);
        %       trainFeats = [trainFeats,feats(:,[subInds, ind])];
        %       trainLabels = [trainLabels,labels([subInds, ind])];
        trainFeats = [trainFeats,feats(:,[subInds, curInds(ind)])];
        trainLabels = [trainLabels,labels([subInds, curInds(ind)])];
        
        newI = length(trainLabels);
    end
end
end

%%
function sucess = knnTSne_old(mappedTrainFeats,trainLabels,i,k)
distances = pdist2(mappedTrainFeats(i,:),mappedTrainFeats);
[sortedDist,sortedInds] = sort(distances);

kInds = sortedInds(2:1+k);
kLabels = trainLabels(kInds);

if sum(strcmp(kLabels,trainLabels{i})) > k/2
    sucess = 1;
else
    sucess = 0;
end

end


%% old and buggy (?) :-(
% Select equal size groups for classification that include curI
function [trainFeats,trainLabels,newI] = selectEqual_oldNBuggy(feats,labels,curI)

uniqueLabels = unique(labels);
nUniqueLabels = length(uniqueLabels);

assert(nUniqueLabels == 2);

n = min(sum(strcmp(labels,uniqueLabels{1})),sum(strcmp(labels,uniqueLabels{2})));

trainFeats = [];
trainLabels = {};
for i = 1 : 2
    curLabel = uniqueLabels{i};
    curInds = find(strcmp(labels,curLabel));
    
    ind = find(curInds==curI); %#ok<EFIND>
    if isempty(ind)
        % not the one under consideration
        subInds = randsample(curInds,n);
        trainFeats = [trainFeats,feats(:,subInds)];
        trainLabels = [trainLabels,labels(subInds)];
    else
        % include i!
        featsNoI = feats(:,curInds([1:ind-1,ind+1:length(curInds)]));
        labelsNoI = labels(:,curInds([1:ind-1,ind+1:length(curInds)]));
        curIndsNoI = curInds([1:ind-1,ind+1:length(curInds)]);
        subInds = randsample(curIndsNoI,n-1);
        trainFeats = [trainFeats,feats(:,[subInds, ind])];
        trainLabels = [trainLabels,labels([subInds, ind])];
        newI = length(trainLabels);
    end
end
end
