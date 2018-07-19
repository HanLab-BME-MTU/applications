
%% tSNE per well
% Assaf Zaritsky, Oct. 2016
% addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition'));
function [] = pcMetaTSNE(feats,labels,outDname,titleStr)
close all;

addpath(genpath('/home2/azaritsky/code/extern/tsne'));

inds = ~strcmp(labels,'');
labels = labels(inds);
feats = feats(:,inds);

uniqueLabels = unique(labels);
uniqueLabels = uniqueLabels(~strcmp(uniqueLabels,''));

init_dims = 15;
perplexity = 30;

mappedFeats = tsne(feats', [], 2, init_dims, perplexity);

figure;
gscatter(mappedFeats(:,1), mappedFeats(:,2), labels');
export_fig([outDname 'tsne_' titleStr '.eps']);
end