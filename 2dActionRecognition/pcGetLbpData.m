%% Returns a list of size 10 x N, were N is the total accumulated number of cells in all experiments in all time points
function [lbpData, lbpLabels, out] = pcGetLbpData(allExps, labels, uniqueLabelsStr, uniqueLabelsValues)

nExps = length(allExps);

lbpData = [];
lbpLabels = [];
for e = 1 : nExps
    curExp = allExps{e};
    curLabel = labels(e); % this is not the index in the unqiue labels, but only an ID for the label - do not get confused!
    
    fname = [curExp.dir '/results/' curExp.name '_lbpPerFrame.mat'];
    
    if ~exist(fname,'file')
        continue;
    end
    
    load(fname); % accumulatedLBP 10 x n
    n = size(accumulatedLBP,2);
    
    lbpData = [lbpData, accumulatedLBP];
    lbpLabels = [lbpLabels; curLabel.*ones(n,1)]; 
end

params.percentile = 2;
params.nBins = 20;
[out] = pcaScatterQuantization(lbpData,lbpLabels,uniqueLabelsStr,uniqueLabelsValues,params); % out.maps, diffMaps

% [coeff,score,latent] = pca(lbpData'); % amount of information cumsum(latent)./sum(latent)

% figure;
% xlabel('First Principal Component','FontSize',24);
% ylabel('Second Principal Component','FontSize',24);
% hold on;
% plot(score(lbpLabels==1,1),score(lbpLabels==1,2),'or','MarkerFaceColor','r','MarkerSize',2);
% plot(score(lbpLabels==2,1),score(lbpLabels==2,2),'o','MarkerFaceColor','c','MarkerSize',2);
% % plot(score(lbpLabels==1,1),score(lbpLabels==1,2),'or','MarkerFaceColor','r','MarkerSize',2);
% set(gca,'XTickLabel',[]);
% set(gca,'YTickLabel',[]);
% hold off;


end