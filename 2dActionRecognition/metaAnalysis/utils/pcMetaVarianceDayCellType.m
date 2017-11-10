%% At the well level: EXPLAIN HERE
% Assaf Zaritsky, May. 2016
% addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition'));
% function pval = pcMetaVarianceDayCellType(featsAll,dateAll,cellTypeAll,outDnameScale)
function pval = pcMetaVarianceDayCellType(D,dateAll,cellTypeAll,outDnameScale)

addpath(genpath('/home2/azaritsky/code/extern'));

close all;

% Dvec = pdist(featsAll','cityblock');
% D = squareform(Dvec);

uniqueDays = unique(dateAll);
nUniqueDays = length(uniqueDays);

allDaySimilarity = [];
allSelfSimilarity = [];

for iDay = 1 : nUniqueDays
    dayInd = find(strcmp(dateAll,uniqueDays{iDay}));
    %     dayInd = strfind([dateAll{:}],uniqueDays{iDay}); % WRONG!
    if length(dayInd) == 1
        continue;
    end
    
    if(length(dayInd) > 2)
        assert(strcmp('140729',uniqueDays{iDay}));
        dayInd = dayInd(2:3);
    end
    
    daySimilarity = D(dayInd(1),dayInd(2));    
    
    selfSimilarity1 = pcMetaCellTypeSelfSimilarity(D,cellTypeAll,dayInd(1)); % self-similarity of 1st cell from that day (to same cells from other days)
    selfSimilarity2 = pcMetaCellTypeSelfSimilarity(D,cellTypeAll,dayInd(2)); % self-similarity of 2nd cell from that day (to same cells from other days)
    selfSimilarity = [selfSimilarity1,selfSimilarity2];
    
    allDaySimilarity = [allDaySimilarity,ones(1,length(selfSimilarity))*daySimilarity]; %#ok<AGROW>
    allSelfSimilarity = [allSelfSimilarity, selfSimilarity];     %#ok<AGROW>
end

pval = plotVarianceDayCellType(allDaySimilarity,allSelfSimilarity,[outDnameScale filesep 'dayCellTypeSimilarity.eps']);

end

%%

% (controlPCA,genePCA,geneStr,pcStr,pclim,pval,outFnamePC)
function pval = plotVarianceDayCellType(allDaySimilarity,allSelfSimilarity,outFname)

pval = ranksum(allDaySimilarity,allSelfSimilarity);

minSS = min([allDaySimilarity,allSelfSimilarity]);
maxSS = max([allDaySimilarity,allSelfSimilarity]);
xylim = [minSS,maxSS];

fontsize = 24;
h = figure;
xlabel('Day','FontSize',fontsize);
ylabel('Cell type','FontSize',fontsize);
hold on;
        
title(sprintf('p = %.4f',pval),'FontSize',fontsize);
plot(allDaySimilarity,allSelfSimilarity,'o','MarkerEdgeColor',[255,165,0]./255,'LineWidth',2,'MarkerSize',8);
plot(xylim,xylim,'--k','LineWidth',2);
xlim(xylim);
ylim(xylim);
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',fontsize);
set(h,'Color','w');
axis square;
axis equal;
%         position = get(h,'position');
%         set(h,'position',[position(1:2) round(1.2*position(3:4))]);
%
export_fig(outFname);
hold off;
end