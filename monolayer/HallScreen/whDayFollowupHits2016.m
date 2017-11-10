function [] = whDayFollowupHits2016(dayGeneData,mainDirname,propertyStr,validateGenes)

followupDir = [mainDirname filesep 'dayGeneControlFollowup' filesep propertyStr filesep];

if ~exist(followupDir,'dir')
    mkdir(followupDir);
end

nValidateGenes = length(validateGenes);
nGeneDay = length(dayGeneData);

loggerFname = [followupDir 'logFollowups_' propertyStr '.txt'];
logger = fopen(loggerFname,'w');

[pc1lim,pc2lim,pc3lim,healingRateLimit] = getPCLimits(dayGeneData); % get limit (max absolute value) for every pc & max healing rate

for ivalid = 1 : nValidateGenes
    valGeneStr = validateGenes{ivalid};
    
    validGeneDir = [followupDir filesep valGeneStr filesep];
    if ~exist(validGeneDir,'dir')
        mkdir(validGeneDir);
    end
    
    seqStrs = {};
    seqStrsRep = {};
    
    controlMeanPC1 = [];
    controlMeanPC2 = [];
    controlMeanPC3 = [];
    controlHealingRate = [];
    
    geneMeanPC1 = [];
    geneMeanPC2 = [];
    geneMeanPC3 = []; 
    geneHealingRate = [];        
    
    controlRepPC1 = [];
    controlRepPC2 = [];
    controlRepPC3 = [];
    controlRepHealingRate = [];
    
    geneRepPC1 = [];
    geneRepPC2 = [];
    geneRepPC3 = [];
    geneRepHealingRate = [];
    
    outFnamePC1 = [validGeneDir propertyStr '_PC1.eps'];
    outFnamePC2 = [validGeneDir propertyStr '_PC2.eps'];
    outFnamePC3 = [validGeneDir propertyStr '_PC3.eps'];    
    outFnameHealingRate = [validGeneDir propertyStr '_healingRate.eps'];
    
    outFnamePC1rep = [validGeneDir propertyStr '_PC1_rep.eps'];
    outFnamePC2rep = [validGeneDir propertyStr '_PC2_rep.eps'];
    outFnamePC3rep = [validGeneDir propertyStr '_PC3_rep.eps'];
    outFnameHealingRateRep = [validGeneDir propertyStr '_healingRate_rep.eps'];
    
    if exist(outFnameHealingRate,'file')
        continue;
    end
    
    % Gene vs. pSup
    for iGeneDay = 1 : nGeneDay
        curGeneStr = dayGeneData{iGeneDay}.geneStr;
        
        %         if (strcmp(curGeneStr,valGeneStr) && (dayGeneData{iGeneDay}.KD > 50 || dayGeneData{iGeneDay}.KD == -1 || isnan(dayGeneData{iGeneDay}.KD)) && ...
        %                 sum(isnan(dayGeneData{iGeneDay}.meanControl)) == 0 && ...
        %                 sum(isnan(dayGeneData{iGeneDay}.meanControl)) == 0)             %#ok<USENS>
        
        if (strcmp(curGeneStr,valGeneStr) && (dayGeneData{iGeneDay}.KD > 0 || dayGeneData{iGeneDay}.KD == -1 || isnan(dayGeneData{iGeneDay}.KD)) && ...
                sum(isnan(dayGeneData{iGeneDay}.meanControl)) == 0 && ...
                sum(isnan(dayGeneData{iGeneDay}.meoutFnamePC1anControl)) == 0)             %#ok<USENS>
            
            curSeqStr = dayGeneData{iGeneDay}.SeqStr;
            
            seqStrs = [seqStrs,{curSeqStr}];            
            
            controlMeanPC1 = [controlMeanPC1 dayGeneData{iGeneDay}.meanControlPCA(1)];
            controlMeanPC2 = [controlMeanPC2 dayGeneData{iGeneDay}.meanControlPCA(2)];
            controlMeanPC3 = [controlMeanPC3 dayGeneData{iGeneDay}.meanControlPCA(3)];
            
            geneMeanPC1 = [geneMeanPC1 dayGeneData{iGeneDay}.meanGenePCA(1)];
            geneMeanPC2 = [geneMeanPC2 dayGeneData{iGeneDay}.meanGenePCA(2)];
            geneMeanPC3 = [geneMeanPC3 dayGeneData{iGeneDay}.meanGenePCA(3)];
            
            controlHealingRate = [controlHealingRaoutFnamePC1te dayGeneData{iGeneDay}.healingRateControl];
            geneHealingRate = [geneHealingRate dayGeneData{iGeneDay}.healingRateGene];
            
            nRep = size(dayGeneData{iGeneDay}.featsGenePCA,1);
            
            seqStrsRep = [seqStrsRep,repmat({curSeqStr},1,nRep)];
            
            controlRepPC1 = [controlRepPC1 ones(1,nRep)*dayGeneData{iGeneDay}.meanControlPCA(1)];
            controlRepPC2 = [controlRepPC2 ones(1,nRep)*dayGeneData{iGeneDay}.meanControlPCA(2)];
            controlRepPC3 = [controlRepPC3 ones(1,nRep)*dayGeneData{iGeneDay}.meanControlPCA(3)];
            
            geneRepPC1 = [geneRepPC1 dayGeneData{iGeneDay}.featsGenePCA(:,1)'];
            geneRepPC2 = [geneRepPC2 dayGeneData{iGeneDay}.featsGenePCA(:,2)'];
            geneRepPC3 = [geneRepPC3 dayGeneData{iGeneDay}.featsGenePCA(:,3)'];                           
            
            controlRepHealingRate = [controlRepHealingRate ones(1,nRep)*mean(dayGeneData{iGeneDay}.controlHealingRates)];
            geneRepHealingRate = [geneRepHealingRate dayGeneData{iGeneDay}.geneHealingRates];
        end
    end
    
    %% Visualization + statistics
                
    pvalPC1 = signrank(controlMeanPC1,geneMeanPC1);
    pvalPC2 = signrank(controlMeanPC2,geneMeanPC2);
    pvalPC3 = signrank(controlMeanPC3,geneMeanPC3);
    plotDayGenePC(controlMeanPC1,geneMeanPC1,valGeneStr,seqStrs,'PC1',[-pc1lim,pc1lim],pvalPC1,outFnamePC1);
    plotDayGenePC(controlMeanPC2,geneMeanPC2,valGeneStr,seqStrs,'PC2',[-pc2lim,pc2lim],pvalPC2,outFnamePC2);
    plotDayGenePC(controlMeanPC3,geneMeanPC3,valGeneStr,seqStrs,'PC3',[-pc3lim,pc3lim],pvalPC3,outFnamePC3);
    
    assert(min([controlHealingRate,geneHealingRate]) > 3);
    pvalHealingRate = signrank(controlHealingRate,geneHealingRate);
    plotDayGenePC(controlHealingRate,geneHealingRate,valGeneStr,seqStrs,'HealingRate',[3,healingRateLimit],pvalHealingRate,outFnameHealingRate);
             
    pvalRepPC1 = signrank(controlRepPC1,geneRepPC1);
    pvalRepPC2 = signrank(controlRepPC2,geneRepPC2);
    pvalRepPC3 = signrank(controlRepPC3,geneRepPC3);
    plotDayGenePC(controlRepPC1,geneRepPC1,valGeneStr,seqStrsRep,'PC1',[-pc1lim,pc1lim],pvalRepPC1,outFnamePC1rep);
    plotDayGenePC(controlRepPC2,geneRepPC2,valGeneStr,seqStrsRep,'PC2',[-pc2lim,pc2lim],pvalRepPC2,outFnamePC2rep);
    plotDayGenePC(controlRepPC3,geneRepPC3,valGeneStr,seqStrsRep,'PC3',[-pc3lim,pc3lim],pvalRepPC3,outFnamePC3rep);
    
    pvalRepHealingRate = signrank(controlRepHealingRate,geneRepHealingRate);
    plotDayGenePC(controlRepHealingRate,geneRepHealingRate,valGeneStr,seqStrsRep,'HealingRate',[10,healingRateLimit],pvalRepHealingRate,outFnameHealingRateRep);
    
    %% log and save
    
    fprintf(logger,sprintf('\n ****** %s (N = %d, n = %d) ********\n',valGeneStr,length(controlMeanPC1),length(controlRepPC1)));
    fprintf(logger,sprintf('Healing rate: %.1f vs. %.1f, pval %.4f (%s)\n',...
        mean(controlRepHealingRate),mean(geneRepHealingRate),pvalRepHealingRate,getStars(pvalRepHealingRate)));
    fprintf(logger,sprintf('PC1: %.1f vs. %.1f, pval %.4f (%s)\n',...
        mean(controlRepPC1),mean(geneRepPC1),pvalRepPC1,getStars(pvalRepPC1)));
    fprintf(logger,sprintf('PC2: %.1f vs. %.1f, pval %.4f (%s)\n',...
        mean(controlRepPC2),mean(geneRepPC2),pvalRepPC2,getStars(pvalRepPC2)));
    fprintf(logger,sprintf('PC3: %.1f vs. %.1f, pval %.4f (%s)\n',...
        mean(controlRepPC3),mean(geneRepPC3),pvalRepPC3,getStars(pvalRepPC3)));                            
    
    save([validGeneDir propertyStr '_stats.mat'],...
        'controlMeanPC1','geneMeanPC1','controlMeanPC2',...
        'geneMeanPC2','controlMeanPC3','geneMeanPC3',...
        'controlHealingRate','geneHealingRate',...
        'controlRepPC1','geneRepPC1','controlRepPC2',...
        'geneRepPC2','controlRepPC3','geneRepPC3',...
        'controlRepHealingRate','geneRepHealingRate',...
        'pvalPC1','pvalPC2','pvalPC3',...
        'pvalRepPC1','pvalRepPC2','pvalRepPC3',...
        'pvalHealingRate','pvalRepHealingRate');
    
    close all;        
end

fclose(logger);
end

%%
function [pc1lim,pc2lim,pc3lim,healingRateLimit] = getPCLimits(dayGeneData)

% outFname = [mainDirname filesep propertyStr '_pcs.mat'];
% 
% if exist(outFname,'file')
%     load(outFname);
%     return;
% end

nGeneDay = length(dayGeneData);
pc1 = zeros(1,nGeneDay);
pc2 = zeros(1,nGeneDay);
pc3 = zeros(1,nGeneDay);
healingRate = zeros(1,nGeneDay);
for iGeneDay = 1 : nGeneDay 
    controlPCA = dayGeneData{iGeneDay}.featsControlPCA;
    genePCA = dayGeneData{iGeneDay}.featsGenePCA;
    pc1(iGeneDay) = max(abs([controlPCA(:,1)' genePCA(:,1)']));
    pc2(iGeneDay) = max(abs([controlPCA(:,2)' genePCA(:,2)']));
    pc3(iGeneDay) = max(abs([controlPCA(:,3)' genePCA(:,3)'])); 
    healingRate(iGeneDay) = max(dayGeneData{iGeneDay}.healingRateControl,dayGeneData{iGeneDay}.healingRateGene);
end

pc1lim = max(pc1);
pc2lim = max(pc2);
pc3lim = max(pc3);
healingRateLimit = max(healingRate);

% save(outFname,'pc1','pc2','pc3','pc1lim','pc2lim','pc3lim');

end

%%
function [] = plotDayGenePC(controlPCA,genePCA,geneStr,seqStrs,pcStr,pclim,pval,outFnamePC)  

uniqueSeq = unique(seqStrs);
nSeq = length(uniqueSeq);

Markers = {'o','+','*','s','v'};

fontsize = 16;outFnamePC1
h = figure;
xlabel('Control ','FontSize',fontsize);
ylabel('Condition','FontSize',fontsize);
hold on;
        
title(sprintf('%s (%s): p = %.4f',strrep([geneStr],'_','\_'),pcStr,pval),'FontSize',16);
for i = 1 : nSeq
    inds = strcmp(seqStrs,uniqueSeq{i});
    plot(controlPCA(inds),genePCA(inds),Markers{i},'MarkerEdgeColor',[255,165,0]./255,'LineWidth',2,'MarkerSize',8);
end
legend(uniqueSeq,'FontSize',fontsize);
plot(pclim,pclim,'--k','LineWidth',2);
xlim(pclim);
ylim(pclim);
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',fontsize);
set(h,'Color','w');
axis square;
%         position = get(h,'position');
%         set(h,'position',[position(1:2) round(1.2*position(3:4))]);
%
export_fig(outFnamePC);
hold off;
end

function starsStr = getStars(pval)
if pval >= 0.01
    starsStr = '';
else if pval >= 0.001
        starsStr = '*';
    else if pval >= 0.0001
            starsStr = '**';
        else if pval < 0.0001
                starsStr = '***';
            end
        end        
    end
end
end