function [] = whVariance(features,geneDayDiff,indsPSup,mainDirname,propertyStr)
close all;
nGeneDaySeq = length(geneDayDiff);

varControl = [];
varGene = [];
dayGenesSeqStr = {};

loggerFname = [mainDirname 'variance/log_' propertyStr '.txt'];
logger = fopen(loggerFname,'w');

fontsize = 24;

h = figure;
xlabel('Control variance','FontSize',fontsize);
ylabel('Gene variance','FontSize',fontsize);
hold on;

% variance of all control experiment & of all experiments

featsAllControls = features(:,indsPSup);

meanAll = mean(features,2)';
meanAllControls = mean(featsAllControls,2)';

distsAll = pdist2(features',meanAll);
distsAllControl = pdist2(featsAllControls',meanAllControls);

varInterDayPSup = var(distsAllControl);
generalVariance = var(distsAll);

fprintf(logger,sprintf('General variace = %.2f \nControls variance = %.2f\n\n',...
    generalVariance,varInterDayPSup));    

% ---------------------

curExp = 0;
for iGeneDay = 1 : nGeneDaySeq      
    nGene = sum(geneDayDiff{iGeneDay}.geneInds);
    nControl = sum(geneDayDiff{iGeneDay}.controlInds);
    
    if nGene < 3 || nControl < 3
        fprintf(logger,sprintf('Excluding %s (gene = %d, control = %d)\n',geneDayDiff{iGeneDay}.geneStr,nGene,nControl));    
        continue;
    end
    
    curExp = curExp + 1;
    
    dayGenesSeqStr{curExp} = geneDayDiff{iGeneDay}.dayGeneSeqStr;
    
    featsGEF = features(:,geneDayDiff{iGeneDay}.geneInds);
    featsControl = features(:,geneDayDiff{iGeneDay}.controlInds);
    
    meanGene = mean(featsGEF,2)';
    meanControl = mean(featsControl,2)';
    
    distsGEF = pdist2(featsGEF',meanGene);
    distsControl = pdist2(featsControl',meanControl);
        
    varControl = [varControl var(distsControl)];
    varGene = [varGene var(distsGEF)];

    plot(varControl(curExp),varGene(curExp),'o','MarkerEdgeColor','k','LineWidth',2,'MarkerSize',7);       
end

fprintf(logger,sprintf('Mean control variance = %.2f\nStd control variance = %.2f\n\n',mean(varControl),std(varControl)));    
fprintf(logger,sprintf('Mean GEFs variance = %.2f\nStd GEFs variance = %.2f\n\n',mean(varGene),std(varGene)));    


haxes = get(h,'CurrentAxes');
% if strcmp(propertyStr,'Speed')
%     maxVar = 2100;
%     assert(max([varControl varGene]) < maxVar);
%     set(haxes,'XLim',[0,maxVar]);
%     set(haxes,'XTick',0:500:maxVar);
%     set(haxes,'XTickLabel',0:500:maxVar);
%     set(haxes,'YLim',[0,maxVar]);
%     set(haxes,'YTick',0:500:maxVar);
%     set(haxes,'YTickLabel',0:500:maxVar);
% else if strcmp(propertyStr,'Directionality')
%         maxVar = 13;
%         assert(max([varControl varGene]) < maxVar);
%         set(haxes,'XLim',[0,maxVar]);
%         set(haxes,'XTick',0:4:maxVar);
%         set(haxes,'XTickLabel',0:4:maxVar);
%         set(haxes,'YLim',[0,maxVar]);
%         set(haxes,'YTick',0:4:maxVar);
%         set(haxes,'YTickLabel',0:4:maxVar);
%     else if strcmp(propertyStr,'Coordination')
%             maxVar = 0.3;
%             assert(max([varControl varGene]) < maxVar);
%             set(haxes,'XLim',[0,maxVar]);
%             set(haxes,'XTick',0:0.1:maxVar);
%             set(haxes,'XTickLabel',0:0.1:maxVar);
%             set(haxes,'YLim',[0,maxVar]);
%             set(haxes,'YTick',0:0.1:maxVar);
%             set(haxes,'YTickLabel',0:0.1:maxVar);
%         end
%     end
% end


set(haxes,'FontSize',fontsize);

maxVar = max([xlim,ylim]);
% step = maxVar/4;
% if step < 1
%     step = round(step*1e2) / 1e2;
% else if step < 2
%         step = round(step*1e1) / 1e1;
%     else
%         step = round(step);
%     end
% end

set(haxes,'XLim',[0,maxVar]);
% set(haxes,'XTick',0:0.1:maxVar);
% set(haxes,'XTickLabel',0:0.1:maxVar);
% set(haxes,'YLim',[0,maxVar]);
% set(haxes,'YTick',0:0.1:maxVar);
% set(haxes,'YTickLabel',0:0.1:maxVar);

set(h,'Color','none');


plot([0,maxVar],[0,maxVar],'--k','LineWidth',3);
plot([varInterDayPSup,varInterDayPSup],[0,maxVar],'--g','LineWidth',3);
plot([0,maxVar],[generalVariance,generalVariance],'--r','LineWidth',3);
plot([generalVariance,generalVariance],[0,maxVar],'--r','LineWidth',3);

controlGeneVarFname = [mainDirname 'variance/controlGeneDayVariance' propertyStr '.eps'];
export_fig(controlGeneVarFname);

hold off;

%% Print sorted variance
[sortedVarControl,sortedIndsControl] = sort(varControl);
sortedDayGeneSeqsStrControl = dayGenesSeqStr(sortedIndsControl);
fprintf(logger,sprintf('\n ****** Variance Controls ********\n'));
for iGeneDay = 1 : length(sortedDayGeneSeqsStrControl)
    fprintf(logger,sprintf('%s: %.2f\n',sortedDayGeneSeqsStrControl{iGeneDay},sortedVarControl(iGeneDay)));
end

[sortedVarGene,sortedIndsGene] = sort(varGene);
sortedDayGeneSeqsStrControl = dayGenesSeqStr(sortedIndsControl);
fprintf(logger,sprintf('\n ****** Variance GEFs ********\n'));
for iGeneDay = 1 : length(sortedDayGeneSeqsStrControl)
    fprintf(logger,sprintf('%s: %.2f\n',sortedDayGeneSeqsStrControl{iGeneDay},sortedVarGene(iGeneDay)));
end

fclose(logger);

end