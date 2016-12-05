function [] = whDayVisPCA2016(dayGeneData,mainDirname,propertyStr)

propOutDir = [mainDirname filesep 'dayGeneControlPCA' filesep propertyStr filesep];

[pc1lim,pc2lim,pc3lim] = getPCLimits(dayGeneData,mainDirname,propertyStr); % get limit (max absolute value) for every pc

nGeneDay = length(dayGeneData);

% Gene vs. pSup
for iGeneDay = 1 : nGeneDay        
    %     curGeneStr = dayGeneData{iGeneDay}.geneStr;
    dayGeneSeqStr = dayGeneData{iGeneDay}.dayGeneSeqStr;
    
    controlPCA = dayGeneData{iGeneDay}.featsControlPCA;
    genePCA = dayGeneData{iGeneDay}.featsGenePCA;
    
    outFnamePC1PC2 = [propOutDir dayGeneSeqStr '_' propertyStr '_PC1PC2.eps'];
    outFnamePC1PC3 = [propOutDir dayGeneSeqStr '_' propertyStr '_PC1PC3.eps'];
        
    if exist(outFnamePC1PC3,'file')
        continue;
    end
    
    plotDayGenePC(controlPCA(:,1:2),genePCA(:,1:2),dayGeneSeqStr,[pc1lim,pc2lim],outFnamePC1PC2);
    plotDayGenePC(controlPCA(:,1:2:3),genePCA(:,1:2:3),dayGeneSeqStr,[pc1lim,pc3lim],outFnamePC1PC3);
    close all;
end
end

%%
function [pc1lim,pc2lim,pc3lim] = getPCLimits(dayGeneData,mainDirname,propertyStr)

outFname = [mainDirname filesep propertyStr '_pcs.mat'];

if exist(outFname,'file')
    load(outFname);
    return;
end

nGeneDay = length(dayGeneData);
pc1 = zeros(1,nGeneDay);
pc2 = zeros(1,nGeneDay);
pc3 = zeros(1,nGeneDay);
for iGeneDay = 1 : nGeneDay 
    controlPCA = dayGeneData{iGeneDay}.featsControlPCA;
    genePCA = dayGeneData{iGeneDay}.featsGenePCA;
    pc1(iGeneDay) = max(abs([controlPCA(:,1)' genePCA(:,1)']));
    pc2(iGeneDay) = max(abs([controlPCA(:,2)' genePCA(:,2)']));
    pc3(iGeneDay) = max(abs([controlPCA(:,3)' genePCA(:,3)']));    
end

pc1lim = max(pc1);
pc2lim = max(pc2);
pc3lim = max(pc3);

save(outFname,'pc1','pc2','pc3','pc1lim','pc2lim','pc3lim');

end

%%
function [] = plotDayGenePC(controlPCA,genePCA,dayGeneSeqStr,pclim,outFnamePCs)    
fontsize = 24;
h = figure;
xlabel('PC1','FontSize',fontsize);
ylabel(['PC' outFnamePCs(end-4)],'FontSize',fontsize); %#ok<COLND>
hold on;
        
title(strrep(dayGeneSeqStr,'_','\_'),'FontSize',fontsize);
plot(controlPCA(:,1)',controlPCA(:,2)','o','MarkerEdgeColor','k','LineWidth',2,'MarkerSize',8);
plot(genePCA(:,1)',genePCA(:,2)','o','MarkerEdgeColor',[255,165,0]./255,'LineWidth',2,'MarkerSize',8);
xlim([-pclim(1),pclim(1)]);
ylim([-pclim(2),pclim(2)]);
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',fontsize);
set(h,'Color','w');
axis square;
%         position = get(h,'position');
%         set(h,'position',[position(1:2) round(1.2*position(3:4))]);
%
export_fig(outFnamePCs);
hold off;
end