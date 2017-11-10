function [] = whDayVisDirectMeasures(dayGeneData,mainDirname,propertyStr)

propOutDir = [mainDirname filesep 'dayGeneControlDirectMeasure' filesep propertyStr filesep];

[maglim,timelim,spacelim] = getDirectLimits(dayGeneData,mainDirname,propertyStr); % get limit (max absolute value) for every measure

nGeneDay = length(dayGeneData);

% Gene vs. pSup
for iGeneDay = 1 : nGeneDay        
    %     curGeneStr = dayGeneData{iGeneDay}.geneStr;
    dayGeneSeqStr = dayGeneData{iGeneDay}.dayGeneSeqStr;
    
    controlDirect = dayGeneData{iGeneDay}.featsControlDirect;
    geneDirect = dayGeneData{iGeneDay}.featsGeneDirect;
    
    outFnameMagTime = [propOutDir dayGeneSeqStr '_' propertyStr '_MagTime.eps'];
    outFnameMagSpace = [propOutDir dayGeneSeqStr '_' propertyStr '_MagSpace.eps'];
        
    if exist(outFnameMagSpace,'file')
        continue;
    end
    
    plotDayGeneDirect(controlDirect(:,1:2),geneDirect(:,1:2),dayGeneSeqStr,[maglim,timelim],outFnameMagTime,'Time Deriv');
    plotDayGeneDirect(controlDirect(:,1:2:3),geneDirect(:,1:2:3),dayGeneSeqStr,[maglim,spacelim],outFnameMagSpace,'Space Deriv');
    close all;
end
end

%%
function [maglim,timelim,spacelim] = getDirectLimits(dayGeneData,mainDirname,propertyStr)

outFname = [mainDirname filesep propertyStr '_direct.mat'];

if exist(outFname,'file')
    load(outFname);
    return;
end

nGeneDay = length(dayGeneData);
mag = zeros(1,nGeneDay);
timee = zeros(1,nGeneDay);
space = zeros(1,nGeneDay);
for iGeneDay = 1 : nGeneDay 
    controlDirect = dayGeneData{iGeneDay}.featsControlDirect;
    geneDirect = dayGeneData{iGeneDay}.featsGeneDirect;
    mag(iGeneDay) = max(abs([controlDirect(:,1)' geneDirect(:,1)']));
    timee(iGeneDay) = max(abs([controlDirect(:,2)' geneDirect(:,2)']));
    space(iGeneDay) = max(abs([controlDirect(:,3)' geneDirect(:,3)']));    
end

maglim = max(mag);
timelim = max(timee);
spacelim = max(space);

save(outFname,'mag','timee','space','maglim','timelim','spacelim');

end

%%
function [] = plotDayGeneDirect(controlDirect,geneDirect,dayGeneSeqStr,directlim,outFnamePCs,yStr)    
fontsize = 24;
h = figure;
xlabel('Magnitude','FontSize',fontsize);
ylabel(yStr,'FontSize',fontsize);
hold on;
        
title(strrep(dayGeneSeqStr,'_','\_'),'FontSize',fontsize);
plot(controlDirect(:,1)',controlDirect(:,2)','o','MarkerEdgeColor','k','LineWidth',2,'MarkerSize',8);
plot(geneDirect(:,1)',geneDirect(:,2)','o','MarkerEdgeColor',[255,165,0]./255,'LineWidth',2,'MarkerSize',8);
xlim([-directlim(1),directlim(1)]);
ylim([-directlim(2),directlim(2)]);
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