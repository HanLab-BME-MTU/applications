function [p,x0,x1] = whShearStrainCntrlCdc42(metaData,mainDirname)
load([mainDirname '/plithotaxisOut/geneDaySpeed.mat']); % geneDayDiff


nDayGeneSeq = length(geneDayDiff);

params.timePerFrame = metaData.timePerFrame;
params.patchSize = 15;


%% Get CDC42 experiments and controls
controlInds = [];
cdc42Inds = [];
allControlInds = [];

for iGeneDay = 1 : nDayGeneSeq  
    geneStr = geneDayDiff{iGeneDay}.geneStr;
    if strcmp(geneStr,'CDC42')        
        controlInds = [controlInds find(geneDayDiff{iGeneDay}.controlInds)];
        cdc42Inds = [cdc42Inds find(geneDayDiff{iGeneDay}.geneInds)];
    end
    allControlInds = [allControlInds controlInds find(geneDayDiff{iGeneDay}.controlInds)];
end
  
%% make sure pixel size is 1.24
controlInds1 = [];
cdc42Inds1 = [];
allControlInds1 = [];
for i = 1 : length(controlInds)
    if abs(metaData.pixelSize{controlInds(i)} - 1.24) < 0.01
        controlInds1 = [controlInds1 controlInds(i)];
    end
end
        
for i = 1 : length(cdc42Inds)
    if abs(metaData.pixelSize{cdc42Inds(i)} - 1.24) < 0.01
        cdc42Inds1 = [cdc42Inds1 cdc42Inds(i)];
    end
end

for i = 1 : length(allControlInds)
    if abs(metaData.pixelSize{allControlInds(i)} - 1.24) < 0.01
        allControlInds1 = [allControlInds1 allControlInds(i)];
    end
end

%% get # motion events, # shear-strain events
% nControl = length(controlInds1);
% nCdc42 = length(cdc42Inds1);

controlInds1 = unique(controlInds1);
cdc42Inds1 = unique(cdc42Inds1);
allControlInds1 = unique(allControlInds1);

[nMotionControl,nShearStrainControl] = getMotionShearStrainEvents(metaData,mainDirname,controlInds1);
[nMotionCdc42,nShearStrainCdc42] = getMotionShearStrainEvents(metaData,mainDirname,cdc42Inds1);
[nMotionAllControl,nShearStrainAllControl] = getMotionShearStrainEvents(metaData,mainDirname,allControlInds1);

%% Visualize
fontsize = 10;
LineWidth = 1;
markerSize = 4;
FPosition = [0 0 300 200];
APosition = [0.2 0.2 0.75 0.75];
%%  Visualize all controls + linear fit
close all;
[rho,pval] = corr(nMotionAllControl',nShearStrainAllControl');
p = polyfit(nMotionAllControl(nMotionAllControl<3000),nShearStrainAllControl(nMotionAllControl<3000),1);

h = figure;
xlabel('Motion','FontSize',fontsize);
ylabel('Shear strain','FontSize',fontsize);
hold on;
plot(nMotionAllControl,nShearStrainAllControl,'o','MarkerEdgeColor','k','LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','Control');   
x0 = min(nMotionAllControl); x1 = max(nMotionAllControl);
plot([x0,x1],[polyval(p,x0),polyval(p,x1)],'--r','LineWidth',2);
haxes = get(h,'CurrentAxes');
set(h,'Color','w','Position',FPosition,'PaperPositionMode','auto');
set(haxes,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize);
set(haxes,'FontSize',fontsize);
hold off;
outFname = [mainDirname '/plithotaxisOut/nShearStrainAllControls.eps'];
export_fig(outFname);


%% Visualize all CDC42
h = figure;
xlabel('Motion','FontSize',fontsize);
ylabel('Shear strain','FontSize',fontsize);
hold on;
plot(nMotionControl,nShearStrainControl,'o','MarkerEdgeColor','k','LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','Control');   
plot(nMotionCdc42,nShearStrainCdc42,'o','MarkerEdgeColor','b','LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','CDC42');   

haxes = get(h,'CurrentAxes');
set(h,'Color','w','Position',[0 0 300 150],'PaperPositionMode','auto');
set(haxes,'Position',[0.2 0.25 0.75 0.75],'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize);
set(haxes,'FontSize',fontsize);

% legend('show','FontSize',fontsize);

% legend('show','Location','NorthEastOutside');

% plot([xlim(1),xlim(2)],[xlim(1),xlim(2)],'-k','LineWidth',2);

% outFname = [fnamePrefix '_legend.eps'];
% export_fig(outFname);
% 
% legend off;

outFname = [mainDirname '/plithotaxisOut/nShearStrainControlCdc42.eps'];
export_fig(outFname);

plot([x0,x1],[polyval(p,x0),polyval(p,x1)],'--r','LineWidth',2);
set(h,'Color','w','Position',FPosition,'PaperPositionMode','auto');
set(haxes,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize);
set(haxes,'FontSize',fontsize);
outFname = [mainDirname '/plithotaxisOut/nShearStrainControlCdc42_cntrl.eps'];
export_fig(outFname);

hold off;

end

%% 
function [nMotion,nShearStrain] = getMotionShearStrainEvents(metaData,mainDirname,inds)
plithotaxisDir = [mainDirname '/plithotaxis/'];

nMotion = [];
nShearStrain = [];

nInds = length(inds);
for i = 1 : nInds    
    load([plithotaxisDir metaData.fnames{(inds(i))} '_plithotaxis.mat']);
    nMotion = [nMotion sum(strainEventsOutput.nMotionEvents)];
    nShearStrain = [nShearStrain sum(strainEventsOutput.nStrainEvents)]; 
end
end