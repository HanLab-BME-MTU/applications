function [] = whMotionShearStrainGeneDaySeq(metaData,mainDirname,pControl,x0Control,x1Control)
load([mainDirname '/plithotaxisOut/geneDaySpeed.mat']); % geneDayDiff


outdir = [mainDirname '/plithotaxisOut/GeneDaySeq/'];
if ~exist(outdir,'dir')
    mkdir(outdir);
end

nDayGeneSeq = length(geneDayDiff);

% params.timePerFrame = metaData.timePerFrame;
% params.patchSize = 15;

fontsize = 10;
LineWidth = 2;
markerSize = 8;
FPosition = [0 0 300 200];
APosition = [0.2 0.2 0.75 0.75];

%% 

for iGeneDay = 1 : nDayGeneSeq  
    close all;
    if abs(geneDayDiff{iGeneDay}.pixelSize - 1.24) > 0.01
        continue;
    end

    geneSeqDayStr = [geneDayDiff{iGeneDay}.geneStr '_' geneDayDiff{iGeneDay}.SeqStr '_' geneDayDiff{iGeneDay}.dayStr];
    controlInds = find(geneDayDiff{iGeneDay}.controlInds);
    geneInds = find(geneDayDiff{iGeneDay}.geneInds);    

    
    %% get # motion events, # shear-strain events
    [nMotionControl,nShearStrainControl] = getMotionShearStrainEvents(metaData,mainDirname,controlInds);
    [nMotionGene,nShearStrainGene] = getMotionShearStrainEvents(metaData,mainDirname,geneInds);
    
    %% Visualize gene
    h = figure;
    xlabel('Motion','FontSize',fontsize);
    ylabel('Shear strain','FontSize',fontsize);
    hold on;
    plot(nMotionControl,nShearStrainControl,'o','MarkerEdgeColor','k','LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','Control');
    plot(nMotionGene,nShearStrainGene,'o','MarkerEdgeColor','b','LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','CDC42');
    plot([x0Control,x1Control],[polyval(pControl,x0Control),polyval(pControl,x1Control)],'--r','LineWidth',2);
    
    haxes = get(h,'CurrentAxes');
    set(h,'Color','w','Position',FPosition,'PaperPositionMode','auto');
    set(haxes,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize);
    set(haxes,'FontSize',fontsize);
    
    hold off;
    outFname = [outdir 'motionShearStrain_' geneSeqDayStr '.eps'];
    export_fig(outFname);
end

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