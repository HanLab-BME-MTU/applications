function [] = gefFig2_pcsHealingRateAssociation()
addpath(genpath('/home2/azaritsky/code/applications/monolayer/Figures/'));
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
close all; clc;

figDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Figures/Fig2/pcsHealingRateAssociation/';
dataDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ScreenFinal/kd50time200/healingRate/control/';

propertiesStr = {'Speed','Directionality','Coordination'};

for ip = 1 : length(propertiesStr)
    propertyStr = propertiesStr{ip};
    prefix = sprintf('healingRateVs%sPC',propertyStr);
    fprintf(sprintf('----- %s ----\n',propertyStr));
    for ipc = 1 : 3
        curFname = [dataDname prefix num2str(ipc) '.mat'];
        outFname = [figDname prefix num2str(ipc) '.eps'];
        load(curFname); % 'healingRate','pcData','pcLabel','maxAbsPC'
        inds = ~isnan(healingRate') & ~isnan(pcData);
        plotPcHealingRateAssociation(healingRate(inds),pcData(inds),pcLabel,maxAbsPC,outFname);
    end
    fprintf(sprintf('%s\n',propertyStr));
end

end


%%
function [] = plotPcHealingRateAssociation(healingRate,pcData,pcLabel,maxAbsPC,outFname)
% [left bottom width height]
FPosition = [0 0 250 250];
APosition = [0.2 0.2 0.75 0.75]; 

fontsize = 10;

h = figure;
xlabel('Healing rate (\mum hour{-1})','FontSize',fontsize);
ylabel(pcLabel,'FontSize',fontsize);
hold all;
plot(healingRate,pcData,'ok','MarkerSize',3,'MarkerFaceColor','k');
ylim([-maxAbsPC,maxAbsPC]);
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',fontsize);
set(h,'Color','w');
set(h,'Position',FPosition,'PaperPositionMode','auto');
axisHandle= findobj(h,'type','axes');
set(axisHandle,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize,'LineWidth',2);
set(get(axisHandle,'XLabel'),'FontSize',fontsize); set(get(axisHandle,'YLabel'),'FontSize',fontsize);
% position = get(h,'position');
% set(h,'position',[position(1:2) round(1.2*position(3:4))]);
hold off;
export_fig(outFname);

[rho,pval] = corr(healingRate',pcData); % 'type','Spearman'
fprintf(sprintf('%s: Rho = %.2f, Pval = %f\n',pcLabel,rho,pval));
end

