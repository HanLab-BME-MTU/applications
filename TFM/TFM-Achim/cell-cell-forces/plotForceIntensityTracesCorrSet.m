function []=plotForceIntensityTracesCorrSet(corrSets,doPrint)
% plots the time traces for the cell-cell force and total intensity at the
% same edge. To bring the time traces to comarable scales, the values of
% both time traces are normalized by the respective means.
%
% INPUT     corrSets needs to have the data structure that is obtained by
%           the output of: [corrSets]=collectEdgeValues(groupedNetworks,goodEdgeSet,'corr');
if nargin<2 || isempty(doPrint)
    doPrint=0;
end
for iEdge=1:length(corrSets)
    figure(iEdge)
        plot(corrSets(iEdge).frames,corrSets(iEdge).fcMag/nanmean(corrSets(iEdge).fcMag),'-k','LineWidth',2)
        hold on;
        plot(corrSets(iEdge).frames,corrSets(iEdge).Itot/nanmean(corrSets(iEdge).Itot),'-r','LineWidth',2)
        title(['normalized Cell-Cell Force Mag. and normalized tot. Ecad-Int. of set= ',num2str(iEdge),' of Cluster= ',num2str(corrSets(iEdge).clusterId),' at edge= ',num2str(corrSets(iEdge).edgeId)])
        xlabel(['frame [',num2str(round(corrSets(iEdge).dt_mean)),' sec]'])
        ylabel('Cell-Cell Force Mag./Mean or tot. Ecad-Intensity/Mean')
        box on
        set(gca,'LineWidth',2,'FontSize',20)
        hold off;
        if doPrint
            currFileName=['corrFI_set_',num2str(iEdge),'_clusterId_',num2str(corrSets(iEdge).clusterId),'_edgeId_',num2str(corrSets(iEdge).edgeId)];
            saveas(gcf,[currFileName,'.fig']);
            saveas(gcf,[currFileName,'.eps'],'psc2');
        end
end