
function [] = visualizePCoefficients2016(cntrlPCA,timePartition,spatialPartition,outDir)
pcPlot(cntrlPCA.speed,'Speed',timePartition,spatialPartition,outDir); close all;
pcPlot(cntrlPCA.directional,'Directionality',timePartition,spatialPartition,outDir); close all;
pcPlot(cntrlPCA.coordination,'Coordination',timePartition,spatialPartition,outDir); close all;
end

function pcPlot(data,measureStr,timePartition,spatialPartition,outDir)
nfeats = timePartition * spatialPartition;
iPc = 1;
while(data.accVariance(iPc) < 0.96)
    h = figure;
    xlabel('Coefficients','FontSize',32);
    ylabel('Weight','FontSize',32);
    title(sprintf('PC #%d: %.2f variance',iPc,data.accVariance(iPc)),'FontSize',32)
    hold all;
    plot(1:nfeats,data.coeff(:,iPc),'sr','MarkerFaceColor','r','MarkerSize',10);
    
    haxes = get(h,'CurrentAxes');
    set(haxes,'XLim',[0,nfeats+1]);
    set(haxes,'XTick',1:3:nfeats);
    set(haxes,'XTickLabel',1:3:nfeats);
    set(haxes,'YLim',[-0.5,0.6]);
    set(haxes,'YTick',-0.5:0.5:0.5);
    set(haxes,'YTickLabel',-0.5:0.5:0.5);        
    set(haxes,'FontSize',32);
    
    hold off;
    
    outFname = [outDir filesep measureStr '_PC_' num2str(iPc) '.eps'];     
    export_fig(outFname);
    
    %     % ----------------- TIME -----------------
    %     h = figure;
    %     xlabel('Space','FontSize',32);
    %     ylabel('Weight','FontSize',32);
    %     title(sprintf('PC #%d: %.2f variance',iPc,data.accVariance(iPc)),'FontSize',32)
    %     hold all;
    %     plot(1:4,data.coeff(1:4,iPc),'s','Color',[0,0.4,0],'MarkerFaceColor',[0,0.4,0],'MarkerSize',10);
    %     plot(1:4,data.coeff(5:8,iPc),'s','Color',[0,0.7,0],'MarkerFaceColor',[0,0.7,0],'MarkerSize',10);
    %     plot(1:4,data.coeff(9:12,iPc),'s','Color',[0,1,0],'MarkerFaceColor',[0,1,0],'MarkerSize',10);
    %     legend('time1','time2','time3','FontSize',32,'Location','NorthEastOutside');
    %
    %
    %     haxes = get(h,'CurrentAxes');
    %     set(haxes,'XLim',[0,5]);
    %     set(haxes,'XTick',1:4);
    %     set(haxes,'XTickLabel',1:4);
    %     set(haxes,'YLim',[-0.5,0.6]);
    %     set(haxes,'YTick',-0.5:0.5:0.5);
    %     set(haxes,'YTickLabel',-0.5:0.5:0.5);
    %     set(haxes,'FontSize',32);
    %
    %     hold off;
    %
    %     pcaFname = [mainDirname 'pcs/' measureStr '_PC_' num2str(iPc) '_Time.bmp'];
    %
    %     eval(sprintf('print -dbmp16m  %s', pcaFname));
    %
    %     % ----------------------------------
    %     h = figure;
    %     xlabel('Time','FontSize',32);
    %     ylabel('Weight','FontSize',32);
    %     title(sprintf('PC #%d: %.2f variance',iPc,data.accVariance(iPc)),'FontSize',32)
    %     hold all;
    %     plot(1:3,data.coeff(1:4:12,iPc),'s','Color',[0,0.4,0],'MarkerFaceColor',[0,0.4,0],'MarkerSize',10);
    %     plot(1:3,data.coeff(2:4:12,iPc),'s','Color',[0,0.6,0],'MarkerFaceColor',[0,0.6,0],'MarkerSize',10);
    %     plot(1:3,data.coeff(3:4:12,iPc),'s','Color',[0,0.8,0],'MarkerFaceColor',[0,0.8,0],'MarkerSize',10);
    %     plot(1:3,data.coeff(4:4:12,iPc),'s','Color',[0,1,0],'MarkerFaceColor',[0,1,0],'MarkerSize',10);
    %     legend('space1','space2','space3','space4','FontSize',32,'Location','NorthEastOutside');
    %
    %     haxes = get(h,'CurrentAxes');
    %     set(haxes,'XLim',[0,4]);
    %     set(haxes,'XTick',1:3);
    %     set(haxes,'XTickLabel',1:3);
    %     set(haxes,'YLim',[-0.5,0.6]);
    %     set(haxes,'YTick',-0.5:0.5:0.5);
    %     set(haxes,'YTickLabel',-0.5:0.5:0.5);
    %     set(haxes,'FontSize',32);
    %
    %     hold off;
    %
    %     pcaFname = [mainDirname 'pcs/' measureStr '_PC_' num2str(iPc) '_Space.bmp'];
    %
    %     eval(sprintf('print -dbmp16m  %s', pcaFname));
    %
    
    iPc = iPc + 1;
end
end