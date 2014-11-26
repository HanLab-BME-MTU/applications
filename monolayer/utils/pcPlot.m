function pcPlot(coeff,latent,measureStr,dirname)

fontsize = 24;

accVariance = cumsum(latent)./sum(latent);

iPc = 1;

while(accVariance(iPc) < 0.96)
    h = figure;
    xlabel('Coefficients','FontSize',32);
    ylabel('Weight','FontSize',32);
    title(sprintf('PC #%d: %.2f variance',iPc,accVariance(iPc)),'FontSize',fontsize)
    hold all;
    plot(1:12,coeff(:,iPc),'sk','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',10);    
    haxes = get(h,'CurrentAxes');
    set(haxes,'YLim',[-0.6,0.6]);
    set(haxes,'YTick',-0.6:0.3:0.6);
    set(haxes,'YTickLabel',-0.6:0.3:0.6);
    set(haxes,'XLim',[0,13]);
    set(haxes,'XTick',1:4:12);
    set(haxes,'XTickLabel',1:4:12);
    set(haxes,'FontSize',fontsize);    
    set(h,'Color','none');    
    hold off;    
    pcaFname = [dirname measureStr '_PC_' num2str(iPc) '.eps'];    
    export_fig(pcaFname);        
    
    iPc = iPc + 1;
end
end