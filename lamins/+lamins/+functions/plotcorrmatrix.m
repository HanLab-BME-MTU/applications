function [fig, ax, h, hcb] = plotcorrmatrix(corrmatrix,datasetname)
    fig = figure;
    % layout correlation matrix of images to sort by channel then by z
    total = sqrt(numel(corrmatrix));
    perChannel = size(corrmatrix,2);
    % permute dimensions so that correlations are grouped by channels first
    twoD = reshape(permute(corrmatrix,[2 1 4 3]),total,total);
    %h = imagesc(zeros(80,80));
    axis image;
    imagesc(twoD);
    % set so that x and y are the same size
    axis image;
    % Center 0 correlation
    set(gca,'CLim',[-1 1]);
    set(gca,'XTick',0.5:perChannel:total);
    set(gca,'XGrid','on');
    set(gca,'XTickLabel','');
    set(gca,'YTick',0.5:perChannel:total);
    set(gca,'YGrid','on');
    set(gca,'YTickLabel','');
    ax = gca;
    hcb = colorbar;
    pos = get(gca,'Position');
    h = axes('Position',pos,'Color','none');
    axis image;
    axis ij;
    set(h,'XLim',[0.5 total + 0.5]);
    set(h,'YLim',[0.5 total + 0.5]);
    set(h,'XTick',perChannel/2:perChannel:total);
    set(h,'YTick',perChannel/2:perChannel:total);
    set(h,'XTickLabel',{'Red','Green'})
    set(h,'YTickLabel',{'Red','Green'})
%     set(h,'XTickLabel',{'DAPI','Lamin A','Lamin B1','anti-Lamin B1, B2'})
%     set(h,'YTickLabel',{'DAPI','Lamin A','Lamin B1','anti-Lamin B1, B2'})
    xlabel('Channel, z position');
    ylabel('Channel, z position');
    %set(get(hcb,'YLabel'),'String','2D Correlation');
    ylabel(hcb,'2D Correlation');
    set(h,'Position',get(ax,'Position'));
    set(fig,'ResizeFcn',@(x,y) set(h,'Position',get(ax,'Position')));
    %title('2D Correlation Matrix: ag\_080712wt\_Reconstructed 3');
    %title('2D Correlation Matrix: ag\_072612\_wt\_Reconstructed 2.mat');
    title(['2D Correlation Matrix: ' strrep(datasetname,'_','\_')]);
    
    filename = ['corrmatrix_' strrep(strrep(datasetname,' ','_'),'.mat','')];
    
    %saveas(fig,[filename '.fig'],'fig')
    %saveas(fig,[filename '.eps'],'eps')
end