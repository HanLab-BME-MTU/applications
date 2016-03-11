function [ h, hq, qh ] = fig3_qq( data , labels, lim)
%fig3_qq Helper function for fig3 qqplots

    import lamins.plot.*;

    h = figure;
%     pvec = 100 - 0.5.^(1:10);
    pvec = 100 - [0.5.^(0:10)*10 0];

    for i=1:6
        hq(i) = subplot(6,1,i);
        qh(i) = qqplot_enhanced(data{i*2+[-1 0]},[0:10:90]);
        hold on;
        scatter_prctile(data{i*2+[-1 0]},pvec,'g','Marker','^');
        xlabel(labels{i*2-1});
        ylabel(labels{i*2});
    %     ylim([0 8]);
        text(0.05,0.9,labels{i*2},'Units','normalized')
        text(0.75,0.1,labels{i*2-1},'Units','normalized')
        text(0.35,0.9,['m = ' num2str(qh(i).scaleFactor,3)],'Units','normalized');
        if(nargin > 2)
            xlim(lim);
            ylim(lim);
        end
    end
    lim.x = [0 max(cellfun(@(x) x(end),arrayfun(@(x) xlim(x),hq,'Unif',false)))];
    lim.y = [0 max(cellfun(@(x) x(end),arrayfun(@(x) ylim(x),hq,'Unif',false)))];
    for i=1:6
%         set(hq(i),'Units','normalized');
%         set(hq(i),'Position',[0.1 1-0.8/6*(i) 0.8 1/6-0.04]);
        xlim(hq(i),lim.x);
        ylim(hq(i),lim.y);
    end
    set(h,'Position',[369 33 290 1064])


end

