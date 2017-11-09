function plotSphericalProjection(ax,azimuths,elevations,colorIndex,markerSize,cmap, varargin)
if(~isempty(azimuths))
[X,Y] = pol2cart(azimuths,elevations); 
% m=min(colorIndex(:));M=max(colorIndex(:));
% colormapIdx=1;
% if(m<M)
% colormapIdx=max(min(ceil((colorIndex-m)*length(cmap)/(M-m)),length(cmap)),1);
% end
colors=cmap(colorIndex,:);
scatter(ax,X,Y,markerSize,colors,varargin{:});

end
axis(ax,'square')
xlim(ax,[-1.6 1.6]);
ylim(ax,[-1.6 1.6]);

% % % axes=cell(1,maxTime);
% % % for t=1:maxTime;
% % %     Idx=(times>t)&(times<=t);
% % %     h=polar(ax,azimuths(Idx),elevations(Idx),'+');
% % %     hold on;
% % %     set(h, 'MarkerFaceColor', cmap(ceil(t*length(cmap)/maxTime),:), 'markeredgecolor',cmap(ceil(t*length(cmap)/maxTime),:));
% % % end
% % % hold off