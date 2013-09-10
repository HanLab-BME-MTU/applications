function syFigureStyle(h,ax,height)
% syFigureStyle(h,height) change the figure size and font for publishable
% format
% Sangyoon Han August 2013
set(gca,'FontSize',7,'LineWidth',0.5);
set(findobj(ax,'Type','line','-and','Marker','o'),'MarkerSize',3)
set(h,'Units','inches')
hSize = get(h, 'Position');
set(h,'PaperPositionMode','auto')
set(h, 'Position', [hSize(1),hSize(2),hSize(3)/hSize(4)*height height])
