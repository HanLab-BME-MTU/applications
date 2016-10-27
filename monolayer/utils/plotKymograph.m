function [] = plotKymograph(kymograph,params)
ntime = size(kymograph,2);
nspace = size(kymograph,1);
maxTime = ntime * params.timePerFrame;

xTick = 1:(60/params.timePerFrame):((maxTime/params.timePerFrame)+1);
xTickLabel = 0:60:maxTime;
yTick = 1:(50/params.patchSize):((180/(params.patchSize))+1);
yTickLabel = 0:50:180;

h = figure;
colormap('jet');
imagescnan(kymograph);
hold on;
caxis(params.caxis); colorbar;
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[1,maxTime/params.timePerFrame]);
set(haxes,'XTick',xTick);
set(haxes,'XTickLabel',xTickLabel);
set(haxes,'YTick',yTick);
set(haxes,'YTickLabel',yTickLabel);
set(haxes,'FontSize',32);
xlabel('Time (minutes)','FontSize',32); ylabel('Distance from edge (\mum)','FontSize',32);
set(h,'Color','w');
hold off;
export_fig(params.fname);
end