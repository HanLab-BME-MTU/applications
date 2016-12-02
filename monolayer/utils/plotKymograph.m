function [] = plotKymograph(kymograph,params)

if ~isfield(params,'xStepMinutes')
    params.xStepMinutes = 60;
end

if ~isfield(params,'yStepUm')
    params.yStepUm = 50;
end

if ~isfield(params,'yMaxUm')
    params.yMaxUm = 180;
end

if ~isfield(params,'fontsize')
    params.fontsize = 24;
end

% new 2016.11.3: limit space of the kymograph
kymograph = kymograph(1:((params.yMaxUm/(params.patchSize))),:);% +1

ntime = size(kymograph,2);
nspace = size(kymograph,1);
maxTime = ntime * params.timePerFrame;

xTick = 1:(params.xStepMinutes/params.timePerFrame):((maxTime/params.timePerFrame)+1); %xTick = 1:(60/params.timePerFrame):((maxTime/params.timePerFrame)+1);
xTickLabel = 0:params.xStepMinutes:maxTime; %xTickLabel = 0:60:maxTime;
yTick = 1:(params.yStepUm/params.patchSize):((params.yMaxUm/(params.patchSize))+1); %yTick = 1:(50/params.patchSize):((180/(params.patchSize))+1);
yTickLabel = 0:params.yStepUm:params.yMaxUm; %yTickLabel = 0:50:180;

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
set(haxes,'FontSize',params.fontsize);
xlabel('Time (minutes)','FontSize',params.fontsize); ylabel('Distance from edge (\mum)','FontSize',params.fontsize);
set(h,'Color','w');
hold off;
% export_fig(params.fname); % _biohpc did not work (Nov. 2016)
print(params.fname, '-dpdf'); % use MATLAB figure image save for packaging purposes
end