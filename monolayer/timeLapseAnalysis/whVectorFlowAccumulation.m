function [pVal] = whVectorFlowAccumulation(params,dirs)

fprintf('starting vector flow accumulation\n');

flowFname = [dirs.plithotaxis dirs.expname '_flow'];
plithotaxisFname = [dirs.plithotaxis dirs.expname '_plithotaxis'];

% if exist([flowFname '.mat'],'file')
%     load([flowFname '.mat']); % flowAccumulation
%     return;
% end

% time = 1 : params.nTime - params.frameJump - 1; % -1 because segmentation is based on motion
time = 1 : 91 - params.frameJump - 1; % -1 because segmentation is based on motion
ntime = length(time);

load([dirs.roiData pad(1,3) '_roi.mat']); % ROI
[ySize,xSize] = size(ROI);

flowAccumulation = zeros(ySize,ntime);
flowAccumulationPatch = zeros(ceil(ySize/params.patchSize),ntime);

nDistPixels = round(180 ./ params.pixelSize);
nParticles = 500;
nIterations = 500;

load([plithotaxisFname '.mat']);% strainEventsOutput.numEventsYAxis
strainEventsOrigResolution = strainEventsOutput.strainEventsOrigResolution;

flowStretchDirname = [dirs.dirname filesep 'flowStretchVisualization'];

if ~exist([flowStretchDirname '_old'],'dir') && exist(flowStretchDirname,'dir')
    movefile(flowStretchDirname,[flowStretchDirname '_old']);
end

if ~exist(flowStretchDirname,'dir')
    mkdir(flowStretchDirname);
end

if ~exist([flowFname '.mat'],'file') || params.always
    for t = time
        load([dirs.mfData sprintf('%.3d',t) '_mf.mat']); % dxs, dys
        load([dirs.roiData sprintf('%.3d',t) '_roi.mat']); % ROI
        I = imread([dirs.images sprintf('%.3d',t) '.tif']); % images data
        
        [ysAcc, ysAccPatch] = getFlowAccumulation(ROI,dxs,dys,nDistPixels,nParticles,nIterations,params.patchSize);
        
        flowAccumulation(:,t) = ysAcc;
        flowAccumulationPatch(:,t) = ysAccPatch;
        
        %% Visualize motion flows and stretching events!
        [qh] = visualizeMotionFields(I,ROI,dxs,dys,params.patchSize,0.5);% last parameter is reduced resolution
        hold on;
        curEvents = strainEventsOrigResolution{t};
        for e = 1: length(curEvents)
            curY = min(curEvents{e}.ys(1),ySize);
            assert(curEvents{e}.ys(1) <= ySize + 5 && curEvents{e}.ys(1) > 1);
            curX = max(find(ROI(curY,:)));
            plot(curX,curY,'*g','MarkerFaceColor','g','MarkerSize',12,'LineWidth',1);
        end
        hold off;
        export_fig_biohpc([flowStretchDirname filesep sprintf('%.3d',t) '_flowStretchVis.eps']);
        close all;
        
        patchSize = params.patchSize;
        save([flowStretchDirname filesep sprintf('%.3d',t) '_flowStretchVis.mat'],'I','ROI','dxs','dys','patchSize','curEvents','ySize');
    end
else
    load([flowFname '.mat']);
end
ySize = size(flowAccumulation,1);

maxTime = size(flowAccumulation,2) * params.timePerFrame;
xTick = 1:(200/params.timePerFrame):((maxTime/params.timePerFrame)+1);
xTickLabel = 0:200:maxTime;
yTick = 1:200:ySize;
yTickLabel = (1:200:ySize)-1;

h = figure;
imagesc(flowAccumulation);
hold on;
caxis([0.04 0.12]);
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[1,maxTime/params.timePerFrame]);
set(haxes,'XTick',xTick);
set(haxes,'XTickLabel',xTickLabel);
set(haxes,'YTick',yTick);
set(haxes,'YTickLabel',yTickLabel);
set(haxes,'FontSize',32);
xlabel('Time (minutes)','FontSize',32); ylabel('Y-axis','FontSize',32);
set(h,'Color','none');
hold off;
eval(sprintf('print -dbmp16m  %s', [flowFname '.bmp']));

flowEvents = sum(flowAccumulation,2);
cumsumFlowOrigRes = cumsum(flowAccumulation,2);
cumsumFlow = cumsum(flowAccumulationPatch,2);

% strechEvents = strainEventsOutput.numEventsYAxis;
strechEvents = strainEventsOutput.seedsVis == 2;
protrudingCellsNoStretch = strainEventsOutput.seedsVis == 1;
cumsumStrech = strainEventsOutput.cumsumVis;
cumsumStrechOrigRes = imresize(cumsumStrech,size(cumsumFlowOrigRes));

h = figure;
imagesc(cumsumStrechOrigRes);
hold on;
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[1,maxTime/params.timePerFrame]);
set(haxes,'XTick',xTick);
set(haxes,'XTickLabel',xTickLabel);
set(haxes,'YTick',yTick);
set(haxes,'YTickLabel',yTickLabel);
set(haxes,'FontSize',32);
xlabel('Time (minutes)','FontSize',32); ylabel('Y-axis','FontSize',32);
set(h,'Color','none');
hold off;
eval(sprintf('print -dbmp16m  %s', [dirs.plithotaxis dirs.expname '_accStretch.bmp']));

h = figure;
imagesc(cumsumFlowOrigRes);
hold on;
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[1,maxTime/params.timePerFrame]);
set(haxes,'XTick',xTick);
set(haxes,'XTickLabel',xTickLabel);
set(haxes,'YTick',yTick);
set(haxes,'YTickLabel',yTickLabel);
set(haxes,'FontSize',32);
xlabel('Time (minutes)','FontSize',32); ylabel('Y-axis','FontSize',32);
set(h,'Color','none');
hold off;
eval(sprintf('print -dbmp16m  %s', [dirs.plithotaxis dirs.expname '_accFlow.bmp']));

maxTimeDelay = 20; % delay in time (frames)
% TRYING TO REDUCE Y
[rData, rSim, pVal, rsDelay] = strechFlowAssociation(imresize(cumsumStrech,[size(cumsumStrech,1)/2,size(cumsumStrech,2)]),...
    imresize(cumsumFlow,[size(cumsumStrech,1)/2,size(cumsumStrech,2)]),maxTimeDelay);
% [rData, rSim, pVal, rsDelay] = strechFlowAssociation(cumsumStrech,cumsumFlow,maxTimeDelay);
[maxR, maxDelayInd] = max(rsDelay);
maxDelay = (maxDelayInd-maxTimeDelay-1) * params.timePerFrame;

h = figure;
hold on;
plot(params.timePerFrame*(-maxTimeDelay:maxTimeDelay),rsDelay,'ok','MarkerSize',5,'MarkerFaceColor','k');
plot(maxDelay,maxR,'or','MarkerSize',10,'MarkerFaceColor','r');
plot(params.timePerFrame*(-maxTimeDelay:maxTimeDelay),ones(1,2*maxTimeDelay+1).*rSim,'g--','LineWidth',3);
haxes = get(h,'CurrentAxes');
% set(haxes,'XLim',params.timePerFrame.*[-maxTimeDelay,maxTimeDelay]);
% set(haxes,'XTick',-maxTimeDelay:maxTimeDelay/2:maxTimeDelay);
% set(haxes,'XTickLabel',-maxTimeDelay:maxTimeDelay/2:maxTimeDelay);
% set(haxes,'YLim',[0,1]);
% set(haxes,'YTick',0:0.5:1);
% set(haxes,'YTickLabel',0:0.5:1);
set(haxes,'FontSize',32);
plot([maxDelay,maxDelay],get(haxes,'YLim'),'b--','LineWidth',3);
xlabel('Time (minutes)','FontSize',32); ylabel('RHO','FontSize',32);
% set(h,'Color','none');
hold off;
eval(sprintf('print -dbmp16m  %s', [dirs.plithotaxis dirs.expname '_accStretchFlowTimeLag.bmp']));


%% Boxplots for flow in different stretching resolution - (not very informative)
% [stretchFlow,protNoStretchFlow,noStretchFlow] = getFlowForStretchingAndNonStretchingEvents(flowAccumulationPatch,strechEvents,protrudingCellsNoStretch);
%
% n1 = length(stretchFlow);
% n2 = length(protNoStretchFlow);
% n3 = length(noStretchFlow);
%
% % Plithotaxis
% data = [stretchFlow;protNoStretchFlow;noStretchFlow];
%
% labels = [...
%     repmat('Stretching   ',n1,1);repmat('Protruding   ',n2,1);repmat('No stretching',n3,1)];
% colors = 'rck';
% fontsize = 24;
% h = figure;
% hold on;
% boxplot(data,labels,'whisker',0.7193); % 0.7193 --> 90 % of the data  {'Speed      ','Speed      ','Stress Mag.','Stress Mag.','Plithotaxis','Plithotaxis','Anisotropy ','Anisotropy '}
% haxes = get(h,'CurrentAxes');
% % set(haxes,'YLim',[-1.1,5]);
% % set(haxes,'YTick',-1:2:5);
% % set(haxes,'YTickLabel',-1:2:5);
% set(gca,'FontSize',fontsize);
% text_h = findobj(gca, 'Type', 'text');
% for cnt = 1:length(text_h)
%     set(text_h(cnt),'FontSize', fontsize,'VerticalAlignment', 'top','HorizontalAlignment', 'center');%,'Rotation',45
% end
% % rotateticklabel(haxes,45);
% h1 = findobj(gca,'Tag','Box');
% h2 = findobj(gca,'Tag','Upper Whisker');
% h3 = findobj(gca,'Tag','Lower Whisker');
% h4 = findobj(gca,'Tag','Median');
% h5 = findobj(gca,'Tag','Upper Adjacent Value');
% h6 = findobj(gca,'Tag','Lower Adjacent Value');
% for j=1:length(h1)
%     patch(get(h1(j),'XData'),get(h1(j),'YData'),colors(j),'FaceAlpha',.5,'LineWidth',2);
%     patch(get(h2(j),'XData'),get(h2(j),'YData'),colors(j),'FaceAlpha',.5,'LineWidth',2);
%     patch(get(h3(j),'XData'),get(h3(j),'YData'),colors(j),'FaceAlpha',.5,'LineWidth',2);
%     patch(get(h4(j),'XData'),get(h4(j),'YData'),colors(j),'FaceAlpha',.5,'LineWidth',2);
%     patch(get(h5(j),'XData'),get(h5(j),'YData'),colors(j),'FaceAlpha',.5,'LineWidth',2);
%     patch(get(h6(j),'XData'),get(h6(j),'YData'),colors(j),'FaceAlpha',.5,'LineWidth',2);
% end
% ylabel('Flow probability');
% oh=findobj(gca,'tag','Outliers');
% set(oh,'Visible','off');
% set(h,'Color','none');
%
% % position = get(h,'position');
% % set(h,'position',[position(1:2) round(1.5*position(3:4))]);
%
% % set(gcf,'Position',[100 100 500 500]) %this changes the size of the figure window
% % set(gca,'Position',[0.1 0.1 0.70 0.7]) %this changes the area in between the figure window and the plot axis
% hold off;
%
% export_fig_biohpc([dirs.plithotaxis dirs.expname '_stretchFlowBoxplot.eps']);



%%
save([flowFname '.mat'],'flowAccumulation','flowAccumulationPatch',...
    'cumsumFlowOrigRes','cumsumFlow',...
    'cumsumStrechOrigRes','cumsumStrech',...
    'rData', 'rSim', 'pVal','rsDelay','maxR', 'maxDelayInd','maxDelay');

fprintf(sprintf('%s: pval = %.3f, delay = %d minutes (R: %.2f, Rsim: %.2f, Rdelay: %.2f)\n',dirs.expname,pVal,maxDelay,rData,rSim,maxR));

end

%% Revision: moved to getFlowAccumulation file
% % %%
% % function [ysAcc, ysAccPatch] = getAccumulation(ROI,dxs,dys,nDistPixels,nParticles,nIterations,patchSize)
% % [ySize,xSize] = size(ROI);
% %
% % % 1. Randomize seed particles x nParticles
% % [xparticles,yparticles] = randomizeParticles(ROI,nDistPixels,nParticles);
% %
% % % 2. Track the particles x nDistPixels times
% % [accParticles] = trackParticles(xparticles,yparticles,dxs,dys,ROI,nIterations);
% %
% % % 3. Final locations --> % at each point
% % [ysAcc, ysAccPatch] = quantizeParticles(accParticles,ROI,nParticles,patchSize);
% %
% % end
% %
% % %%
% %
% % function [xparticles,yparticles] = randomizeParticles(ROI,nDistPixels,nParticles)
% % [ys,xs] = find(bwdist(~ROI) < nDistPixels & ROI);
% % if isempty(ys) || isempty(xs)
% %     error('BUG! ROI is empty');
% % end
% % N = length(ys);
% % inds = round(15 + 1 + (rand(1,nParticles) .* (N-1-30)));
% % yparticles = ys(inds);
% % xparticles = xs(inds);
% % end
% %
% % %%
% %
% % function [accParticles] = trackParticles(xparticles,yparticles,dxs,dys,ROI,nIterations)
% % [sizeY,sizeX] = size(ROI);
% % dxs(~ROI | isnan(dxs)) = 0;
% % dys(~ROI | isnan(dys)) = 0;
% % for i = 1 : nIterations
% %     xtmp = dxs(sub2ind(size(ROI),yparticles,xparticles));
% %     ytmp = dys(sub2ind(size(ROI),yparticles,xparticles));
% %     xparticles = min(max(round(xparticles + xtmp),1),sizeX);
% %     yparticles = min(max(round(yparticles + ytmp),1),sizeY);
% % end
% % accParticles = zeros(size(ROI));
% %
% % for i = 1 : length(yparticles)
% %     accParticles(yparticles(i),xparticles(i)) = accParticles(yparticles(i),xparticles(i)) + 1;
% % end
% % end
% %
% % %%
% %
% % function [ysAcc,ysAccPatch] = quantizeParticles(accParticles,ROI,nParticles,patchSize)
% % ySize = size(ROI,1);
% % step = 20; % pixels
% % yvals = 15 : step : ySize-15;
% % nys = length(yvals) - 1;
% %
% % edgeROI = imdilate(ROI & bwdist(~ROI) < 20,strel('square',17));
% % accParticles(~edgeROI) = 0;
% %
% % ysAcc = zeros(1,ySize);
% %
% % for i = 1 : nys - 1
% %     ysAcc(yvals(i):yvals(i+1)) = sum(sum(accParticles(yvals(i):yvals(i+1),:)))./nParticles;
% % end
% %
% % ysAccPatch = zeros(1,ceil(ySize/patchSize));
% % for y = 1 : length(ysAccPatch)
% %     ysAccPatch(y) = sum(ysAcc(((y-1)*patchSize+1):min((y*patchSize),size(ysAcc,2))));
% % end
% %
% % end

%% Revision: moved to seperate function with the same name strechFlowAssociation
% % %%
% % function [rData, rSim, pVal,rsDelay] = strechFlowAssociation(cumsumStrech,cumsumFlow,maxTimeDelay)
% %
% % [istart, iend] = findRelevantIndices(cumsumFlow);
% % cumsumStrech = cumsumStrech(istart:iend,:);
% % cumsumFlow = cumsumFlow(istart:iend,:);
% %
% % [rData, pvalData] = (corr(cumsumStrech(:),cumsumFlow(:)));
% %
% % rsDelay = getTemporalDelay(cumsumStrech,cumsumFlow,maxTimeDelay);
% %
% % nIter = 5000;
% % rs = zeros(1,nIter);
% % pvals = zeros(1,nIter);
% % for i = 1 : nIter
% %     rndFlow = (cumsumFlow(randperm(size(cumsumFlow,1)),:));
% %     [rs(i), pvals(i)] = corr(cumsumStrech(:),rndFlow(:));
% % end
% %
% % rSim = mean(rs);
% % pVal = sum(rData < rs) / nIter;
% % end
% %
% % function [istart, iend] = findRelevantIndices(cumsumFlow)
% % istart = find(sum(cumsumFlow,2),1,'first');
% % iend = find(sum(cumsumFlow,2),1,'last');
% % end
% %
% % function rsDelay = getTemporalDelay(cumsumStrech,cumsumFlow,maxTimeDelay)
% % rsDelay = nan(1,2*maxTimeDelay+1);
% % for lag = -maxTimeDelay : maxTimeDelay
% %     curStretch = cumsumStrech(:,max(1,1-lag):min(end,end-lag));
% %     curFlow = cumsumFlow(:,max(1,1+lag):min(end,end+lag));
% %
% %     %     if lag > 0
% %     %         curStretch = cumsumStrech(:,1:end-lag);
% %     %         curFlow = cumsumFlow(:,);
% %     %     else if lag < 0
% %     %             curStretch = cumsumStrech
% %     %             curFlow = cumsumFlow
% %     %         else
% %     %             curStretch = cumsumStrech;
% %     %             curFlow = cumsumFlow;
% %     %         end
% %     %     end
% %     rsDelay(lag+maxTimeDelay+1) = corr(curStretch(:),curFlow(:));
% % end
% % end

% % function [stretchFlow,protNoStretchFlow,noStretchFlow] = getFlowForStretchingAndNonStretchingEvents(flowAccumulationPatch,strechEvents,protrudingCellsNoStretch)
% % stretchFlow = flowAccumulationPatch(strechEvents);
% % protNoStretchFlow = flowAccumulationPatch(protrudingCellsNoStretch);
% % noStretchFlow = flowAccumulationPatch(~strechEvents);
% % end