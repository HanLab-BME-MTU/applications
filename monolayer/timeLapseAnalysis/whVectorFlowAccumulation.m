function [pVal] = whVectorFlowAccumulation(params,dirs)

flowFname = [dirs.plithotaxis dirs.expname '_flow'];
plithotaxisFname = [dirs.plithotaxis dirs.expname '_plithotaxis'];

% if exist([flowFname '.mat'],'file')
%     load([flowFname '.mat']); % flowAccumulation
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

for t = time
    load([dirs.mfData pad(t,3) '_mf.mat']); % dxs, dys
    load([dirs.roiData pad(t,3) '_roi.mat']); % ROI
    
    [ysAcc, ysAccPatch] = getAccumulation(ROI,dxs,dys,nDistPixels,nParticles,nIterations,params.patchSize);
    
    flowAccumulation(:,t) = ysAcc;
    flowAccumulationPatch(:,t) = ysAccPatch;
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

load([plithotaxisFname '.mat']);% strainEventsOutput.numEventsYAxis
strechEvents = strainEventsOutput.numEventsYAxis;
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
[rData, rSim, pVal, rsDelay] = strechFlowAssociation(cumsumStrech,cumsumFlow,maxTimeDelay);
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

save([flowFname '.mat'],'flowAccumulation','flowAccumulationPatch',...
    'cumsumFlowOrigRes','cumsumFlow',...
    'cumsumStrechOrigRes','cumsumStrech',...
    'rData', 'rSim', 'pVal','rsDelay','maxR', 'maxDelayInd','maxDelay');

fprintf(sprintf('%s: pval = %.3f, delay = %d minutes (R: %.2f, Rsim: %.2f, Rdelay: %.2f)\n',dirs.expname,pVal,maxDelay,rData,rSim,maxR));

end

%%
function [ysAcc, ysAccPatch] = getAccumulation(ROI,dxs,dys,nDistPixels,nParticles,nIterations,patchSize)
[ySize,xSize] = size(ROI);

% 1. Randomize seed particles x nParticles
[xparticles,yparticles] = randomizeParticles(ROI,nDistPixels,nParticles);

% 2. Track the particles x nDistPixels times
[accParticles] = trackParticles(xparticles,yparticles,dxs,dys,ROI,nIterations);

% 3. Final locations --> % at each point
[ysAcc, ysAccPatch] = quantizeParticles(accParticles,ROI,nParticles,patchSize);

end

%%

function [xparticles,yparticles] = randomizeParticles(ROI,nDistPixels,nParticles)
[ys,xs] = find(bwdist(~ROI) < nDistPixels & ROI);
N = length(ys);
inds = round(15 + 1 + (rand(1,nParticles) .* (N-1-30)));
yparticles = ys(inds);
xparticles = xs(inds);
end

%%

function [accParticles] = trackParticles(xparticles,yparticles,dxs,dys,ROI,nIterations)
dxs(~ROI | isnan(dxs)) = 0;
dys(~ROI | isnan(dys)) = 0;
for i = 1 : nIterations
    xtmp = dxs(sub2ind(size(ROI),yparticles,xparticles));
    ytmp = dys(sub2ind(size(ROI),yparticles,xparticles));
    xparticles = xparticles + xtmp;
    yparticles = yparticles + ytmp;
end
accParticles = zeros(size(ROI));

for i = 1 : length(yparticles)
    accParticles(yparticles(i),xparticles(i)) = accParticles(yparticles(i),xparticles(i)) + 1;
end
end

%%

function [ysAcc,ysAccPatch] = quantizeParticles(accParticles,ROI,nParticles,patchSize)
ySize = size(ROI,1);
step = 20; % pixels
yvals = 15 : step : ySize-15;
nys = length(yvals) - 1;

edgeROI = imdilate(ROI & bwdist(~ROI) < 20,strel('square',17));
accParticles(~edgeROI) = 0;

ysAcc = zeros(1,ySize);

for i = 1 : nys - 1
    ysAcc(yvals(i):yvals(i+1)) = sum(sum(accParticles(yvals(i):yvals(i+1),:)))./nParticles;
end

ysAccPatch = zeros(1,ceil(ySize/patchSize));
for y = 1 : length(ysAccPatch)
    ysAccPatch(y) = sum(ysAcc(((y-1)*patchSize+1):min((y*patchSize),size(ysAcc,2))));
end

end

%% 
function [rData, rSim, pVal,rsDelay] = strechFlowAssociation(cumsumStrech,cumsumFlow,maxTimeDelay)

[istart, iend] = findRelevantIndices(cumsumFlow);
cumsumStrech = cumsumStrech(istart:iend,:);
cumsumFlow = cumsumFlow(istart:iend,:);

[rData, pvalData] = (corr(cumsumStrech(:),cumsumFlow(:)));

rsDelay = getTemporalDelay(cumsumStrech,cumsumFlow,maxTimeDelay);

nIter = 5000;
rs = zeros(1,nIter);
pvals = zeros(1,nIter);
for i = 1 : nIter
    rndFlow = (cumsumFlow(randperm(size(cumsumFlow,1)),:));
    [rs(i), pvals(i)] = corr(cumsumStrech(:),rndFlow(:));
end

rSim = mean(rs);
pVal = sum(rData < rs) / nIter;
end

function [istart, iend] = findRelevantIndices(cumsumFlow)
istart = find(sum(cumsumFlow,2),1,'first');
iend = find(sum(cumsumFlow,2),1,'last');
end

function rsDelay = getTemporalDelay(cumsumStrech,cumsumFlow,maxTimeDelay)
rsDelay = nan(1,2*maxTimeDelay+1);
for lag = -maxTimeDelay : maxTimeDelay
    curStretch = cumsumStrech(:,max(1,1-lag):min(end,end-lag));
    curFlow = cumsumFlow(:,max(1,1+lag):min(end,end+lag));
    
    %     if lag > 0
    %         curStretch = cumsumStrech(:,1:end-lag);
    %         curFlow = cumsumFlow(:,);
    %     else if lag < 0
    %             curStretch = cumsumStrech
    %             curFlow = cumsumFlow
    %         else
    %             curStretch = cumsumStrech;
    %             curFlow = cumsumFlow;
    %         end
    %     end
    rsDelay(lag+maxTimeDelay+1) = corr(curStretch(:),curFlow(:));
end
end