%
function [ysAcc, ysAccPatch] = getFlowAccumulation(ROI,dxs,dys,nDistPixels,nParticles,nIterations,patchSize)
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
if isempty(ys) || isempty(xs)
    xparticles = [];
    yparticles = [];    
    warning('BUG! ROI is empty');    
    return;
end
N = length(ys);
inds = round(15 + 1 + (rand(1,nParticles) .* (N-1-30)));
yparticles = ys(inds);
xparticles = xs(inds);
end

%%

function [accParticles] = trackParticles(xparticles,yparticles,dxs,dys,ROI,nIterations)
accParticles = zeros(size(ROI));

if isempty(xparticles)
    return;
end

[sizeY,sizeX] = size(ROI);
dxs(~ROI | isnan(dxs)) = 0;
dys(~ROI | isnan(dys)) = 0;
for i = 1 : nIterations
    xtmp = dxs(sub2ind(size(ROI),yparticles,xparticles));
    ytmp = dys(sub2ind(size(ROI),yparticles,xparticles));
    xparticles = min(max(round(xparticles + xtmp),1),sizeX);
    yparticles = min(max(round(yparticles + ytmp),1),sizeY);
end


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

function [stretchFlow,protNoStretchFlow,noStretchFlow] = getFlowForStretchingAndNonStretchingEvents(flowAccumulationPatch,strechEvents,protrudingCellsNoStretch)
stretchFlow = flowAccumulationPatch(strechEvents);
protNoStretchFlow = flowAccumulationPatch(protrudingCellsNoStretch);
noStretchFlow = flowAccumulationPatch(~strechEvents);
end