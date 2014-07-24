function testTrendSignificance(particleBehavior,...
    windowDistFromEdge,mode2plot,minNP,figureName,saveLoc,convFact,...
    yAxisLabel,axisLimits,plotWinDist,winFigName,winFigLoc,compArea)

if nargin < 4 || isempty(minNP)
    minNP = 20;
end

if nargin < 5 || isempty(figureName)
    figureName = 'test';
end

if nargin < 6 || isempty(saveLoc)
    saveLoc = [];
end

if nargin < 7 || isempty(convFact)
    convFact = [0.111 10 1]; %conversion factors for pixel size, frame time, property
end
convFactDist = convFact(1);
convFactTime = convFact(2);
convFactProp = convFact(3);

if nargin < 8 || isempty(yAxisLabel)
    yAxisLabel = 'Property';
end

if nargin < 9 || isempty(axisLimits)
    axisLimits = [];
end

if nargin < 10 || isempty(plotWinDist)
    plotWinDist = 0;
end

if nargin < 11 || isempty(winFigName)
    winFigName = 'edge & windows';
end

if nargin < 12 || isempty(winFigLoc)
    winFigLoc = [];
end

if nargin < 13 || isempty(compArea)
    compArea = 0;
end

%get increment range
minInc = -size(particleBehavior.befStatic.mean,1);
maxInc = size(particleBehavior.aftStatic.mean,1);

%get number of bands and modes
[~,numBand,numMode] = size(particleBehavior.onset.mean);

%construct series of behavior to plot

%static series
staticSeriesMean = [particleBehavior.befStatic.mean(end:-1:1,:,:); ...
    particleBehavior.onset.mean; particleBehavior.aftStatic.mean] * convFactProp;
staticSeriesStd = [particleBehavior.befStatic.std(end:-1:1,:,:); ...
    particleBehavior.onset.std; particleBehavior.aftStatic.std] * convFactProp;
staticSeriesNP = [particleBehavior.befStatic.numPoints(end:-1:1,:,:); ...
    particleBehavior.onset.numPoints; particleBehavior.aftStatic.numPoints];
staticSeriesNP(isnan(staticSeriesNP)) = 1;

staticSeriesMean = staticSeriesMean(:,:,mode2plot);
staticSeriesStd = staticSeriesStd(:,:,mode2plot);
staticSeriesNP = staticSeriesNP(:,:,mode2plot);

tmp = [windowDistFromEdge.befStatic(end:-1:1,:,:,1); ...
    windowDistFromEdge.onset(:,:,:,1); windowDistFromEdge.aftStatic(:,:,:,1)] ...
    * convFactDist;
staticSeriesDistMid = (tmp(:,:,1) + tmp(:,:,2)) / 2;
staticSeriesDistHR = (tmp(:,:,2) - tmp(:,:,1)) / 2;

staticSeriesBulkDistMid = (tmp(:,4,1) + tmp(:,end,2)) / 2;
staticSeriesBulkDistHR = (tmp(:,end,2) - tmp(:,4,1)) / 2;

%dynamic series before
dynamicSeriesBefMean = [particleBehavior.befDynamic.mean(end:-1:1,:,:); ...
    particleBehavior.onset.mean; NaN(maxInc,numBand,numMode)] * convFactProp;
dynamicSeriesBefStd = [particleBehavior.befDynamic.std(end:-1:1,:,:); ...
    particleBehavior.onset.std; NaN(maxInc,numBand,numMode)] * convFactProp;
dynamicSeriesBefNP = [particleBehavior.befDynamic.numPoints(end:-1:1,:,:); ...
    particleBehavior.onset.numPoints; ones(maxInc,numBand,numMode)];
dynamicSeriesBefNP(isnan(dynamicSeriesBefNP)) = 1;

dynamicSeriesBefMean = dynamicSeriesBefMean(:,:,mode2plot);
dynamicSeriesBefStd = dynamicSeriesBefStd(:,:,mode2plot);
dynamicSeriesBefNP = dynamicSeriesBefNP(:,:,mode2plot);

tmp = [windowDistFromEdge.befDynamic(end:-1:1,:,:,1); ...
    windowDistFromEdge.onset(:,:,:,1); NaN(maxInc,numBand,2)] ...
    * convFactDist;
dynamicSeriesBefDistMid = (tmp(:,:,1) + tmp(:,:,2)) / 2;
dynamicSeriesBefDistHR = (tmp(:,:,2) - tmp(:,:,1)) / 2;

%dynamic series after
dynamicSeriesAftMean = [NaN(-minInc+1,maxInc,numMode); particleBehavior.aftDynamic.mean] ...
    * convFactProp;
dynamicSeriesAftStd = [NaN(-minInc+1,maxInc,numMode); particleBehavior.aftDynamic.std] ...
    * convFactProp;
dynamicSeriesAftNP = [ones(-minInc+1,maxInc,numMode); particleBehavior.aftDynamic.numPoints];
dynamicSeriesAftNP(isnan(dynamicSeriesAftNP)) = 1;

dynamicSeriesAftMean = dynamicSeriesAftMean(:,:,mode2plot);
dynamicSeriesAftStd = dynamicSeriesAftStd(:,:,mode2plot);
dynamicSeriesAftNP = dynamicSeriesAftNP(:,:,mode2plot);

tmp = [NaN(-minInc+1,maxInc,2); windowDistFromEdge.aftDynamic(:,:,:,1)] ...
    * convFactDist;
dynamicSeriesAftDistMid = (tmp(:,:,1) + tmp(:,:,2)) / 2;
dynamicSeriesAftDistHR = (tmp(:,:,2) - tmp(:,:,1)) / 2;

%dynamic series right behind edge before and after
edgeSeriesMean = [dynamicSeriesBefMean(1:-minInc+1,1); diag(dynamicSeriesAftMean(-minInc+2:end,:))];
edgeSeriesStd = [dynamicSeriesBefStd(1:-minInc+1,1); diag(dynamicSeriesAftStd(-minInc+2:end,:))];
edgeSeriesNP = [dynamicSeriesBefNP(1:-minInc+1,1); diag(dynamicSeriesAftNP(-minInc+2:end,:))];

edgeSeriesDistMid = [dynamicSeriesBefDistMid(1:-minInc+1,1); diag(dynamicSeriesAftDistMid(-minInc+2:end,:))];
edgeSeriesDistHR = [dynamicSeriesBefDistHR(1:-minInc+1,1); diag(dynamicSeriesAftDistHR(-minInc+2:end,:))];

%combined dynamic series after
combDynamicSeriesAftMean = [NaN(-minInc,1,numMode); particleBehavior.aftDynamicComb.mean; ...
    NaN(1,1,numMode)] * convFactProp;
combDynamicSeriesAftStd = [NaN(-minInc,1,numMode); particleBehavior.aftDynamicComb.std; ...
    NaN(1,1,numMode)] * convFactProp;
combDynamicSeriesAftNP = [ones(-minInc,1,numMode); particleBehavior.aftDynamicComb.numPoints; ones(1,1,numMode)];
combDynamicSeriesAftNP(isnan(combDynamicSeriesAftNP)) = 1;

combDynamicSeriesAftMean = combDynamicSeriesAftMean(:,:,mode2plot);
combDynamicSeriesAftStd = combDynamicSeriesAftStd(:,:,mode2plot);
combDynamicSeriesAftNP = combDynamicSeriesAftNP(:,:,mode2plot);

tmp = [NaN(-minInc,1,2); windowDistFromEdge.aftDynamicComb(:,:,:,1); ...
    NaN(1,1,2)] * convFactDist;
combDynamicSeriesAftDistMid = (tmp(:,:,1) + tmp(:,:,2)) / 2;
combDynamicSeriesAftDistHR = (tmp(:,:,2) - tmp(:,:,1)) / 2;

%remove measurements with less than minNP points
staticSeriesMean(staticSeriesNP<minNP) = NaN;
staticSeriesStd(staticSeriesNP<minNP) = NaN;
staticSeriesDistMid(staticSeriesNP<minNP) = NaN;
staticSeriesDistHR(staticSeriesNP<minNP) = NaN;

staticSeriesBulkDistMid(staticSeriesNP(:,4)<minNP) = NaN;
staticSeriesBulkDistHR(staticSeriesNP(:,4)<minNP) = NaN;

dynamicSeriesBefMean(dynamicSeriesBefNP<minNP) = NaN;
dynamicSeriesBefStd(dynamicSeriesBefNP<minNP) = NaN;
dynamicSeriesBefDistMid(dynamicSeriesBefNP<minNP) = NaN;
dynamicSeriesBefDistHR(dynamicSeriesBefNP<minNP) = NaN;

dynamicSeriesAftMean(dynamicSeriesAftNP<minNP) = NaN;
dynamicSeriesAftStd(dynamicSeriesAftNP<minNP) = NaN;
dynamicSeriesAftDistMid(dynamicSeriesAftNP<minNP) = NaN;
dynamicSeriesAftDistHR(dynamicSeriesAftNP<minNP) = NaN;

edgeSeriesMean(edgeSeriesNP<minNP) = NaN;
edgeSeriesStd(edgeSeriesNP<minNP) = NaN;
edgeSeriesDistMid(edgeSeriesNP<minNP) = NaN;
edgeSeriesDistHR(edgeSeriesNP<minNP) = NaN;

combDynamicSeriesAftMean(combDynamicSeriesAftNP<minNP) = NaN;
combDynamicSeriesAftStd(combDynamicSeriesAftNP<minNP) = NaN;
combDynamicSeriesAftDistMid(combDynamicSeriesAftNP<minNP) = NaN;
combDynamicSeriesAftDistHR(combDynamicSeriesAftNP<minNP) = NaN;

%get edge position over time
edgePos = staticSeriesDistMid(:,1) - staticSeriesDistHR(:,1);

%compensate for window area over-estimation if requested
if compArea
    edgeDispBef = abs(nanmean(diff(edgePos(1:-minInc+1))));
    edgeDispAft = nanmean(diff(edgePos(-minInc+1:end)));
    winLengthBef = 2 * nanmean(edgeSeriesDistHR(1:-minInc+1));
    winLengthAft = 2 * nanmean(edgeSeriesDistHR(-minInc+2:end));
    factorBef = 1 - edgeDispBef/winLengthBef/8;
    factorAft = 1 - edgeDispAft/winLengthAft/8;
    edgeSeriesMean = edgeSeriesMean ./ [factorBef*ones(-minInc,1); 1; factorAft*ones(maxInc,1)];
    edgeSeriesStd = edgeSeriesStd ./ [factorBef*ones(-minInc,1); 1; factorAft*ones(maxInc,1)];
end

%define color vectors for plotting
%green to magenta in 11 steps
incColor = repmat((1:-0.1:0)',1,3).*repmat([0 1 0],11,1) + repmat((0:0.1:1)',1,3).*repmat([1 0 1],11,1);
incColor = min(incColor,1);
%black, blue, green, magenta
bandColor = [0 0 0; 0 0 1; 0 1 0; 1 0 1];

%define maximum number of increments and bands to plot
maxIncPlot = min(maxInc,size(incColor,1));
maxBandPlot = min(numBand,size(bandColor,1));

%calcualte SEMs
dynamicSeriesBefSem = dynamicSeriesBefStd./sqrt(dynamicSeriesBefNP);
staticSeriesSem = staticSeriesStd./sqrt(staticSeriesNP);
dynamicSeriesAftSem = dynamicSeriesAftStd./sqrt(dynamicSeriesAftNP);
combDynamicSeriesAftSem = combDynamicSeriesAftStd./sqrt(combDynamicSeriesAftNP);
edgeSeriesSem = edgeSeriesStd./sqrt(edgeSeriesNP);


%bands
for iBand = 1 : 3
    inputDataTmp.observations = [staticSeriesMean(:,iBand) staticSeriesSem(:,iBand)];
    inputDataTmp.time = [(minInc:maxInc)'*convFactTime zeros(maxInc-minInc+1,1)];
    inputData = convertTrajectoryData(inputDataTmp);
    trajectoryAnalysis(inputData,struct('saveTxt',0));
end
inputDataTmp.observations = [mean(staticSeriesMean(:,4:end),2) mean(staticSeriesSem(:,4:end),2)];
inputDataTmp.time = [(minInc:maxInc)'*convFactTime zeros(maxInc-minInc+1,1)];
inputData = convertTrajectoryData(inputDataTmp);
trajectoryAnalysis(inputData,struct('saveTxt',0));


%right behind edge
inputDataTmp.observations = [edgeSeriesMean(:,1) edgeSeriesSem(:,1)];
inputDataTmp.time = [(minInc:maxInc)'*convFactTime zeros(maxInc-minInc+1,1)];
inputData = convertTrajectoryData(inputDataTmp);
trajectoryAnalysis(inputData,struct('saveTxt',0));

%new cell area aligned & combined
inputDataTmp.observations = [combDynamicSeriesAftMean(:,1) combDynamicSeriesAftSem(:,1)];
inputDataTmp.time = [(minInc:maxInc)'*convFactTime zeros(maxInc-minInc+1,1)];
inputData = convertTrajectoryData(inputDataTmp);
trajectoryAnalysis(inputData,struct('saveTxt',0));



