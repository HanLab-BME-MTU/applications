function plotSptRelToActivityOnsetAdaptiveWindows(particleBehavior,...
    windowDistFromEdge,mode2plot,minNP,figureName,saveLoc,convFact,...
    yAxisLabel,axisLimits)

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

dynamicSeriesBefMean(dynamicSeriesBefNP<minNP) = NaN;
dynamicSeriesBefStd(dynamicSeriesBefNP<minNP) = NaN;
dynamicSeriesBefDistMid(dynamicSeriesBefNP<minNP) = NaN;
dynamicSeriesBefDistHR(dynamicSeriesBefNP<minNP) = NaN;

dynamicSeriesAftMean(dynamicSeriesAftNP<minNP) = NaN;
dynamicSeriesAftStd(dynamicSeriesAftNP<minNP) = NaN;
dynamicSeriesAftDistMid(dynamicSeriesAftNP<minNP) = NaN;
dynamicSeriesAftDistHR(dynamicSeriesAftNP<minNP) = NaN;

combDynamicSeriesAftMean(combDynamicSeriesAftNP<minNP) = NaN;
combDynamicSeriesAftStd(combDynamicSeriesAftNP<minNP) = NaN;
combDynamicSeriesAftDistMid(combDynamicSeriesAftNP<minNP) = NaN;
combDynamicSeriesAftDistHR(combDynamicSeriesAftNP<minNP) = NaN;

%define color vectors for plotting
%red, green, brown, gray, magenta, cyan, yellow, purple
incColor = [1 0 0; 0 1 0; 0.8 0.5 0.2; 0.7 0.7 0.7; 1 0 1; 0 1 1; 1 1 0; 0.6 0.5 1];
bandColor = [0 0 0; 0 0 1; incColor];

%define maximum number of increments and bands to plot
maxIncPlot = min(maxInc,size(incColor,1));
maxBandPlot = min(numBand,size(bandColor,1));

%open figure
hFig = figure('Name',figureName);
hold on

%calcualte SEMs
dynamicSeriesBefSem = dynamicSeriesBefStd./sqrt(dynamicSeriesBefNP);
staticSeriesSem = staticSeriesStd./sqrt(staticSeriesNP);
dynamicSeriesAftSem = dynamicSeriesAftStd./sqrt(dynamicSeriesAftNP);
combDynamicSeriesAftSem = combDynamicSeriesAftStd./sqrt(combDynamicSeriesAftNP);

%determine y-axis limits
if isempty(axisLimits)
    yMax = max([dynamicSeriesBefMean(:)+dynamicSeriesBefSem(:); ...
        staticSeriesMean(:)+staticSeriesSem(:); ...
        dynamicSeriesAftMean(:)+dynamicSeriesAftSem(:); ...
        combDynamicSeriesAftMean(:)+combDynamicSeriesAftSem(:)]);
    yMin = min([dynamicSeriesBefMean(:)-dynamicSeriesBefSem(:); ...
        staticSeriesMean(:)-staticSeriesSem(:); ...
        dynamicSeriesAftMean(:)-dynamicSeriesAftSem(:); ...
        combDynamicSeriesAftMean(:)-combDynamicSeriesAftSem(:)]);
    if isnan(yMax)
        yMax = max([dynamicSeriesBefMean(:); ...
            staticSeriesMean(:); ...
            dynamicSeriesAftMean(:); ...
            combDynamicSeriesAftMean(:)]);
        yMin = min([dynamicSeriesBefMean(:); ...
            staticSeriesMean(:); ...
            dynamicSeriesAftMean(:); ...
            combDynamicSeriesAftMean(:)]);
    end
    if yMax==yMin
        yMax = yMin*1.01 + eps;
    end
    xMin = minInc*convFactTime;
    xMax = maxInc*convFactTime;
else
    xMin = axisLimits(1);
    xMax = axisLimits(2);
    yMin = axisLimits(3);
    yMax = axisLimits(4);
end

if ~isnan(yMax)
    
    %plot 1 - values at the edge
    subplot(4,2,[1 3]), hold on
    
    %before dynamic
    plot((minInc:maxInc)*convFactTime,dynamicSeriesBefMean(:,1),'k--','Marker','.')
    myErrorbar((minInc:maxInc)*convFactTime,dynamicSeriesBefMean(:,1),dynamicSeriesBefSem(:,1))
    
    %static series
    plot((minInc:maxInc)*convFactTime,staticSeriesMean(:,1),'k','Marker','.','LineWidth',2)
    myErrorbar((minInc:maxInc)*convFactTime,staticSeriesMean(:,1),staticSeriesSem(:,1))
    
    %after dynamic
    legendEntries = cell(1,maxIncPlot);
    for iInc = 1 : maxIncPlot
        plot((minInc:maxInc)*convFactTime,dynamicSeriesAftMean(:,iInc),'color',incColor(iInc,:),'Marker','.')
        myErrorbar((minInc:maxInc)*convFactTime,dynamicSeriesAftMean(:,iInc),dynamicSeriesAftSem(:,iInc))
        legendEntries{iInc} = ['Dynamic after Inc ' num2str(iInc)];
    end
    
    %after dynamic combined
    plot((minInc:maxInc)*convFactTime,combDynamicSeriesAftMean(:,1),'b','Marker','.','LineWidth',2)
    myErrorbar((minInc:maxInc)*convFactTime,combDynamicSeriesAftMean(:,1),combDynamicSeriesAftSem(:,1))
    
    axis([xMin xMax yMin yMax]);
    xlabel('Time from protrusion onset (s)')
    ylabel(yAxisLabel)
    
    legendEntries = [{'Dynamic before'} {'Static'} legendEntries {'Dynamic after aligned & combined'}];
    legend(legendEntries);
    
    %line indicating activity onset
    plot([0 0],[yMin yMax],'k:')
    
    %plot 2 - values in bands
    subplot(4,2,[2 4]), hold on
    
    legendEntries = cell(1,maxBandPlot);
    for iBand = 1 : 2
        plot((minInc:maxInc)*convFactTime,staticSeriesMean(:,iBand),'Color',bandColor(iBand,:),'Marker','.','LineWidth',2)
        myErrorbar((minInc:maxInc)*convFactTime,staticSeriesMean(:,iBand),staticSeriesSem(:,iBand))
        legendEntries{iBand} = ['Band ' num2str(iBand)];
    end
    for iBand = 3 : maxBandPlot
        plot((minInc:maxInc)*convFactTime,staticSeriesMean(:,iBand),'Color',bandColor(iBand,:),'Marker','.')
        myErrorbar((minInc:maxInc)*convFactTime,staticSeriesMean(:,iBand),staticSeriesSem(:,iBand))
        legendEntries{iBand} = ['Band ' num2str(iBand)];
    end
    axis([xMin xMax yMin yMax]);
    xlabel('Time from protrusion onset (s)')
    ylabel(yAxisLabel)
    legend(legendEntries);
    for iBand = 1 : maxBandPlot
        plot((minInc:maxInc)*convFactTime,dynamicSeriesBefMean(:,iBand),'Color',bandColor(iBand,:),'Marker','.','LineStyle','--')
        myErrorbar((minInc:maxInc)*convFactTime,dynamicSeriesBefMean(:,iBand),dynamicSeriesBefSem(:,iBand))
    end
    
    %line indicating activity onset
    plot([0 0],[yMin yMax],'k:')
    
    %plot 3a - distances at the edge, view 1
    subplot(4,2,5,'YDir','reverse'), hold on
    
    tmp = [dynamicSeriesBefDistMid(:,1)+dynamicSeriesBefDistHR(:,1) ...
        staticSeriesDistMid(:,1)+staticSeriesDistHR(:,1) ...
        dynamicSeriesAftDistMid+dynamicSeriesAftDistHR ...
        combDynamicSeriesAftDistMid(:,1)+combDynamicSeriesAftDistHR(:,1)];
    yMin = -0.1;
    yMax2 = max(tmp(:));
    yMax2(isnan(yMax2)) = eps;
    
    %line indicating activity onset
    plot([0 0],[yMin yMax2],'k:')
    
    %static series
    plot((minInc:maxInc)*convFactTime,staticSeriesDistMid(:,1),'k','LineStyle','none','LineWidth',2)
    myErrorbar((minInc:maxInc)*convFactTime,staticSeriesDistMid(:,1),staticSeriesDistHR(:,1))
    
    %after dynamic
    for iInc = 1 : maxIncPlot
        plot((minInc:maxInc)*convFactTime,dynamicSeriesAftDistMid(:,iInc),'color',incColor(iInc,:),'LineStyle','none')
        myErrorbar((minInc:maxInc)*convFactTime,dynamicSeriesAftDistMid(:,iInc),dynamicSeriesAftDistHR(:,iInc))
    end
    
    axis([xMin xMax yMin yMax2]);
    xlabel('Time from protrusion onset (s)')
    ylabel('Window distance from cell edge (um)')
    
    %plot 3b - distances at the edge, view 2
    subplot(4,2,7,'YDir','reverse'), hold on
    
    %line indicating activity onset
    plot([0 0],[yMin yMax2],'k:')
    
    %static series
    plot((minInc:maxInc)*convFactTime,staticSeriesDistMid(:,1),'k','LineStyle','none','LineWidth',2)
    myErrorbar((minInc:maxInc)*convFactTime,staticSeriesDistMid(:,1),staticSeriesDistHR(:,1))
    
    %after dynamic combined
    plot((minInc:maxInc)*convFactTime,combDynamicSeriesAftDistMid(:,1),'b','LineStyle','none','LineWidth',2)
    myErrorbar((minInc:maxInc)*convFactTime,combDynamicSeriesAftDistMid(:,1),combDynamicSeriesAftDistHR(:,1))
    
    axis([xMin xMax yMin yMax2]);
    xlabel({'Time from protrusion onset (black)/','new cell-matrix contact (blue) (s)'})
    ylabel('Window distance from cell edge (um)')
    
    %plot 4 - distances in bands
    subplot(4,2,[6 8],'YDir','reverse'), hold on
    
    tmp = staticSeriesDistMid+staticSeriesDistHR;
    yMax2 = max(tmp(:));
    yMax2(isnan(yMax2)) = eps;
    
    %line indicating activity onset
    plot([0 0],[yMin yMax2],'k:')
    
    for iBand = 1 : maxBandPlot
        plot((minInc:maxInc)*convFactTime,staticSeriesDistMid(:,iBand),'Color',bandColor(iBand,:),'LineStyle','none')
        myErrorbar((minInc:maxInc)*convFactTime,staticSeriesDistMid(:,iBand),staticSeriesDistHR(:,iBand))
    end
    
    axis([xMin xMax yMin yMax2]);
    xlabel('Time from protrusion onset (s)')
    ylabel('Window distance from cell edge (um)')
    
end %(if ~isnan(yMax))

if ~isempty(saveLoc)
    saveas(hFig,saveLoc)
end

