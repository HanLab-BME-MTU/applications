function plotSptRelToActivityOnsetAdaptiveWindows(particleBehavior,minNP,figureName,saveLoc)

if nargin < 2 || iesmpty(minNP)
    minNP = 20;
end

if nargin < 3 || isempty(figureName)
    figureName = 'test';
end

if nargin < 4 || isempty(saveLoc)
    saveLoc = [];
end

%get increment range
minInc = -size(particleBehavior.befStatic.mean,1);
maxInc = size(particleBehavior.aftStatic.mean,1);

%construct series of behavior to plot

%static series
staticSeriesMean = [particleBehavior.befStatic.mean(end:-1:1,:,:); ...
    particleBehavior.onset.mean; particleBehavior.aftStatic.mean];
staticSeriesStd = [particleBehavior.befStatic.std(end:-1:1,:,:); ...
    particleBehavior.onset.std; particleBehavior.aftStatic.std];
staticSeriesNP = [particleBehavior.befStatic.numPoints(end:-1:1,:,:); ...
    particleBehavior.onset.numPoints; particleBehavior.aftStatic.numPoints];

%dynamic series before
dynamicSeriesBefMean = [particleBehavior.befDynamic.mean(end:-1:1,:,:); ...
    particleBehavior.onset.mean];
dynamicSeriesBefStd = [particleBehavior.befDynamic.std(end:-1:1,:,:); ...
    particleBehavior.onset.std];
dynamicSeriesBefNP = [particleBehavior.befDynamic.numPoints(end:-1:1,:,:); ...
    particleBehavior.onset.numPoints];

%dynamic series after
dynamicSeriesAftMean = particleBehavior.aftDynamic.mean;
dynamicSeriesAftStd = particleBehavior.aftDynamic.std;
dynamicSeriesAftNP = particleBehavior.aftDynamic.numPoints;

%remove measurement with less than minNP points
staticSeriesMean(staticSeriesNP<minNP) = NaN;
staticSeriesStd(staticSeriesNP<minNP) = NaN;
dynamicSeriesBefMean(dynamicSeriesBefNP<minNP) = NaN;
dynamicSeriesBefStd(dynamicSeriesBefNP<minNP) = NaN;
dynamicSeriesAftMean(dynamicSeriesAftNP<minNP) = NaN;
dynamicSeriesAftStd(dynamicSeriesAftNP<minNP) = NaN;

%get number of modes
numMode = size(staticSeriesMean,3);

%define color vector for plotting
incColor = [1 0 0; 0 1 0; 1 0 1; 0 1 1; 1 1 0; 0.5 0.5 0.5];

%open figure
hFig = figure('Name',figureName);
hold on

%plot for each mode
for iMode = 1 : numMode
    
    subplot(numMode,1,iMode), hold on
    
    plot(minInc:maxInc,staticSeriesMean(:,1,iMode),'k','Marker','.')
    myErrorbar(minInc:maxInc,staticSeriesMean(:,1,iMode),staticSeriesStd(:,1,iMode)./sqrt(staticSeriesNP(:,1,iMode)))
    
    plot(minInc:0,dynamicSeriesBefMean(:,1,iMode),'b','Marker','.')
    myErrorbar(minInc:0,dynamicSeriesBefMean(:,1,iMode),dynamicSeriesBefStd(:,1,iMode)./sqrt(dynamicSeriesBefNP(:,1,iMode)))
    
    legendEntries = cell(1,6);
    for iInc = 1 : 6
        plot(1:maxInc,dynamicSeriesAftMean(:,iInc,iMode),'color',incColor(iInc,:),'Marker','.')
        myErrorbar(1:maxInc,dynamicSeriesAftMean(:,iInc,iMode),dynamicSeriesAftStd(:,iInc,iMode)./sqrt(dynamicSeriesAftNP(:,iInc,iMode)))
        legendEntries{iInc} = ['Dynamic after Inc ' num2str(iInc)];
    end
    
    legendEntries = [{'Static'} {'Dynamic before'} legendEntries]; %#ok<AGROW>
    legend(legendEntries);
    
end

if ~isempty(saveLoc)
    saveas(hFig,saveLoc)
end


