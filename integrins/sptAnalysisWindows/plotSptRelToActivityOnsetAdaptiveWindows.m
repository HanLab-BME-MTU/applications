function plotSptRelToActivityOnsetAdaptiveWindows(particleBehavior,...
    mode2plot,minNP,figureName,saveLoc)

if nargin < 3 || isempty(minNP)
    minNP = 20;
end

if nargin < 4 || isempty(figureName)
    figureName = 'test';
end

if nargin < 5 || isempty(saveLoc)
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

%get number of modes
numMode = size(staticSeriesMean,3);

%define which modes to plot
if nargin < 2 || isempty(mode2plot)
    mode2plot = 1 : numMode;
end
mode2plot = mode2plot(:)';
numMode2plot = length(mode2plot);
    
%dynamic series before
dynamicSeriesBefMean = [particleBehavior.befDynamic.mean(end:-1:1,:,:); ...
    particleBehavior.onset.mean; NaN(maxInc,1,numMode)];
dynamicSeriesBefStd = [particleBehavior.befDynamic.std(end:-1:1,:,:); ...
    particleBehavior.onset.std; NaN(maxInc,1,numMode)];
dynamicSeriesBefNP = [particleBehavior.befDynamic.numPoints(end:-1:1,:,:); ...
    particleBehavior.onset.numPoints; ones(maxInc,1,numMode)];

%dynamic series after
dynamicSeriesAftMean = [NaN(-minInc+1,maxInc,numMode); particleBehavior.aftDynamic.mean];
dynamicSeriesAftStd = [NaN(-minInc+1,maxInc,numMode); particleBehavior.aftDynamic.std];
dynamicSeriesAftNP = [ones(-minInc+1,maxInc,numMode); particleBehavior.aftDynamic.numPoints];

%combined dynamic series after
combDynamicSeriesAftMean = [NaN(-minInc,1,numMode); particleBehavior.aftDynamicComb.mean; NaN(1,1,numMode)];
combDynamicSeriesAftStd = [NaN(-minInc,1,numMode); particleBehavior.aftDynamicComb.std; NaN(1,1,numMode)];
combDynamicSeriesAftNP = [ones(-minInc,1,numMode); particleBehavior.aftDynamicComb.numPoints; ones(1,1,numMode)];

%remove measurements with less than minNP points
staticSeriesMean(staticSeriesNP<minNP) = NaN;
staticSeriesStd(staticSeriesNP<minNP) = NaN;
dynamicSeriesBefMean(dynamicSeriesBefNP<minNP) = NaN;
dynamicSeriesBefStd(dynamicSeriesBefNP<minNP) = NaN;
dynamicSeriesAftMean(dynamicSeriesAftNP<minNP) = NaN;
dynamicSeriesAftStd(dynamicSeriesAftNP<minNP) = NaN;
combDynamicSeriesAftMean(combDynamicSeriesAftNP<minNP) = NaN;
combDynamicSeriesAftStd(combDynamicSeriesAftNP<minNP) = NaN;

%define color vector for plotting
incColor = [1 0 0; 0 1 0; 1 0 1; 0 1 1; 1 1 0; 0.7 0.7 0.7];

%open figure
hFig = figure('Name',figureName);
hold on

%plot for each mode
for jMode = 1 : numMode2plot
    
    %specify subplot
    subplot(numMode2plot,1,jMode), hold on
    
    %get mode to plot
    iMode = mode2plot(jMode);
    
    %before dynamic
    plot(minInc:maxInc,dynamicSeriesBefMean(:,1,iMode),'b','Marker','.')
    myErrorbar(minInc:maxInc,dynamicSeriesBefMean(:,1,iMode),dynamicSeriesBefStd(:,1,iMode)./sqrt(dynamicSeriesBefNP(:,1,iMode)))
    
    %static series
    plot(minInc:maxInc,staticSeriesMean(:,1,iMode),'k','Marker','.','LineWidth',2)
    myErrorbar(minInc:maxInc,staticSeriesMean(:,1,iMode),staticSeriesStd(:,1,iMode)./sqrt(staticSeriesNP(:,1,iMode)))
    
    %after dynamic
    legendEntries = cell(1,6);
    for iInc = 1 : 6
        plot(minInc:maxInc,dynamicSeriesAftMean(:,iInc,iMode),'color',incColor(iInc,:),'Marker','.')
        myErrorbar(minInc:maxInc,dynamicSeriesAftMean(:,iInc,iMode),dynamicSeriesAftStd(:,iInc,iMode)./sqrt(dynamicSeriesAftNP(:,iInc,iMode)))
        legendEntries{iInc} = ['Dynamic after Inc ' num2str(iInc)];
    end
    
    %after dynamic combined
    plot(minInc:maxInc,combDynamicSeriesAftMean(:,1,iMode),'b','Marker','.','LineWidth',2)
    myErrorbar(minInc:maxInc,combDynamicSeriesAftMean(:,1,iMode),combDynamicSeriesAftStd(:,1,iMode)./sqrt(combDynamicSeriesAftNP(:,1,iMode)))    
    
    legendEntries = [{'Dynamic before'} {'Static'} legendEntries {'Dynamic after aligned & combined'}]; %#ok<AGROW>
    legend(legendEntries);
    
end

if ~isempty(saveLoc)
    saveas(hFig,saveLoc)
end


