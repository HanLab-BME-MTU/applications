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
[numRow,numBand,numMode] = size(staticSeriesMean);

%define which modes to plot
if nargin < 2 || isempty(mode2plot)
    mode2plot = 1 : numMode;
end
mode2plot = mode2plot(:)';
numMode2plot = length(mode2plot);
    
%dynamic series before
dynamicSeriesBefMean = [particleBehavior.befDynamic.mean(end:-1:1,:,:); ...
    particleBehavior.onset.mean; NaN(maxInc,numBand,numMode)];
dynamicSeriesBefStd = [particleBehavior.befDynamic.std(end:-1:1,:,:); ...
    particleBehavior.onset.std; NaN(maxInc,numBand,numMode)];
dynamicSeriesBefNP = [particleBehavior.befDynamic.numPoints(end:-1:1,:,:); ...
    particleBehavior.onset.numPoints; ones(maxInc,numBand,numMode)];

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

%define color vectors for plotting
incColor = [1 0 0; 0 1 0; 1 0 1; 0 1 1; 1 1 0; 0.7 0.7 0.7];
bandColor = [0 0 0; 0 0 1; incColor];

%open figure
hFig = figure('Name',figureName);
hold on

%plot for each mode
for jMode = 1 : numMode2plot
    
    %extract mode information
    iMode = mode2plot(jMode);
    dynamicSeriesBefMeanMode = dynamicSeriesBefMean(:,:,iMode);
    dynamicSeriesBefSemMode = dynamicSeriesBefStd(:,:,iMode)./sqrt(dynamicSeriesBefNP(:,:,iMode));
    staticSeriesMeanMode = staticSeriesMean(:,:,iMode);
    staticSeriesSemMode = staticSeriesStd(:,:,iMode)./sqrt(staticSeriesNP(:,:,iMode));
    dynamicSeriesAftMeanMode = dynamicSeriesAftMean(:,:,iMode);
    dynamicSeriesAftSemMode = dynamicSeriesAftStd(:,:,iMode)./sqrt(dynamicSeriesAftNP(:,:,iMode));
    combDynamicSeriesAftMeanMode = combDynamicSeriesAftMean(:,:,iMode);
    combDynamicSeriesAftSemMode = combDynamicSeriesAftStd(:,:,iMode)./sqrt(combDynamicSeriesAftNP(:,:,iMode));
    
    %determine y-axis limits
    yMax = max([dynamicSeriesBefMeanMode(:)+dynamicSeriesBefSemMode(:); ...
        staticSeriesMeanMode(:)+staticSeriesSemMode(:); ...
        dynamicSeriesAftMeanMode(:)+dynamicSeriesAftSemMode(:); ...
        combDynamicSeriesAftMeanMode(:)+combDynamicSeriesAftSemMode(:)]);
    yMin = min([dynamicSeriesBefMeanMode(:)-dynamicSeriesBefSemMode(:); ...
        staticSeriesMeanMode(:)-staticSeriesSemMode(:); ...
        dynamicSeriesAftMeanMode(:)-dynamicSeriesAftSemMode(:); ...
        combDynamicSeriesAftMeanMode(:)-combDynamicSeriesAftSemMode(:)]);
    if isnan(yMax)
        yMax = max([dynamicSeriesBefMeanMode(:); ...
            staticSeriesMeanMode(:); ...
            dynamicSeriesAftMeanMode(:); ...
            combDynamicSeriesAftMeanMode(:)]);
        yMin = min([dynamicSeriesBefMeanMode(:); ...
            staticSeriesMeanMode(:); ...
            dynamicSeriesAftMeanMode(:); ...
            combDynamicSeriesAftMeanMode(:)]);
    end
    if yMax==yMin
        yMax = yMin*1.01 + eps;
    end
    
    if ~isnan(yMax)
        
        %first plot
        subplot(numMode2plot,2,(jMode-1)*2+1), hold on
        
        %before dynamic
        plot(minInc:maxInc,dynamicSeriesBefMeanMode(:,1),'k--','Marker','.')
        myErrorbar(minInc:maxInc,dynamicSeriesBefMeanMode(:,1),dynamicSeriesBefSemMode(:,1))
        
        %static series
        plot(minInc:maxInc,staticSeriesMeanMode(:,1),'k','Marker','.','LineWidth',2)
        myErrorbar(minInc:maxInc,staticSeriesMeanMode(:,1),staticSeriesSemMode(:,1))
        
        %after dynamic
        legendEntries = cell(1,6);
        for iInc = 1 : 6
            plot(minInc:maxInc,dynamicSeriesAftMeanMode(:,iInc),'color',incColor(iInc,:),'Marker','.')
            myErrorbar(minInc:maxInc,dynamicSeriesAftMeanMode(:,iInc),dynamicSeriesAftSemMode(:,iInc))
            legendEntries{iInc} = ['Dynamic after Inc ' num2str(iInc)];
        end
        
        %after dynamic combined
        plot(minInc:maxInc,combDynamicSeriesAftMeanMode(:,1),'b','Marker','.','LineWidth',2)
        myErrorbar(minInc:maxInc,combDynamicSeriesAftMeanMode(:,1),combDynamicSeriesAftSemMode(:,1))
        
        axis([minInc maxInc yMin yMax]);
        xlabel('Time from protruion onset (frames)')
        
        legendEntries = [{'Dynamic before'} {'Static'} legendEntries {'Dynamic after aligned & combined'}]; %#ok<AGROW>
        legend(legendEntries);
        
        %second plot
        subplot(numMode2plot,2,jMode*2), hold on
        
        %bands
        legendEntries = cell(1,8);
        for iBand = 1 : 2
            plot(minInc:maxInc,staticSeriesMeanMode(:,iBand),'Color',bandColor(iBand,:),'Marker','.','LineWidth',2)
            myErrorbar(minInc:maxInc,staticSeriesMeanMode(:,iBand),staticSeriesSemMode(:,iBand))
            legendEntries{iBand} = ['Band ' num2str(iBand)];
        end
        for iBand = 3 : min(numBand,8)
            plot(minInc:maxInc,staticSeriesMeanMode(:,iBand),'Color',bandColor(iBand,:),'Marker','.')
            myErrorbar(minInc:maxInc,staticSeriesMeanMode(:,iBand),staticSeriesSemMode(:,iBand))
            legendEntries{iBand} = ['Band ' num2str(iBand)];
        end
        axis([minInc maxInc yMin yMax]);
        xlabel('Time from protruion onset (frames)')
        legend(legendEntries);
        for iBand = 1 : min(numBand,8)
            plot(minInc:maxInc,dynamicSeriesBefMeanMode(:,iBand),'Color',bandColor(iBand,:),'Marker','.','LineStyle','--')
            myErrorbar(minInc:maxInc,dynamicSeriesBefMeanMode(:,iBand),dynamicSeriesBefSemMode(:,iBand))
        end
        
    end %(if ~isnan(yMax))
    
end %(for jMode = 1 : numMode2plot)

if ~isempty(saveLoc)
    saveas(hFig,saveLoc)
end
