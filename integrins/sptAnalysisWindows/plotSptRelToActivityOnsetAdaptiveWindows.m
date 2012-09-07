function plotSptRelToActivityOnsetAdaptiveWindows(particleBehavior,...
    windowDistFromEdge,mode2plot,minNP,figureName,saveLoc,convFact,yAxisLabel)

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

%get increment range
minInc = -size(particleBehavior.befStatic.mean,1);
maxInc = size(particleBehavior.aftStatic.mean,1);

%construct series of behavior to plot

%static series
staticSeriesMean = [particleBehavior.befStatic.mean(end:-1:1,:,:); ...
    particleBehavior.onset.mean; particleBehavior.aftStatic.mean] * convFactProp;
staticSeriesStd = [particleBehavior.befStatic.std(end:-1:1,:,:); ...
    particleBehavior.onset.std; particleBehavior.aftStatic.std] * convFactProp;
staticSeriesNP = [particleBehavior.befStatic.numPoints(end:-1:1,:,:); ...
    particleBehavior.onset.numPoints; particleBehavior.aftStatic.numPoints];
staticSeriesDist = [windowDistFromEdge.befStatic(end:-1:1,:,:,1); ...
    windowDistFromEdge.onset(:,:,:,1); windowDistFromEdge.aftStatic(:,:,:,1)] ...
    * convFactDist;
tmp = staticSeriesDist;
staticSeriesDist(:,:,1) = (tmp(:,:,1) + tmp(:,:,2)) / 2;
staticSeriesDist(:,:,2) = (tmp(:,:,2) - tmp(:,:,1)) / 2;

%get number of modes
[numRow,numBand,numMode] = size(staticSeriesMean);

%define which modes to plot
if nargin < 3 || isempty(mode2plot)
    mode2plot = 1 : numMode;
end
mode2plot = mode2plot(:)';
numMode2plot = length(mode2plot);

%dynamic series before
dynamicSeriesBefMean = [particleBehavior.befDynamic.mean(end:-1:1,:,:); ...
    particleBehavior.onset.mean; NaN(maxInc,numBand,numMode)] * convFactProp;
dynamicSeriesBefStd = [particleBehavior.befDynamic.std(end:-1:1,:,:); ...
    particleBehavior.onset.std; NaN(maxInc,numBand,numMode)] * convFactProp;
dynamicSeriesBefNP = [particleBehavior.befDynamic.numPoints(end:-1:1,:,:); ...
    particleBehavior.onset.numPoints; ones(maxInc,numBand,numMode)];
dynamicSeriesBefDist = [windowDistFromEdge.befDynamic(end:-1:1,:,:,1); ...
    windowDistFromEdge.onset(:,:,:,1); NaN(maxInc,numBand,2)] ...
    * convFactDist;
tmp = dynamicSeriesBefDist;
dynamicSeriesBefDist(:,:,1) = (tmp(:,:,1) + tmp(:,:,2)) / 2;
dynamicSeriesBefDist(:,:,2) = (tmp(:,:,2) - tmp(:,:,1)) / 2;

%dynamic series after
dynamicSeriesAftMean = [NaN(-minInc+1,maxInc,numMode); particleBehavior.aftDynamic.mean] ...
    * convFactProp;
dynamicSeriesAftStd = [NaN(-minInc+1,maxInc,numMode); particleBehavior.aftDynamic.std] ...
    * convFactProp;
dynamicSeriesAftNP = [ones(-minInc+1,maxInc,numMode); particleBehavior.aftDynamic.numPoints];
dynamicSeriesAftDist = [NaN(-minInc+1,maxInc,2); windowDistFromEdge.aftDynamic(:,:,:,1)] ...
    * convFactDist;
tmp = dynamicSeriesAftDist;
dynamicSeriesAftDist(:,:,1) = (tmp(:,:,1) + tmp(:,:,2)) / 2;
dynamicSeriesAftDist(:,:,2) = (tmp(:,:,2) - tmp(:,:,1)) / 2;

%combined dynamic series after
combDynamicSeriesAftMean = [NaN(-minInc,1,numMode); particleBehavior.aftDynamicComb.mean; ...
    NaN(1,1,numMode)] * convFactProp;
combDynamicSeriesAftStd = [NaN(-minInc,1,numMode); particleBehavior.aftDynamicComb.std; ...
    NaN(1,1,numMode)] * convFactProp;
combDynamicSeriesAftNP = [ones(-minInc,1,numMode); particleBehavior.aftDynamicComb.numPoints; ones(1,1,numMode)];
combDynamicSeriesAftDist = [NaN(-minInc,1,2); windowDistFromEdge.aftDynamicComb(:,:,:,1); ...
    NaN(1,1,2)] * convFactDist;
tmp = combDynamicSeriesAftDist;
combDynamicSeriesAftDist(:,:,1) = (tmp(:,:,1) + tmp(:,:,2)) / 2;
combDynamicSeriesAftDist(:,:,2) = (tmp(:,:,2) - tmp(:,:,1)) / 2;

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
incColor = [1 0 0; 0 1 0; 0.8 0.5 0.2; 0.7 0.7 0.7; 1 0 1; 0 1 1];
% incColor = [(1:-0.2:0)' (0:0.2:1)' zeros(6,1)]; %[(0:0.4:1)'; (1:-0.4:0)'] ];
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
        
        %plot 1 - values at the edge
        subplot(numMode2plot+1,2,(jMode-1)*2+1), hold on
        
        %line indicating activity onset
        plot([0 0],[yMin yMax],'k--')
        
        %before dynamic
        plot((minInc:maxInc)*convFactTime,dynamicSeriesBefMeanMode(:,1),'k--','Marker','.')
        myErrorbar((minInc:maxInc)*convFactTime,dynamicSeriesBefMeanMode(:,1),dynamicSeriesBefSemMode(:,1))
        
        %static series
        plot((minInc:maxInc)*convFactTime,staticSeriesMeanMode(:,1),'k','Marker','.','LineWidth',2)
        myErrorbar((minInc:maxInc)*convFactTime,staticSeriesMeanMode(:,1),staticSeriesSemMode(:,1))
        
        %after dynamic
        legendEntries = cell(1,6);
        for iInc = 1 : 6
            plot((minInc:maxInc)*convFactTime,dynamicSeriesAftMeanMode(:,iInc),'color',incColor(iInc,:),'Marker','.')
            myErrorbar((minInc:maxInc)*convFactTime,dynamicSeriesAftMeanMode(:,iInc),dynamicSeriesAftSemMode(:,iInc))
            legendEntries{iInc} = ['Dynamic after Inc ' num2str(iInc)];
        end
        
        %after dynamic combined
        plot((minInc:maxInc)*convFactTime,combDynamicSeriesAftMeanMode(:,1),'b','Marker','.','LineWidth',2)
        myErrorbar((minInc:maxInc)*convFactTime,combDynamicSeriesAftMeanMode(:,1),combDynamicSeriesAftSemMode(:,1))
        
        axis([minInc*convFactTime maxInc*convFactTime yMin yMax]);
        xlabel('Time from protruion onset (s)')
        ylabel(yAxisLabel)
        
        legendEntries = [{'Dynamic before'} {'Static'} legendEntries {'Dynamic after aligned & combined'}]; %#ok<AGROW>
        legend(legendEntries);
        
        
        %plot 2 - values in bands
        subplot(numMode2plot+1,2,jMode*2), hold on
        
        %line indicating activity onset
        plot([0 0],[yMin yMax],'k--')
        
        legendEntries = cell(1,8);
        for iBand = 1 : 2
            plot((minInc:maxInc)*convFactTime,staticSeriesMeanMode(:,iBand),'Color',bandColor(iBand,:),'Marker','.','LineWidth',2)
            myErrorbar((minInc:maxInc)*convFactTime,staticSeriesMeanMode(:,iBand),staticSeriesSemMode(:,iBand))
            legendEntries{iBand} = ['Band ' num2str(iBand)];
        end
        for iBand = 3 : min(numBand,8)
            plot((minInc:maxInc)*convFactTime,staticSeriesMeanMode(:,iBand),'Color',bandColor(iBand,:),'Marker','.')
            myErrorbar((minInc:maxInc)*convFactTime,staticSeriesMeanMode(:,iBand),staticSeriesSemMode(:,iBand))
            legendEntries{iBand} = ['Band ' num2str(iBand)];
        end
        axis([minInc*convFactTime maxInc*convFactTime yMin yMax]);
        xlabel('Time from protruion onset (s)')
        ylabel(yAxisLabel)
        legend(legendEntries);
        for iBand = 1 : min(numBand,8)
            plot((minInc:maxInc)*convFactTime,dynamicSeriesBefMeanMode(:,iBand),'Color',bandColor(iBand,:),'Marker','.','LineStyle','--')
            myErrorbar((minInc:maxInc)*convFactTime,dynamicSeriesBefMeanMode(:,iBand),dynamicSeriesBefSemMode(:,iBand))
        end
        
    end %(if ~isnan(yMax))
    
end %(for jMode = 1 : numMode2plot)

%plot 3 - distances at the edge
subplot(numMode2plot+1,2,numMode2plot*2+1,'YDir','reverse'), hold on

yMin = 0;
yMax = 30*convFactDist;

%line indicating activity onset
plot([0 0],[yMin yMax],'k--')

%static series
plot((minInc:maxInc)*convFactTime,staticSeriesDist(:,1,1),'k','LineStyle','none')
myErrorbar((minInc:maxInc)*convFactTime,staticSeriesDist(:,1,1),staticSeriesDist(:,1,2))

%after dynamic
for iInc = 1 : 6
    plot((minInc:maxInc)*convFactTime,dynamicSeriesAftDist(:,iInc,1),'color',incColor(iInc,:),'LineStyle','none')
    myErrorbar((minInc:maxInc)*convFactTime,dynamicSeriesAftDist(:,iInc,1),dynamicSeriesAftDist(:,iInc,2))
end

%         %after dynamic combined
%         plot((minInc:maxInc)*convFactTime,combDynamicSeriesAftDist(:,1,1),'b','LineStyle','none')
%         myErrorbar((minInc:maxInc)*convFactTime,combDynamicSeriesAftDist(:,1,1),combDynamicSeriesAftDist(:,1,2))

axis([minInc*convFactTime maxInc*convFactTime yMin yMax]);
xlabel('Time from protruion onset (s)')
ylabel('Window distance from cell edge (um)')

%plot 4 - distances in bands
subplot(numMode2plot+1,2,numMode2plot*2+2,'YDir','reverse'), hold on

yMin = 0;
yMax = 100*convFactDist;

%line indicating activity onset
plot([0 0],[yMin yMax],'k--')

for iBand = 1 : min(numBand,8)
    plot((minInc:maxInc)*convFactTime,staticSeriesDist(:,iBand,1),'Color',bandColor(iBand,:),'LineStyle','none')
    myErrorbar((minInc:maxInc)*convFactTime,staticSeriesDist(:,iBand,1),staticSeriesDist(:,iBand,2))
end

axis([minInc*convFactTime maxInc*convFactTime yMin yMax]);
xlabel('Time from protruion onset (s)')
ylabel('Window distance from cell edge (um)')

if ~isempty(saveLoc)
    saveas(hFig,saveLoc)
end

