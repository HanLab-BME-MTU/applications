    function plotSptRelToActivityOnsetAdaptiveWindowsV2(particleBehavior,...
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

%dynamic series right behind edge before and after
edgeSeriesMean = [dynamicSeriesBefMean(1:-minInc+1,1); diag(dynamicSeriesAftMean(-minInc+2:end,:))];
edgeSeriesStd = [dynamicSeriesBefStd(1:-minInc+1,1); diag(dynamicSeriesAftStd(-minInc+2:end,:))];
edgeSeriesNP = [dynamicSeriesBefNP(1:-minInc+1,1); diag(dynamicSeriesAftNP(-minInc+2:end,:))];

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

edgeSeriesMean(edgeSeriesNP<minNP) = NaN;
edgeSeriesStd(edgeSeriesNP<minNP) = NaN;

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
edgeSeriesSem = edgeSeriesStd./sqrt(edgeSeriesNP);

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
    
    %plot 1 - right behind edge, before and after
    subplot(2,2,1), hold on
    
    plot((minInc:maxInc)*convFactTime,edgeSeriesMean(:,1),'k','Marker','.')
    %     myErrorbar((minInc:maxInc)*convFactTime,edgeSeriesMean(:,1),edgeSeriesSem(:,1))
    
    axis([xMin xMax yMin yMax]);
    xlabel('Time from protrusion onset (s)')
    ylabel(yAxisLabel)
    
    legend('Right behind cell edge');
    
    plot((minInc:maxInc)*convFactTime,edgeSeriesMean(:,1)+edgeSeriesSem(:,1),'k:')
    plot((minInc:maxInc)*convFactTime,edgeSeriesMean(:,1)-edgeSeriesSem(:,1),'k:')

    %line indicating activity onset
    plot([0 0],[yMin yMax],'k:')
    
    %plot 2 - at point of protrusion onset, before and after
    subplot(2,2,2), hold on
    
    plot((minInc:maxInc)*convFactTime,staticSeriesMean(:,1),'k','Marker','.')
    %     myErrorbar((minInc:maxInc)*convFactTime,staticSeriesMean(:,1),staticSeriesSem(:,1))
    
    axis([xMin xMax yMin yMax]);
    xlabel('Time from protrusion onset (s)')
    ylabel(yAxisLabel)
    
    legend('Point of protrusion onset');
    
    plot((minInc:maxInc)*convFactTime,staticSeriesMean(:,1)+staticSeriesSem(:,1),'k:')
    plot((minInc:maxInc)*convFactTime,staticSeriesMean(:,1)-staticSeriesSem(:,1),'k:')

    %line indicating activity onset
    plot([0 0],[yMin yMax],'k:')
    
    %plot 3 - at point of protrusion onset and new cell areas, after
    subplot(2,2,3), hold on
    
    %point of protruion onset
    plot((0:maxInc)*convFactTime,staticSeriesMean(-minInc+1:end,1),'k','Marker','.')
    %     myErrorbar((0:maxInc)*convFactTime,staticSeriesMean(-minInc+1:end,1),staticSeriesSem(-minInc+1:end,1))
    
    %new cell areas
    legendEntries = cell(1,maxIncPlot);
    for iInc = 1 : maxIncPlot
        plot((minInc:maxInc)*convFactTime,dynamicSeriesAftMean(:,iInc),'color',incColor(iInc,:),'Marker','.')
        %         myErrorbar((minInc:maxInc)*convFactTime,dynamicSeriesAftMean(:,iInc),dynamicSeriesAftSem(:,iInc))
        legendEntries{iInc} = ['Dynamic after Inc ' num2str(iInc)];
    end
    
    %right behind edge
    plot((0:maxInc)*convFactTime,edgeSeriesMean(-minInc+1:end,1),'color',[0.5 0.5 0.5],'Marker','.')
    %     myErrorbar((0:maxInc)*convFactTime,edgeSeriesMean(-minInc+1:end,1),edgeSeriesSem(-minInc+1:end,1))
    
    axis([xMin xMax yMin yMax]);
    xlabel('Time from protrusion onset (s)')
    ylabel(yAxisLabel)
    
    legendEntries = [{'Point of protrusion onset'} legendEntries {'Behind cell edge'}];
    legend(legendEntries);
    
    plot((0:maxInc)*convFactTime,staticSeriesMean(-minInc+1:end,1)+staticSeriesSem(-minInc+1:end,1),'k:')
    plot((0:maxInc)*convFactTime,staticSeriesMean(-minInc+1:end,1)-staticSeriesSem(-minInc+1:end,1),'k:')

    for iInc = 1 : maxIncPlot
        plot((minInc:maxInc)*convFactTime,dynamicSeriesAftMean(:,iInc)+dynamicSeriesAftSem(:,iInc),'color',incColor(iInc,:),'LineStyle',':')
        plot((minInc:maxInc)*convFactTime,dynamicSeriesAftMean(:,iInc)-dynamicSeriesAftSem(:,iInc),'color',incColor(iInc,:),'LineStyle',':')
    end
    
    plot((0:maxInc)*convFactTime,edgeSeriesMean(-minInc+1:end,1)+edgeSeriesSem(-minInc+1:end,1),'color',[0.5 0.5 0.5],'LineStyle',':')
    plot((0:maxInc)*convFactTime,edgeSeriesMean(-minInc+1:end,1)-edgeSeriesSem(-minInc+1:end,1),'color',[0.5 0.5 0.5],'LineStyle',':')
   
    %line indicating activity onset
    plot([0 0],[yMin yMax],'k:')
    
    %plot 4 - values in bands
    subplot(2,2,4), hold on
    
    for iBand = 1 : 3
        plot((minInc:maxInc)*convFactTime,staticSeriesMean(:,iBand),'Color',bandColor(iBand,:),'Marker','.')
        %         myErrorbar((minInc:maxInc)*convFactTime,staticSeriesMean(:,iBand),staticSeriesSem(:,iBand))
    end
    plot((minInc:maxInc)*convFactTime,mean(staticSeriesMean(:,4:end),2),'Color',bandColor(4,:),'Marker','.')
    %     myErrorbar((minInc:maxInc)*convFactTime,mean(staticSeriesMean(:,4:end),2),mean(staticSeriesSem(:,4:end),2))
    
    axis([xMin xMax yMin yMax]);
    xlabel('Time from protrusion onset (s)')
    ylabel(yAxisLabel)
    
    legend({'Point of protrusion onset (0-0.67 um at onset)','Band 2 (0.67-1.33 um at onset)','Band 3 (1.33-2 um at onset)','> 2 um at onset'});
    
    for iBand = 1 : 3
        plot((minInc:maxInc)*convFactTime,staticSeriesMean(:,iBand)+staticSeriesSem(:,iBand),'Color',bandColor(iBand,:),'LineStyle',':')
        plot((minInc:maxInc)*convFactTime,staticSeriesMean(:,iBand)-staticSeriesSem(:,iBand),'Color',bandColor(iBand,:),'LineStyle',':')
    end
    plot((minInc:maxInc)*convFactTime,mean(staticSeriesMean(:,4:end),2)+mean(staticSeriesSem(:,4:end),2),'Color',bandColor(4,:),'LineStyle',':')
    plot((minInc:maxInc)*convFactTime,mean(staticSeriesMean(:,4:end),2)-mean(staticSeriesSem(:,4:end),2),'Color',bandColor(4,:),'LineStyle',':')

    %line indicating activity onset
    plot([0 0],[yMin yMax],'k:')
        
end %(if ~isnan(yMax))

if ~isempty(saveLoc)
    saveas(hFig,saveLoc)
end

