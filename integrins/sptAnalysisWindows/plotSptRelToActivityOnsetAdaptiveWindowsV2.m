    function plotSptRelToActivityOnsetAdaptiveWindowsV2(particleBehavior,...
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
    
    %     plot((minInc:maxInc)*convFactTime,edgeSeriesMean(:,1),'k','Marker','.')
    
    axis([xMin xMax yMin yMax]);
    xlabel('Time from protrusion onset (s)')
    ylabel(yAxisLabel)
    
    %     legend('Right behind cell edge');
    %
    %     plot((minInc:maxInc)*convFactTime,edgeSeriesMean(:,1)+edgeSeriesSem(:,1),'k--')
    %     plot((minInc:maxInc)*convFactTime,edgeSeriesMean(:,1)-edgeSeriesSem(:,1),'k--')
    
    plot((minInc:maxInc)*convFactTime,edgeSeriesMean(:,1)+edgeSeriesSem(:,1),'k')
    legend('Right behind cell edge');
    plot((minInc:maxInc)*convFactTime,edgeSeriesMean(:,1)-edgeSeriesSem(:,1),'k')
    
    %line indicating activity onset
    plot([0 0],[yMin yMax],'k:')
    
    %plot 2 - at point of protrusion onset, before and after
    subplot(2,2,2), hold on
    
    %     plot((minInc:maxInc)*convFactTime,staticSeriesMean(:,1),'k','Marker','.')
    
    axis([xMin xMax yMin yMax]);
    xlabel('Time from protrusion onset (s)')
    ylabel(yAxisLabel)
    
    %     legend('Band 1 (0-0.67 um at onset)');
    %
    %     plot((minInc:maxInc)*convFactTime,staticSeriesMean(:,1)+staticSeriesSem(:,1),'k--')
    %     plot((minInc:maxInc)*convFactTime,staticSeriesMean(:,1)-staticSeriesSem(:,1),'k--')
    
    plot((minInc:maxInc)*convFactTime,staticSeriesMean(:,1)+staticSeriesSem(:,1),'k')
    legend('Band 1 (0-0.67 um at onset)');
    plot((minInc:maxInc)*convFactTime,staticSeriesMean(:,1)-staticSeriesSem(:,1),'k')
    
    %line indicating activity onset
    plot([0 0],[yMin yMax],'k:')
    
    %plot 3 - at point of protrusion onset and new cell areas, after
    subplot(2,2,3), hold on
    
    %     %point of protruion onset
    %     plot((0:maxInc)*convFactTime,staticSeriesMean(-minInc+1:end,1),'k','Marker','.')
    %
    %     %new cell areas
    %     legendEntries = cell(1,maxIncPlot);
    %     for iInc = maxIncPlot : -1 : 1
    %         plot((minInc:maxInc)*convFactTime,dynamicSeriesAftMean(:,iInc),'color',incColor(iInc,:),'Marker','.')
    %         legendEntries{maxIncPlot-iInc+1} = ['Dynamic after Inc ' num2str(iInc)];
    %     end
    %
    %     %right behind edge
    %     plot((0:maxInc)*convFactTime,edgeSeriesMean(-minInc+1:end,1),'b','Marker','.')
    
    axis([xMin xMax yMin yMax]);
    xlabel('Time from protrusion onset (s)')
    ylabel(yAxisLabel)
    
    %     legendEntries = [{'Band 1'} legendEntries {'Right behind cell edge'}];
    %     legend(legendEntries);
    %
    %     plot((0:maxInc)*convFactTime,staticSeriesMean(-minInc+1:end,1)+staticSeriesSem(-minInc+1:end,1),'k--')
    %     plot((0:maxInc)*convFactTime,staticSeriesMean(-minInc+1:end,1)-staticSeriesSem(-minInc+1:end,1),'k--')
    %
    %     for iInc = 1 : maxIncPlot
    %         plot((minInc:maxInc)*convFactTime,dynamicSeriesAftMean(:,iInc)+dynamicSeriesAftSem(:,iInc),'color',incColor(iInc,:),'LineStyle','--')
    %         plot((minInc:maxInc)*convFactTime,dynamicSeriesAftMean(:,iInc)-dynamicSeriesAftSem(:,iInc),'color',incColor(iInc,:),'LineStyle','--')
    %     end
    %
    %     plot((0:maxInc)*convFactTime,edgeSeriesMean(-minInc+1:end,1)+edgeSeriesSem(-minInc+1:end,1),'color',[0.5 0.5 0.5],'LineStyle','--')
    %     plot((0:maxInc)*convFactTime,edgeSeriesMean(-minInc+1:end,1)-edgeSeriesSem(-minInc+1:end,1),'color',[0.5 0.5 0.5],'LineStyle','--')
   
    plot((0:maxInc)*convFactTime,staticSeriesMean(-minInc+1:end,1)+staticSeriesSem(-minInc+1:end,1),'k')
    legendEntries = cell(1,maxIncPlot);
    for iInc = 1 : maxIncPlot
        plot((minInc:maxInc)*convFactTime,dynamicSeriesAftMean(:,iInc)+dynamicSeriesAftSem(:,iInc),'color',incColor(iInc,:))
        legendEntries{maxIncPlot-iInc+1} = ['Dynamic after Inc ' num2str(iInc)];
    end
    plot((0:maxInc)*convFactTime,edgeSeriesMean(-minInc+1:end,1)+edgeSeriesSem(-minInc+1:end,1),'b')

    legendEntries = [{'Band 1'} legendEntries {'Right behind cell edge'}];
    legend(legendEntries);
    
    plot((0:maxInc)*convFactTime,staticSeriesMean(-minInc+1:end,1)-staticSeriesSem(-minInc+1:end,1),'k')
    for iInc = 1 : maxIncPlot
        plot((minInc:maxInc)*convFactTime,dynamicSeriesAftMean(:,iInc)-dynamicSeriesAftSem(:,iInc),'color',incColor(iInc,:))
    end
    plot((0:maxInc)*convFactTime,edgeSeriesMean(-minInc+1:end,1)-edgeSeriesSem(-minInc+1:end,1),'b')
    
    %line indicating activity onset
    plot([0 0],[yMin yMax],'k:')
    
    %plot 4 - values in bands
    subplot(2,2,4), hold on
    
    %     plot((minInc:maxInc)*convFactTime,mean(staticSeriesMean(:,4:end),2),'Color',bandColor(4,:),'Marker','.')
    %     for iBand = 3 : -1 : 1
    %         plot((minInc:maxInc)*convFactTime,staticSeriesMean(:,iBand),'Color',bandColor(iBand,:),'Marker','.')
    %     end
    
    axis([xMin xMax yMin yMax]);
    xlabel('Time from protrusion onset (s)')
    ylabel(yAxisLabel)
    
    %     legend({'2-6 um at onset','Band 3 (1.33-2 um at onset)','Band 2 (0.67-1.33 um at onset)','Band 1 (0-0.67 um at onset)'});
    %
    %     for iBand = 1 : 3
    %         plot((minInc:maxInc)*convFactTime,staticSeriesMean(:,iBand)+staticSeriesSem(:,iBand),'Color',bandColor(iBand,:),'LineStyle','--')
    %         plot((minInc:maxInc)*convFactTime,staticSeriesMean(:,iBand)-staticSeriesSem(:,iBand),'Color',bandColor(iBand,:),'LineStyle','--')
    %     end
    %     plot((minInc:maxInc)*convFactTime,mean(staticSeriesMean(:,4:end),2)+mean(staticSeriesSem(:,4:end),2),'Color',bandColor(4,:),'LineStyle','--')
    %     plot((minInc:maxInc)*convFactTime,mean(staticSeriesMean(:,4:end),2)-mean(staticSeriesSem(:,4:end),2),'Color',bandColor(4,:),'LineStyle','--')

    for iBand = 1 : 3
        plot((minInc:maxInc)*convFactTime,staticSeriesMean(:,iBand)+staticSeriesSem(:,iBand),'Color',bandColor(iBand,:))
    end
    plot((minInc:maxInc)*convFactTime,mean(staticSeriesMean(:,4:end),2)+mean(staticSeriesSem(:,4:end),2),'Color',bandColor(4,:))

    legend({'2-6 um at onset','Band 3 (1.33-2 um at onset)','Band 2 (0.67-1.33 um at onset)','Band 1 (0-0.67 um at onset)'});
    
    for iBand = 1 : 3
        plot((minInc:maxInc)*convFactTime,staticSeriesMean(:,iBand)-staticSeriesSem(:,iBand),'Color',bandColor(iBand,:))
    end
    plot((minInc:maxInc)*convFactTime,mean(staticSeriesMean(:,4:end),2)-mean(staticSeriesSem(:,4:end),2),'Color',bandColor(4,:))
    
    %line indicating activity onset
    plot([0 0],[yMin yMax],'k:')
    
    %edge and window distances from edge
    if plotWinDist
        
        %open figure
        hFig2 = figure('Name',winFigName);
        hold on
        
        %shift window positions to follow moving edge instead of the
        %original assumption that edge position is always at zero
        staticSeriesDistMid = repmat(edgePos,1,numBand) - staticSeriesDistMid;
        staticSeriesBulkDistMid = edgePos - staticSeriesBulkDistMid;
        dynamicSeriesBefDistMid = repmat(edgePos,1,numBand) - dynamicSeriesBefDistMid;
        dynamicSeriesAftDistMid = repmat(edgePos,1,maxInc) - dynamicSeriesAftDistMid;
        edgeSeriesDistMid = edgePos - edgeSeriesDistMid;
        combDynamicSeriesAftDistMid = edgePos - combDynamicSeriesAftDistMid;
        
        yMin = min(staticSeriesDistMid(:,1)-staticSeriesDistHR(:,1)) - 0.5;
        yMax = max(edgePos) + 0.5;
        
        %plot 1 - right behind edge, before and after
        subplot(2,2,1), hold on
        
        %         plot((minInc:maxInc)*convFactTime,edgeSeriesDistMid,'k','LineStyle','none')
        %         myErrorbar((minInc:maxInc)*convFactTime,edgeSeriesDistMid,edgeSeriesDistHR)
        plot((minInc:maxInc)*convFactTime,edgeSeriesDistMid+edgeSeriesDistHR,'k')
        plot((minInc:maxInc)*convFactTime,edgeSeriesDistMid-edgeSeriesDistHR,'k')
        
        %         plot([(minInc-0.5:maxInc-0.5)*convFactTime; (minInc+0.5:maxInc+0.5)*convFactTime],...
        %             repmat(edgePos',2,1),'r','LineWidth',2)
        plot((minInc:maxInc)*convFactTime,edgePos,'r')
        
        %line indicating activity onset
        plot([0 0],[yMin yMax],'k:')
        
        axis([xMin xMax yMin yMax]);
        xlabel('Time from protrusion onset (s)')
        ylabel('Edge and window position (um)')
        
        %plot 2 - at point of protrusion onset, before and after
        subplot(2,2,2), hold on
        
        %         plot((minInc:maxInc)*convFactTime,staticSeriesDistMid(:,1),'k','LineStyle','none')
        %         myErrorbar((minInc:maxInc)*convFactTime,staticSeriesDistMid(:,1),staticSeriesDistHR(:,1))
        plot((minInc:maxInc)*convFactTime,staticSeriesDistMid(:,1)+staticSeriesDistHR(:,1),'k')
        plot((minInc:maxInc)*convFactTime,staticSeriesDistMid(:,1)-staticSeriesDistHR(:,1),'k')
        
        %         plot([(minInc-0.5:maxInc-0.5)*convFactTime; (minInc+0.5:maxInc+0.5)*convFactTime],...
        %             repmat(edgePos',2,1),'r','LineWidth',2)
        plot((minInc:maxInc)*convFactTime,edgePos,'r')
        
        %line indicating activity onset
        plot([0 0],[yMin yMax],'k:')
        
        axis([xMin xMax yMin yMax]);
        xlabel('Time from protrusion onset (s)')
        ylabel('Edge and window position (um)')
        
        %plot 3 - at point of protrusion onset and new cell areas, after
        subplot(2,2,3), hold on
        
        %         plot((0:maxInc)*convFactTime,staticSeriesDistMid(-minInc+1:end,1),'k','LineStyle','none')
        %         myErrorbar((0:maxInc)*convFactTime,staticSeriesDistMid(-minInc+1:end,1),staticSeriesDistHR(-minInc+1:end,1))
        plot((0:maxInc)*convFactTime,staticSeriesDistMid(-minInc+1:end,1)+staticSeriesDistHR(-minInc+1:end,1),'k')
        plot((0:maxInc)*convFactTime,staticSeriesDistMid(-minInc+1:end,1)-staticSeriesDistHR(-minInc+1:end,1),'k')
        
        for iInc = maxIncPlot : -1 : 1
            %             plot((minInc:maxInc)*convFactTime,dynamicSeriesAftDistMid(:,iInc),'color',incColor(iInc,:),'LineStyle','none')
            %             myErrorbar((minInc:maxInc)*convFactTime,dynamicSeriesAftDistMid(:,iInc),dynamicSeriesAftDistHR(:,iInc))
            plot((minInc:maxInc)*convFactTime,dynamicSeriesAftDistMid(:,iInc)+dynamicSeriesAftDistHR(:,iInc),'color',incColor(iInc,:))
            plot((minInc:maxInc)*convFactTime,dynamicSeriesAftDistMid(:,iInc)-dynamicSeriesAftDistHR(:,iInc),'color',incColor(iInc,:))
        end
        
        %         plot((0:maxInc)*convFactTime,edgeSeriesDistMid(-minInc+1:end),'b','LineStyle','none')
        %         myErrorbar((0:maxInc)*convFactTime,edgeSeriesDistMid(-minInc+1:end),edgeSeriesDistHR(-minInc+1:end))
        plot((0:maxInc)*convFactTime,edgeSeriesDistMid(-minInc+1:end)+edgeSeriesDistHR(-minInc+1:end),'b')
        plot((0:maxInc)*convFactTime,edgeSeriesDistMid(-minInc+1:end)-edgeSeriesDistHR(-minInc+1:end),'b')
        
        %         plot([(minInc-0.5:maxInc-0.5)*convFactTime; (minInc+0.5:maxInc+0.5)*convFactTime],...
        %             repmat(edgePos',2,1),'r','LineWidth',2)
        plot((minInc:maxInc)*convFactTime,edgePos,'r')
        
        %line indicating activity onset
        plot([0 0],[yMin yMax],'k:')
        
        axis([xMin xMax yMin yMax]);
        xlabel('Time from protrusion onset (s)')
        ylabel('Edge and window position (um)')
        
        %plot 4 - values in bands
        subplot(2,2,4), hold on
        
        %         plot((minInc:maxInc)*convFactTime,staticSeriesBulkDistMid,'Color',bandColor(4,:),'LineStyle','none')
        %         myErrorbar((minInc:maxInc)*convFactTime,staticSeriesBulkDistMid,staticSeriesBulkDistHR);
        plot((minInc:maxInc)*convFactTime,staticSeriesBulkDistMid+staticSeriesBulkDistHR,'Color',bandColor(4,:))
        plot((minInc:maxInc)*convFactTime,staticSeriesBulkDistMid-staticSeriesBulkDistHR,'Color',bandColor(4,:))
        for iBand = 3 : -1 : 1
            %             plot((minInc:maxInc)*convFactTime,staticSeriesDistMid(:,iBand),'Color',bandColor(iBand,:),'LineStyle','none')
            %             myErrorbar((minInc:maxInc)*convFactTime,staticSeriesDistMid(:,iBand),staticSeriesDistHR(:,iBand));
            plot((minInc:maxInc)*convFactTime,staticSeriesDistMid(:,iBand)+staticSeriesDistHR(:,iBand),'Color',bandColor(iBand,:))
            plot((minInc:maxInc)*convFactTime,staticSeriesDistMid(:,iBand)-staticSeriesDistHR(:,iBand),'Color',bandColor(iBand,:))
        end
        
        %         plot([(minInc-0.5:maxInc-0.5)*convFactTime; (minInc+0.5:maxInc+0.5)*convFactTime],...
        %             repmat(edgePos',2,1),'r','LineWidth',2)
        plot((minInc:maxInc)*convFactTime,edgePos,'r')
        
        yMin = min(staticSeriesBulkDistMid-staticSeriesBulkDistHR) - 0.5;

        %line indicating activity onset
        plot([0 0],[yMin yMax],'k:')
        
        axis([xMin xMax yMin yMax]);
        xlabel('Time from protrusion onset (s)')
        ylabel('Edge and window position (um)')
        
    end
    
end %(if ~isnan(yMax))

if ~isempty(saveLoc)
    saveas(hFig,saveLoc)
end
if ~isempty(winFigLoc)
    saveas(hFig2,winFigLoc)
end

