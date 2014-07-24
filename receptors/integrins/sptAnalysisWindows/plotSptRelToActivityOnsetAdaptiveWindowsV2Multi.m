    function plotSptRelToActivityOnsetAdaptiveWindowsV2Multi(particleBehavior,...
    windowDistFromEdge,mode2plot,minNP,figureName,saveLoc,convFact,...
    yAxisLabel,axisLimits,plotWinDist,winFigName,winFigLoc,condList)

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

if nargin < 13 || isempty(condList)
    condList = [];
end

%get number of conditions
numCond = length(particleBehavior);

%get increment range
for iCond = 1 : numCond
    minInc(iCond) = -size(particleBehavior(iCond).befStatic.mean,1);
    maxInc(iCond) = size(particleBehavior(iCond).aftStatic.mean,1);
end

%get number of bands and modes
for iCond = 1 : numCond
    [~,numBand(iCond),numMode(iCond)] = size(particleBehavior(iCond).onset.mean);
end

condColor = [0 0 0; 0 1 0; 0 0 1; 1 0 0; 0 1 1; 1 0 1];

%open figure
hFig = figure('Name',figureName);
hold on

for iCond = 1 : numCond
    
    %construct series of behavior to plot
    
    %static series
    staticSeriesMean = [particleBehavior(iCond).befStatic.mean(end:-1:1,:,:); ...
        particleBehavior(iCond).onset.mean; particleBehavior(iCond).aftStatic.mean] * convFactProp;
    staticSeriesStd = [particleBehavior(iCond).befStatic.std(end:-1:1,:,:); ...
        particleBehavior(iCond).onset.std; particleBehavior(iCond).aftStatic.std] * convFactProp;
    staticSeriesNP = [particleBehavior(iCond).befStatic.numPoints(end:-1:1,:,:); ...
        particleBehavior(iCond).onset.numPoints; particleBehavior(iCond).aftStatic.numPoints];
    staticSeriesNP(isnan(staticSeriesNP)) = 1;
    
    staticSeriesMean = staticSeriesMean(:,:,mode2plot);
    staticSeriesStd = staticSeriesStd(:,:,mode2plot);
    staticSeriesNP = staticSeriesNP(:,:,mode2plot);
    
    tmp = [windowDistFromEdge(iCond).befStatic(end:-1:1,:,:,1); ...
        windowDistFromEdge(iCond).onset(:,:,:,1); windowDistFromEdge(iCond).aftStatic(:,:,:,1)] ...
        * convFactDist;
    staticSeriesDistMid = (tmp(:,:,1) + tmp(:,:,2)) / 2;
    staticSeriesDistHR = (tmp(:,:,2) - tmp(:,:,1)) / 2;
    
    staticSeriesBulkDistMid = (tmp(:,4,1) + tmp(:,end,2)) / 2;
    staticSeriesBulkDistHR = (tmp(:,end,2) - tmp(:,4,1)) / 2;
    
    %dynamic series before
    dynamicSeriesBefMean = [particleBehavior(iCond).befDynamic.mean(end:-1:1,:,:); ...
        particleBehavior(iCond).onset.mean; NaN(maxInc(iCond),numBand(iCond),numMode(iCond))] * convFactProp;
    dynamicSeriesBefStd = [particleBehavior(iCond).befDynamic.std(end:-1:1,:,:); ...
        particleBehavior(iCond).onset.std; NaN(maxInc(iCond),numBand(iCond),numMode(iCond))] * convFactProp;
    dynamicSeriesBefNP = [particleBehavior(iCond).befDynamic.numPoints(end:-1:1,:,:); ...
        particleBehavior(iCond).onset.numPoints; ones(maxInc(iCond),numBand(iCond),numMode(iCond))];
    dynamicSeriesBefNP(isnan(dynamicSeriesBefNP)) = 1;
    
    dynamicSeriesBefMean = dynamicSeriesBefMean(:,:,mode2plot);
    dynamicSeriesBefStd = dynamicSeriesBefStd(:,:,mode2plot);
    dynamicSeriesBefNP = dynamicSeriesBefNP(:,:,mode2plot);
    
    tmp = [windowDistFromEdge(iCond).befDynamic(end:-1:1,:,:,1); ...
        windowDistFromEdge(iCond).onset(:,:,:,1); NaN(maxInc(iCond),numBand(iCond),2)] ...
        * convFactDist;
    dynamicSeriesBefDistMid = (tmp(:,:,1) + tmp(:,:,2)) / 2;
    dynamicSeriesBefDistHR = (tmp(:,:,2) - tmp(:,:,1)) / 2;
    
    %dynamic series after
    dynamicSeriesAftMean = [NaN(-minInc(iCond)+1,maxInc(iCond),numMode(iCond)); particleBehavior(iCond).aftDynamic.mean] ...
        * convFactProp;
    dynamicSeriesAftStd = [NaN(-minInc(iCond)+1,maxInc(iCond),numMode(iCond)); particleBehavior(iCond).aftDynamic.std] ...
        * convFactProp;
    dynamicSeriesAftNP = [ones(-minInc(iCond)+1,maxInc(iCond),numMode(iCond)); particleBehavior(iCond).aftDynamic.numPoints];
    dynamicSeriesAftNP(isnan(dynamicSeriesAftNP)) = 1;
    
    dynamicSeriesAftMean = dynamicSeriesAftMean(:,:,mode2plot);
    dynamicSeriesAftStd = dynamicSeriesAftStd(:,:,mode2plot);
    dynamicSeriesAftNP = dynamicSeriesAftNP(:,:,mode2plot);
    
    tmp = [NaN(-minInc(iCond)+1,maxInc(iCond),2); windowDistFromEdge(iCond).aftDynamic(:,:,:,1)] ...
        * convFactDist;
    dynamicSeriesAftDistMid = (tmp(:,:,1) + tmp(:,:,2)) / 2;
    dynamicSeriesAftDistHR = (tmp(:,:,2) - tmp(:,:,1)) / 2;
    
    %dynamic series right behind edge before and after
    edgeSeriesMean = [dynamicSeriesBefMean(1:-minInc(iCond)+1,1); diag(dynamicSeriesAftMean(-minInc(iCond)+2:end,:))];
    edgeSeriesStd = [dynamicSeriesBefStd(1:-minInc(iCond)+1,1); diag(dynamicSeriesAftStd(-minInc(iCond)+2:end,:))];
    edgeSeriesNP = [dynamicSeriesBefNP(1:-minInc(iCond)+1,1); diag(dynamicSeriesAftNP(-minInc(iCond)+2:end,:))];
    
    edgeSeriesDistMid = [dynamicSeriesBefDistMid(1:-minInc(iCond)+1,1); diag(dynamicSeriesAftDistMid(-minInc(iCond)+2:end,:))];
    edgeSeriesDistHR = [dynamicSeriesBefDistHR(1:-minInc(iCond)+1,1); diag(dynamicSeriesAftDistHR(-minInc(iCond)+2:end,:))];
    
    %combined dynamic series after
    combDynamicSeriesAftMean = [NaN(-minInc(iCond),1,numMode(iCond)); particleBehavior(iCond).aftDynamicComb.mean; ...
        NaN(1,1,numMode(iCond))] * convFactProp;
    combDynamicSeriesAftStd = [NaN(-minInc(iCond),1,numMode(iCond)); particleBehavior(iCond).aftDynamicComb.std; ...
        NaN(1,1,numMode(iCond))] * convFactProp;
    combDynamicSeriesAftNP = [ones(-minInc(iCond),1,numMode(iCond)); particleBehavior(iCond).aftDynamicComb.numPoints; ones(1,1,numMode(iCond))];
    combDynamicSeriesAftNP(isnan(combDynamicSeriesAftNP)) = 1;
    
    combDynamicSeriesAftMean = combDynamicSeriesAftMean(:,:,mode2plot);
    combDynamicSeriesAftStd = combDynamicSeriesAftStd(:,:,mode2plot);
    combDynamicSeriesAftNP = combDynamicSeriesAftNP(:,:,mode2plot);
    
    tmp = [NaN(-minInc(iCond),1,2); windowDistFromEdge(iCond).aftDynamicComb(:,:,:,1); ...
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
    
    % %define color vectors for plotting
    % %green to magenta in 11 steps
    % incColor = repmat((1:-0.1:0)',1,3).*repmat([0 1 0],11,1) + repmat((0:0.1:1)',1,3).*repmat([1 0 1],11,1);
    % incColor = min(incColor,1);
    % %black, blue, green, magenta
    % bandColor = [0 0 0; 0 0 1; 0 1 0; 1 0 1];
    %
    % %define maximum number of increments and bands to plot
    % maxInc(iCond)Plot = min(maxInc(iCond),size(incColor,1));
    % maxBandPlot = min(numBand(iCond),size(bandColor,1));
    
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
        xMin = minInc(iCond)*convFactTime;
        xMax = maxInc(iCond)*convFactTime;
    else
        xMin = axisLimits(1);
        xMax = axisLimits(2);
        yMin = axisLimits(3);
        yMax = axisLimits(4);
    end
    
    if ~isnan(yMax)
        
        %plot 1 - right behind edge, before and after
        subplot(1,2,1), hold on
        
        plot((minInc(iCond):maxInc(iCond))*convFactTime,edgeSeriesMean(:,1),'Color',condColor(iCond,:),'Marker','.')
        
        axis([xMin xMax yMin yMax]);
        xlabel('Time from protrusion onset (s)')
        ylabel(yAxisLabel)
        
        plot((minInc(iCond):maxInc(iCond))*convFactTime,edgeSeriesMean(:,1)+edgeSeriesSem(:,1),'Color',condColor(iCond,:),'LineStyle','--')
        plot((minInc(iCond):maxInc(iCond))*convFactTime,edgeSeriesMean(:,1)-edgeSeriesSem(:,1),'Color',condColor(iCond,:),'LineStyle','--')
        
        %plot 2 - at point of protrusion onset, before and after
        subplot(1,2,2), hold on
        
        plot((minInc(iCond):maxInc(iCond))*convFactTime,staticSeriesMean(:,1),'Color',condColor(iCond,:),'Marker','.')
        
        axis([xMin xMax yMin yMax]);
        xlabel('Time from protrusion onset (s)')
        ylabel(yAxisLabel)
        
        plot((minInc(iCond):maxInc(iCond))*convFactTime,staticSeriesMean(:,1)+staticSeriesSem(:,1),'Color',condColor(iCond,:),'LineStyle','--')
        plot((minInc(iCond):maxInc(iCond))*convFactTime,staticSeriesMean(:,1)-staticSeriesSem(:,1),'Color',condColor(iCond,:),'LineStyle','--')
        
        %     %plot 3 - at point of protrusion onset and new cell areas, after
        %     subplot(2,2,3), hold on
        %
        %     %point of protruion onset
        %     plot((0:maxInc(iCond))*convFactTime,staticSeriesMean(-minInc(iCond)+1:end,1),'k','Marker','.')
        %
        %     %new cell areas
        %     legendEntries = cell(1,maxInc(iCond)Plot);
        %     for iInc = maxIncPlot : -1 : 1
        %         plot((minInc(iCond):maxInc(iCond))*convFactTime,dynamicSeriesAftMean(:,iInc),'color',incColor(iInc,:),'Marker','.')
        %         legendEntries{maxIncPlot-iInc+1} = ['Dynamic after Inc ' num2str(iInc)];
        %     end
        %
        %     %right behind edge
        %     plot((0:maxInc(iCond))*convFactTime,edgeSeriesMean(-minInc(iCond)+1:end,1),'b','Marker','.')
        %
        %     axis([xMin xMax yMin yMax]);
        %     xlabel('Time from protrusion onset (s)')
        %     ylabel(yAxisLabel)
        %
        %     legendEntries = [{'Band 1'} legendEntries {'Right behind cell edge'}];
        %     legend(legendEntries);
        %
        %     plot((0:maxInc(iCond))*convFactTime,staticSeriesMean(-minInc(iCond)+1:end,1)+staticSeriesSem(-minInc(iCond)+1:end,1),'k--')
        %     plot((0:maxInc(iCond))*convFactTime,staticSeriesMean(-minInc(iCond)+1:end,1)-staticSeriesSem(-minInc(iCond)+1:end,1),'k--')
        %
        %     for iInc = 1 : maxIncPlot
        %         plot((minInc(iCond):maxInc(iCond))*convFactTime,dynamicSeriesAftMean(:,iInc)+dynamicSeriesAftSem(:,iInc),'color',incColor(iInc,:),'LineStyle','--')
        %         plot((minInc(iCond):maxInc(iCond))*convFactTime,dynamicSeriesAftMean(:,iInc)-dynamicSeriesAftSem(:,iInc),'color',incColor(iInc,:),'LineStyle','--')
        %     end
        %
        %     plot((0:maxInc(iCond))*convFactTime,edgeSeriesMean(-minInc(iCond)+1:end,1)+edgeSeriesSem(-minInc(iCond)+1:end,1),'color',[0.5 0.5 0.5],'LineStyle','--')
        %     plot((0:maxInc(iCond))*convFactTime,edgeSeriesMean(-minInc(iCond)+1:end,1)-edgeSeriesSem(-minInc(iCond)+1:end,1),'color',[0.5 0.5 0.5],'LineStyle','--')
        %
        %     %line indicating activity onset
        %     plot([0 0],[yMin yMax],'k:')
        %
        %     %plot 4 - values in bands
        %     subplot(2,2,4), hold on
        %
        %     plot((minInc(iCond):maxInc(iCond))*convFactTime,mean(staticSeriesMean(:,4:end),2),'Color',bandColor(4,:),'Marker','.')
        %     for iBand = 3 : -1 : 1
        %         plot((minInc(iCond):maxInc(iCond))*convFactTime,staticSeriesMean(:,iBand),'Color',bandColor(iBand,:),'Marker','.')
        %     end
        %
        %     axis([xMin xMax yMin yMax]);
        %     xlabel('Time from protrusion onset (s)')
        %     ylabel(yAxisLabel)
        %
        %     legend({'2-6 um at onset','Band 3 (1.33-2 um at onset)','Band 2 (0.67-1.33 um at onset)','Band 1 (0-0.67 um at onset)'});
        %
        %     for iBand = 1 : 3
        %         plot((minInc(iCond):maxInc(iCond))*convFactTime,staticSeriesMean(:,iBand)+staticSeriesSem(:,iBand),'Color',bandColor(iBand,:),'LineStyle','--')
        %         plot((minInc(iCond):maxInc(iCond))*convFactTime,staticSeriesMean(:,iBand)-staticSeriesSem(:,iBand),'Color',bandColor(iBand,:),'LineStyle','--')
        %     end
        %     plot((minInc(iCond):maxInc(iCond))*convFactTime,mean(staticSeriesMean(:,4:end),2)+mean(staticSeriesSem(:,4:end),2),'Color',bandColor(4,:),'LineStyle','--')
        %     plot((minInc(iCond):maxInc(iCond))*convFactTime,mean(staticSeriesMean(:,4:end),2)-mean(staticSeriesSem(:,4:end),2),'Color',bandColor(4,:),'LineStyle','--')
        %
        %     %line indicating activity onset
        %     plot([0 0],[yMin yMax],'k:')
        %
        %     %edge and window distances from edge
        %     if plotWinDist
        %
        %         %open figure
        %         hFig2 = figure('Name',winFigName);
        %         hold on
        %
        %         %get edge position over time
        %         edgePos = staticSeriesDistMid(:,1) - staticSeriesDistHR(:,1);
        %
        %         %shift window positions to follow moving edge instead of the
        %         %original assumption that edge position is always at zero
        %         staticSeriesDistMid = repmat(edgePos,1,numBand(iCond)) - staticSeriesDistMid;
        %         staticSeriesBulkDistMid = edgePos - staticSeriesBulkDistMid;
        %         dynamicSeriesBefDistMid = repmat(edgePos,1,numBand(iCond)) - dynamicSeriesBefDistMid;
        %         dynamicSeriesAftDistMid = repmat(edgePos,1,maxInc(iCond)) - dynamicSeriesAftDistMid;
        %         edgeSeriesDistMid = edgePos - edgeSeriesDistMid;
        %         combDynamicSeriesAftDistMid = edgePos - combDynamicSeriesAftDistMid;
        %
        %         yMin = min(staticSeriesDistMid(:,1)-staticSeriesDistHR(:,1)) - 0.5;
        %         yMax = max(edgePos) + 0.5;
        %
        %         %plot 1 - right behind edge, before and after
        %         subplot(2,2,1), hold on
        %
        %         plot((minInc(iCond):maxInc(iCond))*convFactTime,edgeSeriesDistMid,'k','LineStyle','none')
        %         myErrorbar((minInc(iCond):maxInc(iCond))*convFactTime,edgeSeriesDistMid,edgeSeriesDistHR)
        %
        %         plot([(minInc(iCond)-0.5:maxInc(iCond)-0.5)*convFactTime; (minInc(iCond)+0.5:maxInc(iCond)+0.5)*convFactTime],...
        %             repmat(edgePos',2,1),'r','LineWidth',2)
        %
        %         %line indicating activity onset
        %         plot([0 0],[yMin yMax],'k:')
        %
        %         axis([xMin xMax yMin yMax]);
        %         xlabel('Time from protrusion onset (s)')
        %         ylabel('Edge and window position (um)')
        %
        %         %plot 2 - at point of protrusion onset, before and after
        %         subplot(2,2,2), hold on
        %
        %         plot((minInc(iCond):maxInc(iCond))*convFactTime,staticSeriesDistMid(:,1),'k','LineStyle','none')
        %         myErrorbar((minInc(iCond):maxInc(iCond))*convFactTime,staticSeriesDistMid(:,1),staticSeriesDistHR(:,1))
        %
        %         plot([(minInc(iCond)-0.5:maxInc(iCond)-0.5)*convFactTime; (minInc(iCond)+0.5:maxInc(iCond)+0.5)*convFactTime],...
        %             repmat(edgePos',2,1),'r','LineWidth',2)
        %
        %         %line indicating activity onset
        %         plot([0 0],[yMin yMax],'k:')
        %
        %         axis([xMin xMax yMin yMax]);
        %         xlabel('Time from protrusion onset (s)')
        %         ylabel('Edge and window position (um)')
        %
        %         %plot 3 - at point of protrusion onset and new cell areas, after
        %         subplot(2,2,3), hold on
        %
        %         plot((0:maxInc(iCond))*convFactTime,staticSeriesDistMid(-minInc(iCond)+1:end,1),'k','LineStyle','none')
        %         myErrorbar((0:maxInc(iCond))*convFactTime,staticSeriesDistMid(-minInc(iCond)+1:end,1),staticSeriesDistHR(-minInc(iCond)+1:end,1))
        %
        %         for iInc = maxInc(iCond)Plot : -1 : 1
        %             plot((minInc(iCond):maxInc(iCond))*convFactTime,dynamicSeriesAftDistMid(:,iInc),'color',incColor(iInc,:),'LineStyle','none')
        %             myErrorbar((minInc(iCond):maxInc(iCond))*convFactTime,dynamicSeriesAftDistMid(:,iInc),dynamicSeriesAftDistHR(:,iInc))
        %         end
        %
        %         plot((0:maxInc(iCond))*convFactTime,edgeSeriesDistMid(-minInc(iCond)+1:end),'b','LineStyle','none')
        %         myErrorbar((0:maxInc(iCond))*convFactTime,edgeSeriesDistMid(-minInc(iCond)+1:end),edgeSeriesDistHR(-minInc(iCond)+1:end))
        %
        %         plot([(minInc(iCond)-0.5:maxInc(iCond)-0.5)*convFactTime; (minInc(iCond)+0.5:maxInc(iCond)+0.5)*convFactTime],...
        %             repmat(edgePos',2,1),'r','LineWidth',2)
        %
        %         %line indicating activity onset
        %         plot([0 0],[yMin yMax],'k:')
        %
        %         axis([xMin xMax yMin yMax]);
        %         xlabel('Time from protrusion onset (s)')
        %         ylabel('Edge and window position (um)')
        %
        %         %plot 4 - values in bands
        %         subplot(2,2,4), hold on
        %
        %         plot((minInc(iCond):maxInc(iCond))*convFactTime,staticSeriesBulkDistMid,'Color',bandColor(4,:),'LineStyle','none')
        %         myErrorbar((minInc(iCond):maxInc(iCond))*convFactTime,staticSeriesBulkDistMid,staticSeriesBulkDistHR);
        %         for iBand = 3 : -1 : 1
        %             plot((minInc(iCond):maxInc(iCond))*convFactTime,staticSeriesDistMid(:,iBand),'Color',bandColor(iBand,:),'LineStyle','none')
        %             myErrorbar((minInc(iCond):maxInc(iCond))*convFactTime,staticSeriesDistMid(:,iBand),staticSeriesDistHR(:,iBand));
        %         end
        %
        %         plot([(minInc(iCond)-0.5:maxInc(iCond)-0.5)*convFactTime; (minInc(iCond)+0.5:maxInc(iCond)+0.5)*convFactTime],...
        %             repmat(edgePos',2,1),'r','LineWidth',2)
        %
        %         yMin = min(staticSeriesBulkDistMid-staticSeriesBulkDistHR) - 0.5;
        %
        %         %line indicating activity onset
        %         plot([0 0],[yMin yMax],'k:')
        %
        %         axis([xMin xMax yMin yMax]);
        %         xlabel('Time from protrusion onset (s)')
        %         ylabel('Edge and window position (um)')
        %
        %     end
        
    end %(if ~isnan(yMax))
    
end %(for iCond = 1 : numCond)

%plot 1 - right behind edge, before and after
subplot(1,2,1), hold on

%line indicating activity onset
plot([0 0],[yMin yMax],'k:')

%legend
legendList = repmat(condList,3,1);
legend(legendList(:))

%plot 2 - at point of protrusion onset, before and after
subplot(1,2,2), hold on

%line indicating activity onset
plot([0 0],[yMin yMax],'k:')

%legend
legendList = repmat(condList,3,1);
legend(legendList(:))

if ~isempty(saveLoc)
    saveas(hFig,saveLoc)
end
if ~isempty(winFigLoc)
    saveas(hFig,winFigLoc)
end

