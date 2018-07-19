function [h] = gcaPlotIndividualProcessedTimeSeries(analysisResults,varargin)
%gcaPlotIndividualProcessedTimeSeries 
% 
% analysisResults : (struct) output of marcos edgeVelocityAnalysis which
%               finds significant protrusions and retraction events above noise
%               for each velocity time series. 
% 
% 
% windsToPlot : (numeric) 
% 
% OutputDirectory 
% 
% 
%% Input Parser 
ip = inputParser;

ip.CaseSensitive = false;
ip.KeepUnmatched = true;

ip.addRequired('analysisResults');
ip.addParameter('windsToPlot',1);
ip.addParameter('OutputDirectory',[]); 
ip.addParameter('saveResults',true); 
ip.addParameter('yLims',[]); 

ip.parse('analysisResults',varargin{:});
% 
windsToPlot = ip.Results.windsToPlot; 

%% 
if ~isempty(ip.Results.OutputDirectory); 
outputDirectory = ip.Results.OutputDirectory; 
else 
    outputDirectory = [pwd filesep 'ProcessedIndividualTimeSeries'];
    if ~isdir(outputDirectory)
        mkdir(outputDirectory) ; 
    end 
    
end 

for iWind = 1:length(windsToPlot)
    h(iWind) = setAxis('on'); 
    % TS = analysisResults.procExcEdgeMotion{windsToPlot};
    TS = analysisResults.data.rawTimeSeries(windsToPlot(iWind),:);
    TS = TS.*analysisResults.data.scaling;
    winIntC = analysisResults.data.winInterval{windsToPlot(iWind)};
    TS = TS(winIntC);
    time = 1:length(TS);
    time = (time.*analysisResults.data.timeInterval)-analysisResults.data.timeInterval;
    plot(time,TS,'Color','k','LineWidth',2);
    
    protBlock = analysisResults.protrusionAnalysis.windows(windsToPlot(iWind)).blockOut; 
    retBlock = analysisResults.retractionAnalysis.windows(windsToPlot(iWind)).blockOut; 

%     hold on
%     protBlock = analysisResults.protrusionAnalysis.blockOut; 
%     retBlock =  analysisResults.retractionAnalysis.blockOut; 
%     
    
    nProtB = numel(protBlock);
    nRetrB = numel(retBlock);
    
    for i = 1:nProtB
        plot(protBlock{i}.*analysisResults.data.timeInterval-analysisResults.data.timeInterval,TS(protBlock{i}),'r','LineWidth',2)
        scatter(protBlock{i}.*analysisResults.data.timeInterval-analysisResults.data.timeInterval,TS(protBlock{i}),50,'r','filled')
    end
    for i = 1:nRetrB
        plot(retBlock{i}.*analysisResults.data.timeInterval-analysisResults.data.timeInterval,TS(retBlock{i}),'b','LineWidth',2)
         scatter(retBlock{i}.*analysisResults.data.timeInterval-analysisResults.data.timeInterval,TS(retBlock{i}),50,'b','filled')
    end

limitProt = analysisResults.protrusionAnalysis.windows(windsToPlot(iWind)).limit; 
limitRet = analysisResults.retractionAnalysis.windows(windsToPlot(iWind)).limit; 

% Plot the limits from the EMD (empirical mode decomposition test) 
line([0,time(end)],[limitProt,limitProt],'LineWidth',2,'LineStyle','--','Color','k'); 
line([0,time(end)],[limitRet,limitRet],'LineWidth',2,'LineStyle','--','Color','k'); 
line([0,time(end)],[0,0],'Linewidth',2,'Color','k'); 
% plot the outlier data. 
procTimeSeries = analysisResults.data.procExcEdgeMotion{windsToPlot(iWind)}; 
outliersx = find(isnan(procTimeSeries)); 
outliersy = TS(outliersx); 
scatter(outliersx.*analysisResults.data.timeInterval-analysisResults.data.timeInterval,outliersy,100,'k','x'); 
xlabel('Time (s)'); 
ylabel('Velocity (nm/sec)'); 

title(['EMD Signal Detection for Window ' num2str(windsToPlot(iWind))]); 
if ~isempty(ip.Results.yLims)
    axis([time(1),time(end),ip.Results.yLims(1),ip.Results.yLims(2)]); 
      
end 
if ip.Results.saveResults
    saveas(h(iWind),[outputDirectory filesep 'Window_Number' num2str(windsToPlot(iWind),'%03d'),'.png']);  
    saveas(h(iWind),[outputDirectory filesep 'Window_Number' num2str(windsToPlot(iWind),'%03d') '.eps'],'psc2'); 

end

end % iWind 
