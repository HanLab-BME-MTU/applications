function [avgFreqAll,WindowsSampled] = plotMaxValuesInMaps(outEdgeL1,outSpeedL1,curPath)
%plotMaxValuesInMaps(outEdgeL1,outSpeedL1) plots time series in the
%activity map where its max edge velocity is maximum.
%   Detailed explanation goes here
%% setting
savePath = [curPath filesep 'FlowAndEdge'];
if ~exist(savePath,'dir')
    mkdir(savePath)
end
%% Sample the windows that experience the protrusion and retraction relatively equally
protrudingWindows = sum(outEdgeL1>0,2);
LT = size(outEdgeL1,2);
% Will use the windows that experience at least 20% of protrusion
thresProt = 0.3;
protrudingWindFraction = protrudingWindows/LT;
windsOfInterest = protrudingWindFraction>thresProt;
disp(['Protruding window fraction: ' num2str(sum(windsOfInterest)/length(windsOfInterest))])
maxOutCols = quantile(outEdgeL1,0.9,2);
% meanOutCols = mean(outEdgeL1,2);
% h=histogram(meanOutCols(protrudingWindows),100);
h=histogram(maxOutCols(windsOfInterest),100);
% a = gca;
% a.YScale = 'log'

%% Go through each window
% maxOutCols = max(outEdgeL1,[],2);
% [~,iMaxOutOne] = max(maxOutCols);
p=0;
avgFreqAll = zeros(sum(windsOfInterest),1);
for ii = find(windsOfInterest)'
    p=p+1;
    hf = figure(100);
    hf.Units='inches';
    hf.Position(3:4) = [2 3];
    % plot
    subplot(3,1,1),plot(outEdgeL1(ii,:)),title('Edge Speed')
    subplot(3,1,2),plot(outSpeedL1(ii,:)), title('Flow Vel')

    x = outSpeedL1(ii,:);
    x = fillmissing(x,'linear');
    n = 2^nextpow2(length(x));
    y = fft(x-mean(x),n);
    tInterval = 6; %sec
    fs = 1/tInterval; % sample frequency (Hz)
    f = (0:n-1)*(fs/n); %frequency range
%     power = abs(y).^2/n; % power of DFT
    P2 = abs(y/n); % power of DFT
    P1 = P2(1:n/2+1);
    % P1 = power(1:n/2);
    P1(2:end) = 2*P1(2:end);
    f1 = f(1:n/2+1);

    subplot(3,1,3),
    plot(f1,P1,"-o")
    xlabel('Frequency (Hz)')
    ylabel('|P1(f)|')
    title(['Single-sided spectrum ' num2str(ii) 'th window'])
    avgFreq = mean(P1*f1'/sum(P1));
    line([avgFreq avgFreq],[0 max(P1)],'LineStyle',':')
    text(avgFreq,max(P1),['Avg freq: ' num2str(avgFreq,2) ' Hz'],'FontSize',5)
    disp(['The average frequency: ' num2str(avgFreq) ' Hz in ' num2str(ii) 'th window.'])
    avgFreqAll(p) = avgFreq;
    [~,iMax] = max(P1);
    peakFreq = f1(iMax);
    disp(['The peak frequency: ' num2str(peakFreq) ' Hz'])
    savefig(hf,[curPath filesep 'edgeFlowFreq' num2str(ii) '.fig'])
    saveas(hf,[curPath filesep 'edgeFlowFreq' num2str(ii) '.png'])
end
WindowsSampled = find(windsOfInterest);
disp(['The avg of the average frequencies: ' num2str(mean(avgFreqAll)) ' Hz'])
disp(['The std of the average frequencies: ' num2str(std(avgFreqAll)) ' Hz'])
end

