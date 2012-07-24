function plotParticleBehaviorRelToActivityOnset(particleBehavior,figureName,saveLoc)

if nargin < 2 || isempty(figureName)
    figureName = 'test';
end

if nargin < 3 || isempty(saveLoc)
    saveLoc = [];
end

%get mean, standard deviation and number of data points contributing to
%each measurement
meanTmp = particleBehavior.mean;
stdTmp  = particleBehavior.std;
numTmp  = particleBehavior.numPoints;
meanTmp(numTmp<25) = NaN;
stdTmp(numTmp<25) = NaN;

%get number of modes
numMode = size(meanTmp,3);

%get number of frames
minInc = -3;
numFrames = size(meanTmp,2) + (minInc - 1);

%open figure
hFig = figure('Name',figureName);
hold on

%plot for each mode
for iMode = 1 : numMode
    
    %protrusion
    %     subplot(numMode,3,(iMode-1)*3+1), hold on
    subplot(numMode,1,iMode), hold on
    plot(minInc:numFrames,meanTmp(1,:,iMode),'b','Marker','.')
    myErrorbar(minInc:numFrames,meanTmp(1,:,iMode),stdTmp(1,:,iMode)./sqrt(numTmp(1,:,iMode)))
    legend({'pr > re or pa'})
    %     plot(minInc:numFrames,meanTmp(2,:,iMode),'c','Marker','.')
    %     myErrorbar(minInc:numFrames,meanTmp(2,:,iMode),stdTmp(2,:,iMode)./sqrt(numTmp(2,:,iMode)))
    %     plot(minInc:numFrames,meanTmp(3,:,iMode),'r','Marker','.')
    %     myErrorbar(minInc:numFrames,meanTmp(3,:,iMode),stdTmp(3,:,iMode)./sqrt(numTmp(3,:,iMode)))
    %     plot(minInc:numFrames,meanTmp(4,:,iMode),'r--','Marker','.')
    %     myErrorbar(minInc:numFrames,meanTmp(4,:,iMode),stdTmp(4,:,iMode)./sqrt(numTmp(4,:,iMode)))
    %     plot(minInc:numFrames,meanTmp(5,:,iMode),'r:','Marker','.')
    %     myErrorbar(minInc:numFrames,meanTmp(5,:,iMode),stdTmp(5,:,iMode)./sqrt(numTmp(5,:,iMode)))
    %     plot(minInc:numFrames,meanTmp(6,:,iMode),'g','Marker','.')
    %     myErrorbar(minInc:numFrames,meanTmp(6,:,iMode),stdTmp(6,:,iMode)./sqrt(numTmp(6,:,iMode)))
    %     plot(minInc:numFrames,meanTmp(7,:,iMode),'g--','Marker','.')
    %     myErrorbar(minInc:numFrames,meanTmp(7,:,iMode),stdTmp(7,:,iMode)./sqrt(numTmp(7,:,iMode)))
    %     plot(minInc:numFrames,meanTmp(8,:,iMode),'g:','Marker','.')
    %     myErrorbar(minInc:numFrames,meanTmp(8,:,iMode),stdTmp(8,:,iMode)./sqrt(numTmp(8,:,iMode)))
    %     plot(minInc:numFrames,meanTmp(9,:,iMode),'k','Marker','.','LineWidth',2,'MarkerSize',10)
    %     myErrorbar(minInc:numFrames,meanTmp(9,:,iMode),stdTmp(9,:,iMode)./sqrt(numTmp(9,:,iMode)))
    %     legend({'pr>re','pr>un','pr>pa(s)>pr','pr>pa(s)>re','pr>pa(s)>un',...
    %         'pr>pa(l)>pr','pr>pa(l)>re','pr>pa(l)>un','all'})
    
    %     %retraction
    %     subplot(numMode,3,(iMode-1)*3+2), hold on
    %     plot(minInc:numFrames,meanTmp(10,:,iMode),'b','Marker','.')
    %     myErrorbar(minInc:numFrames,meanTmp(10,:,iMode),stdTmp(10,:,iMode)./sqrt(numTmp(10,:,iMode)))
    %     legend({'re > pr or pa'})
    %     %     plot(minInc:numFrames,meanTmp(11,:,iMode),'c','Marker','.')
    %     %     myErrorbar(minInc:numFrames,meanTmp(11,:,iMode),stdTmp(11,:,iMode)./sqrt(numTmp(11,:,iMode)))
    %     %     plot(minInc:numFrames,meanTmp(12,:,iMode),'r','Marker','.')
    %     %     myErrorbar(minInc:numFrames,meanTmp(12,:,iMode),stdTmp(12,:,iMode)./sqrt(numTmp(12,:,iMode)))
    %     %     plot(minInc:numFrames,meanTmp(13,:,iMode),'r--','Marker','.')
    %     %     myErrorbar(minInc:numFrames,meanTmp(13,:,iMode),stdTmp(13,:,iMode)./sqrt(numTmp(13,:,iMode)))
    %     %     plot(minInc:numFrames,meanTmp(14,:,iMode),'r:','Marker','.')
    %     %     myErrorbar(minInc:numFrames,meanTmp(14,:,iMode),stdTmp(14,:,iMode)./sqrt(numTmp(14,:,iMode)))
    %     %     plot(minInc:numFrames,meanTmp(15,:,iMode),'g','Marker','.')
    %     %     myErrorbar(minInc:numFrames,meanTmp(15,:,iMode),stdTmp(15,:,iMode)./sqrt(numTmp(15,:,iMode)))
    %     %     plot(minInc:numFrames,meanTmp(16,:,iMode),'g--','Marker','.')
    %     %     myErrorbar(minInc:numFrames,meanTmp(16,:,iMode),stdTmp(16,:,iMode)./sqrt(numTmp(16,:,iMode)))
    %     %     plot(minInc:numFrames,meanTmp(17,:,iMode),'g:','Marker','.')
    %     %     myErrorbar(minInc:numFrames,meanTmp(17,:,iMode),stdTmp(17,:,iMode)./sqrt(numTmp(17,:,iMode)))
    %     %     plot(minInc:numFrames,meanTmp(18,:,iMode),'k','Marker','.','LineWidth',2,'MarkerSize',10)
    %     %     myErrorbar(minInc:numFrames,meanTmp(18,:,iMode),stdTmp(18,:,iMode)./sqrt(numTmp(18,:,iMode)))
    %     %     legend({'re>pr','re>un','re>pa(s)>pr','re>pa(s)>re','re>pa(s)>un',...
    %     %         're>pa(l)>pr','re>pa(l)>re','re>pa(l)>un','all'})
    %
    %     %pause
    %     subplot(numMode,3,(iMode-1)*3+3), hold on
    %     %     plot(minInc:numFrames,meanTmp(19,:,iMode),'b','Marker','.')
    %     %     myErrorbar(minInc:numFrames,meanTmp(19,:,iMode),stdTmp(19,:,iMode)./sqrt(numTmp(19,:,iMode)))
    %     plot(minInc:numFrames,meanTmp(20,:,iMode),'b','Marker','.')
    %     myErrorbar(minInc:numFrames,meanTmp(20,:,iMode),stdTmp(20,:,iMode)./sqrt(numTmp(20,:,iMode)))
    %     plot(minInc:numFrames,meanTmp(21,:,iMode),'b:','Marker','.')
    %     myErrorbar(minInc:numFrames,meanTmp(21,:,iMode),stdTmp(21,:,iMode)./sqrt(numTmp(21,:,iMode)))
    %     plot(minInc:numFrames,meanTmp(22,:,iMode),'r','Marker','.')
    %     myErrorbar(minInc:numFrames,meanTmp(22,:,iMode),stdTmp(22,:,iMode)./sqrt(numTmp(22,:,iMode)))
    %     %     plot(minInc:numFrames,meanTmp(23,:,iMode),'r--','Marker','.')
    %     %     myErrorbar(minInc:numFrames,meanTmp(23,:,iMode),stdTmp(23,:,iMode)./sqrt(numTmp(23,:,iMode)))
    %     plot(minInc:numFrames,meanTmp(24,:,iMode),'r:','Marker','.')
    %     myErrorbar(minInc:numFrames,meanTmp(24,:,iMode),stdTmp(24,:,iMode)./sqrt(numTmp(24,:,iMode)))
    %     plot(minInc:numFrames,meanTmp(25,:,iMode),'g','Marker','.')
    %     myErrorbar(minInc:numFrames,meanTmp(25,:,iMode),stdTmp(25,:,iMode)./sqrt(numTmp(25,:,iMode)))
    %     plot(minInc:numFrames,meanTmp(26,:,iMode),'g--','Marker','.')
    %     myErrorbar(minInc:numFrames,meanTmp(26,:,iMode),stdTmp(26,:,iMode)./sqrt(numTmp(26,:,iMode)))
    %     plot(minInc:numFrames,meanTmp(27,:,iMode),'g:','Marker','.')
    %     myErrorbar(minInc:numFrames,meanTmp(27,:,iMode),stdTmp(27,:,iMode)./sqrt(numTmp(27,:,iMode)))
    %     %     legend({'pa>pr&<pr','pa>pr&<re','pa>pr&<un','pa>re&<pr','pa>re&<re',...
    %     %         'pa>re&<un','pa>un&<pr','pa>un&<re','pa>un&<un'})
    %     legend({'pa>pr&<re','pa>pr&<un','pa>re&<pr','pa>re&<un','pa>un&<pr','pa>un&<re','pa>un&<un'})
    
end

if ~isempty(saveLoc)
    saveas(hFig,saveLoc)
end


