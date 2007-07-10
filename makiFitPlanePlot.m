function makiFitPlanePlot(dataStruct)
%MAKIFITPLANEPLOT is the plotting function for makiFitPlane
%
% SYNOPSIS: makiFitPlanePlot(dataStruct)
%
% INPUT dataStruct: maki data structure (see makiMakeDataStruct for details)
%
% OUTPUT 
%
% REMARKS
%
% created with MATLAB ver.: 7.3.0.267 (R2006b) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 10-Jul-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input for empty
if nargin == 0 || isempty(dataStruct)
    dataStruct = makiLoadDataFile;
end
% check whether analysis has been done
if isempty(dataStruct.planeFit)
    dataStruct = makiFitPlane(dataStruct,0);
end

% read variables
nTimepoints = dataStruct.dataProperties.movieSize(end);
planeFit = dataStruct.planeFit;
eigenValues = cat(1,planeFit.eigenValues);

% plot (ratio of) eigenvalues, boundaries of metaphase
figure('Name',sprintf('EigenRatio and distribution for %s',dataStruct.projectName))
    time = (1:nTimepoints)'.*dataStruct.dataProperties.timeLapse/60; % time in minutes
    subplot(2,1,1)
    plot(time,10*eigenValues(:,1)./mean(eigenValues(:,2:3),2),'-ro',...
        time,eigenValues(:,1),'--m',...
        time,mean(eigenValues(:,2:3),2),'--b+',...
        time,eigenValues(:,2),'.b',...
        time,eigenValues(:,3),'.b');
    legend('10x eigenRatio','different eigenValue','meanSimilar','similar1','similar2')
    xMinMax = [0,(nTimepoints+1)*dataStruct.dataProperties.timeLapse/60];
    xlim(xMinMax)
    ylim([0 50])
    xlabel('Time (min)')
    hold on
    % plot lines showing cutoff and mitosis
    plot(xMinMax,[10 10],':k')
    plot(planeFit(1).metaphaseFrames([1 1])*dataStruct.dataProperties.timeLapse/60,[0 20],'k',...
        planeFit(1).metaphaseFrames([end end])*dataStruct.dataProperties.timeLapse/60,[0 20],'k')
    text(planeFit(1).metaphaseFrames(1)*dataStruct.dataProperties.timeLapse/60+1,13,'Mitosis')

    % distribution parameters
    distParms = catStruct(3,'planeFit.distParms');
    subplot(2,1,2)
    plot(repmat(time(planeFit(1).metaphaseFrames),1,3),squeeze(distParms(1:3,1,:))')
    hold on
    plot(time(planeFit(1).metaphaseFrames),10*squeeze(distParms(4,1,:))','k')
    xlim(xMinMax)
    legend('variance','skewness','kurtosis','10*p(normal)')

    figure('Name',sprintf('P-value and orientation changes for %s',dataStruct.projectName))
    time = (1.5:nTimepoints-0.5)'.*dataStruct.dataProperties.timeLapse/60; % time in minutes
    subplot(2,1,1)
    plot(time(planeFit(1).metaphaseFrames(1:end-1)),-log10(cat(1,planeFit.deltaP)))
    title('p-value for change in distribution (-log(p))')
    % plot 5% and 0.1% thresholds
    hold on
    plot(xMinMax,-log10([0.05 0.05]),':k')
    text(time(1),-log10(0.05)+0.15,'5% threshold')
    plot(xMinMax,-log10([0.001 0.001]),':k')
    text(time(1),-log10(0.001)+0.15,'0.1% threshold')
    xlim(xMinMax)
    ylim([0,-log10(0.001)+0.5])
    subplot(2,1,2)
    plot(time(planeFit(1).metaphaseFrames(1:end-1)),(cat(1,planeFit.deltaAngle)))
    title('orientation change (deg)')
    xlim(xMinMax)

    % histogram of spot distribution
    c=-5.75:0.5:5.75;
    nBins = length(c);
    z = zeros(nBins,nTimepoints);
    for t=planeFit(1).metaphaseFrames',
        z(:,t)=hist(planeFit(t).planeCoord(:,1),c);
    end
    uiViewPanel,
    imshow(z,[]),
    cm=isomorphicColormap('b');
    colormap(cm)
