function [] = makiOverlayPlots(analStrucArray,labelArray,plotConfidence,...
    colorArray,maxLag,timeInterval)

nStrucs = length(analStrucArray);

if nargin <2
    labelArray = {};
    for i = 1:nStrucs
        labelArray = [ labelArray analStrucArray(i).fileName];
    end
end

if nargin < 3 || isempty(plotConfidence)
    plotConfidence = 0;
end

if nargin < 4 || isempty(colorArray)
    if nStrucs <= 8
        colorArray = ['k' 'r' 'b' 'c' 'm' 'g' 'y' 'w']';
    else
        colorArray = rand(nStrucs,3);
    end
end

if nargin < 5 || isempty(maxLag)
    maxLag = 20;
end

if nargin < 6 || isempty(timeInterval)
    timeInterval = 7.5;
end



%% Separation change autocorrelation, displacement projection cross-correlation and center displacement autocorrelation

figure

for i = 1:nStrucs

    %separation change autocorrelation
    subplot(1,3,1)
    hold on
    plot((0:maxLag)*timeInterval,analStrucArray(i).sisterConnection.Inlier.autocorr.all.rateChangeDist(:,1),'Color',colorArray(i,:));


    if i == nStrucs

        xlabel('Time Lag (s)')
        ylabel('Sister separation change autocorrelation')

        if ~isempty(labelArray)
            legend(labelArray)
        end
    end

    %displacement projection crosscorrelation
    subplot(1,3,2)
    hold on
    plot((-maxLag:maxLag)*timeInterval,analStrucArray(i).sisterMotionCoupling.Inlier.crosscorr.all.projections(:,1),'Color',colorArray(i,:));

    if i == nStrucs

        xlabel('Time Lag (s)')
        ylabel('Sister displacement projection crosscorrelation')

        if ~isempty(labelArray)
            legend(labelArray)
        end
    end
    
    %center displacement along normal autocorrelation
    subplot(1,3,3)
    hold on
    plot((0:maxLag)*timeInterval,analStrucArray(i).sepDispSpaceTime.Inlier.autocorr.all.centerPosChange(:,1),'Color',colorArray(i,:));

    if i == nStrucs

        xlabel('Time Lag (s)')
        ylabel('Center normal displacement autocorrelation')

        if ~isempty(labelArray)
            legend(labelArray)
        end
    end
    
end

%Plot the 95% confidence intervals afterwards to prevent legend confusion
if plotConfidence

    for i = 1:nStrucs

        %separation change autocorrelation
        subplot(1,3,1)
        hold on
        plot((0:maxLag)*timeInterval,analStrucArray(i).sisterConnection.Inlier.autocorr.all.rateChangeDist(:,1) + ...
            analStrucArray(i).sisterConnection.Inlier.autocorr.all.rateChangeDist(:,2)*2,'--','Color',colorArray(i,:))

        plot((0:maxLag)*timeInterval,analStrucArray(i).sisterConnection.Inlier.autocorr.all.rateChangeDist(:,1) - ...
            analStrucArray(i).sisterConnection.Inlier.autocorr.all.rateChangeDist(:,2)*2,'--','Color',colorArray(i,:))


        %displacement projection crosscorrelation
        subplot(1,3,2)
        hold on
        plot((-maxLag:maxLag)*timeInterval,analStrucArray(i).sisterMotionCoupling.Inlier.crosscorr.all.projections(:,1) + ...
            analStrucArray(i).sisterMotionCoupling.Inlier.crosscorr.all.projections(:,2)*2,'--','Color',colorArray(i,:));

        plot((-maxLag:maxLag)*timeInterval,analStrucArray(i).sisterMotionCoupling.Inlier.crosscorr.all.projections(:,1) - ...
            analStrucArray(i).sisterMotionCoupling.Inlier.crosscorr.all.projections(:,2)*2,'--','Color',colorArray(i,:));

        %center displacement along normal autocorrelation
        subplot(1,3,3)
        hold on
        plot((0:maxLag)*timeInterval,analStrucArray(i).sepDispSpaceTime.Inlier.autocorr.all.centerPosChange(:,1) + ...
            analStrucArray(i).sepDispSpaceTime.Inlier.autocorr.all.centerPosChange(:,2)*2,'--','Color',colorArray(i,:))

        plot((0:maxLag)*timeInterval,analStrucArray(i).sepDispSpaceTime.Inlier.autocorr.all.centerPosChange(:,1) - ...
            analStrucArray(i).sepDispSpaceTime.Inlier.autocorr.all.centerPosChange(:,2)*2,'--','Color',colorArray(i,:))

    end

end %if plotconfidence

%% Histograms for sister separation and displacement along normal

% %get the distributions
% for iStruc = 1 : nStrucs   
%     sisterSep(iStruc).observations = analStrucArray(iStruc).sepDispSpaceTime.Inlier.distribution.sisterSeparation;
%     sisterSepChange(iStruc).observations = analStrucArray(iStruc).sepDispSpaceTime.Inlier.distribution.sisterSepChange;
%     sepPChangeInt(iStruc).observations = analStrucArray(iStruc).sepDispSpaceTime.Inlier.distribution.sepPChangeInterval;
%     sepNChangeInt(iStruc).observations = analStrucArray(iStruc).sepDispSpaceTime.Inlier.distribution.sepNChangeInterval;
%     centerPos(iStruc).observations = analStrucArray(iStruc).sepDispSpaceTime.Inlier.distribution.centerPosition;
%     centerPosChange(iStruc).observations = analStrucArray(iStruc).sepDispSpaceTime.Inlier.distribution.centerPosChange;
%     centerPosPChangeInt(iStruc).observations = analStrucArray(iStruc).sepDispSpaceTime.Inlier.distribution.centerPosPChangeInterval;
%     centerPosNChangeInt(iStruc).observations = analStrucArray(iStruc).sepDispSpaceTime.Inlier.distribution.centerPosNChangeInterval;
% end
% 
% %get minima and maxima to to calculate distribution upper and lower bounds
% tmp = vertcat(sisterSep.observations);
% sisterSepMin = min(tmp);
% sisterSepMax = min(tmp);
% tmp = vertcat(sisterSepChange.observations);
% sisterSepChangeMin = min(tmp);
% sisterSepChangeMax = max(tmp);
% tmp = [vertcat(sepPChangeInt.observations); vertcat(sepNChangeInt.observations)];
% sepPNChangeIntMin = min(tmp);
% sepPNChangeIntMax = max(tmp);
% 
% tmp = vertcat(centerPos.observations);
% centerPosMin = min(tmp);
% centerPosMax = min(tmp);
% tmp = vertcat(centerPosChange.observations);
% centerPosChangeMin = min(tmp);
% centerPosChangeMax = max(tmp);
% tmp = [vertcat(centerPosPChangeInt.observations); vertcat(centerPosNChangeInt.observations)];
% centerPosPNChangeIntMin = min(tmp);
% centerPosPNChangeIntMax = max(tmp);


