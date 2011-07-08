function [autoCorrProtrusion,crossCorrSptProp] = ...
    particleBehaviorVsEdgeMotion(protPerWindow,sptPropInWindow,...
    prop2analyze,band2analyze,window2analyze,frameRange) %#ok<STOUT,INUSL>
%PARTICLEBEHAVIORVSEDGEMOTION looks for correlation between particle behavior and cell edge protrusion activity
%
%SYNOPSIS
%
%INPUT  protPerWindow  : Average protrusion vector per window (from
%                        Hunter).
%       sptPropInWindow: Output of particleBehaviorInWindows.
%       prop2analyze   : Cell array indicating properties to analyze. The
%                        names must exactly match the field names in
%                        sptPropInWindow.
%                        Optional. Default: All properties.
%       band2analyze   : Vector indicating which bands to analyze. Bands
%                        arise from the division of a cell into windows and
%                        run parallel to the cell edge. Band #1 is at the
%                        cell edge.
%                        Optional. Default: 1.
%      window2analyze  : Vector indicating which windows to analyze. 
%                        Optional. Default: all windows.
%      frameRange      : Row vector with first and last frame to include
%                        in analysis.
%                        Optional. Default: all frames.
%
%OUTPUT
%
%
%Khuloud Jaqaman, September 2010

%% Input

if nargin < 3 || isempty(prop2analyze)
    prop2analyze = {'spDensity','f2fDispMag2D','angleMean','angleStd',...
        'f2fDispMagParaDir','f2fDispMagPerpDir','f2fDispSignParaDir',...
        'f2fDispSignPerpDir','f2fDispMagParaProt','f2fDispMagPerpProt',...
        'f2fDispSignParaProt','f2fDispSignPerpProt','ratioDispMagDir',...
        'ratioDispSignDir','ratioDispMagProt','ratioDispSignProt',...
        'asymParam','f2fDispMag2DLin','f2fDispMagParaDirLin',...
        'f2fDispMarPerpDirLin','f2fDispSignParaDirLin',...
        'f2fDispSignPerpDirLin','f2fDispMagParaProtLin',...
        'f2fDispMagPerpProtLin','f2fDispSignParaProtLin',...
        'f2fDispSignPerpProtLin','ratioDispMagDirLin',...
        'ratioDispSignDirLin','ratioDispMagProtLin',...
        'ratioDispSignProtLin','asymParamLin',...
        'fracUnclass','fracLin','fracIso','fracIsoUnclass',...
        'fracConf','fracBrown','fracDir','diffCoef','confRad'};
end

%get bands to analyze
if nargin < 4 || isempty(band2analyze)
    band2analyze = 1;
end

%calculate number of properties and number of bands
numProp2analyze = length(prop2analyze);
numBands = length(band2analyze);

%get protrusion vector matrix
protNormVecMag = protPerWindow.avgNormal(:,1:end-1);

%get number of windows and time points
[numWindows,numTP] = size(protNormVecMag);

%get which windows to analyze and which to ignore
if nargin < 5 || isempty(window2analyze)
    window2analyze = 1 : numWindows;
end
window2ignore = setdiff(1:numWindows,window2analyze);

%get which frames to include in analysis
if nargin < 6 || isempty(frameRange)
    frameRange = [1 numTP];
else
    frameRange(1) = max(frameRange(1),1);
    frameRange(2) = min(frameRange(2),numTP);
end
numFramesAnalysis = diff(frameRange) + 1;

%determine maximum lag for correlation analysis
maxLag = floor(numFramesAnalysis/4);

%define extended color map
cmapExt = [[0 0 0]; colormap];
close

%% Analysis

%Protrusion auto-correlation ...

%initialize variables
protNormAllWindows = repmat(struct('observations',[]),numWindows,1);
autoCorrProtInd = NaN(maxLag+1,2,numWindows);

%set protrusion vector to NaN in ignored windows
protNormVecMag(window2ignore,:) = NaN;

%calculate auto-correlation for each individual window
for iWindow = window2analyze
    protNormAllWindows(iWindow).observations = protNormVecMag(iWindow,frameRange(1):frameRange(2))';
    autoCorrProtInd(:,:,iWindow) = autoCorr(protNormAllWindows(iWindow),maxLag);
end

%calculate for all windows combined
autoCorrProtComb = autoCorr(protNormAllWindows(window2analyze),maxLag,-1,1);

%estimate number of points used in auto-correlation calculation
tmp = vertcat(protNormAllWindows(window2analyze).observations);
numPointsAuto = length(find(~isnan(tmp)));

%store in output variable
autoCorrProtrusion = struct('timeLag',(0:maxLag)',...
    'indWindows',autoCorrProtInd,...
    'combined',autoCorrProtComb);

%Particle property cross-correlations and plotting ...

%go over requested particle properties
for iProp = 1 : numProp2analyze
    
    %get particle property name
    propName = prop2analyze{iProp};
    
    %get current particle property
    eval(['propCurrent = sptPropInWindow.' propName '.values;']);
    
    %initialize variables
    crossCorrBandInd = NaN(2*maxLag+1,2,numWindows,numBands);
    crossCorrBandComb = NaN(2*maxLag+1,2,numBands);
    
    %go over requested bands
    for iBand = 1 : numBands
        
        %get current particle property in current band for all windows and all times
        propCurrentBand = squeeze(propCurrent(band2analyze(iBand),:,:)); %#ok<NODEF>
        
        %set property value to NaN in ignored windows
        propCurrentBand(window2ignore,:) = NaN;
        
        %initialize one more variable
        partPropAllWindows = repmat(struct('observations',[]),numWindows,1);
        
        %go over all windows and cross-correlate particle property with
        %protrusion vector
        for iWindow = window2analyze
            partPropAllWindows(iWindow).observations = propCurrentBand(iWindow,frameRange(1):frameRange(2))';
            crossCorrBandInd(:,:,iWindow,iBand) = crossCorr(protNormAllWindows(iWindow),...
                partPropAllWindows(iWindow),maxLag);
        end
        
        %calculate cross-correlation for all windows combined
        crossCorrBandComb(:,:,iBand) = crossCorr(protNormAllWindows(window2analyze),partPropAllWindows(window2analyze),maxLag,1);
        
        %estimate number of points used in cross-correlation calculation
        tmp = vertcat(partPropAllWindows(window2analyze).observations);
        numPointsCross = length(find(~isnan(tmp)));
        
        %Plotting ...
        
        %open new figure
        figure('units','normalized','position',[0 0 1 1])
        
        %plot protrusion vectors
        subplot(2,3,1);
        protNormVecMagRange = protNormVecMag(:,frameRange(1):frameRange(2));
        imagesc(protNormVecMag);
        maxValue = max(abs(protNormVecMagRange(:)));
        nanValue = -maxValue - 2*maxValue/64;
        caxis([nanValue maxValue])
        colormap(cmapExt);
        colorbar
        title('Protrusion Vectors')
        xlim([frameRange(1)-0.5 frameRange(2)+0.5])
        xlabel('Time (frames)')
        ylabel('Position along edge (window #)')
        
        %plot auto-correlation of protrusion vectors in each window
        subplot(2,3,2);
        autoCorr2plot = squeeze(autoCorrProtInd(:,1,:))';
        nanValue = -1 - 2/64;
        autoCorr2plot(isnan(autoCorr2plot)) = nanValue;
        imagesc((0:maxLag),(1:numWindows),autoCorr2plot);
        caxis([nanValue 1])
        colormap(cmapExt);
        colorbar
        title('Auto-correlation of protrusion vector - Individual windows')
        xlabel('Time lag (frames)')
        ylabel('Position along edge (window #)')
        
        %plot auto-correlation of protrusion vectors combined over all
        %windows
        subplot(2,3,3)
        plot((0:maxLag),autoCorrProtComb(:,1));
        myErrorbar((0:maxLag),autoCorrProtComb(:,1),autoCorrProtComb(:,2));
        hold on
        plot([-maxLag-0.1 maxLag+0.1],1.96/sqrt(numPointsAuto)*[1 1],'k')
        plot([-maxLag-0.1 maxLag+0.1],-1.96/sqrt(numPointsAuto)*[1 1],'k')
        xlim([-maxLag-0.1 maxLag+0.1])
        ylim([-1.2 1.2])
        title('Auto-correlation of protrusion vector - All windows combined')
        xlabel('Time lag (frames)')
        ylabel('Auto-correlation')
        
        %plot current particle property
        subplot(2,3,4);
        propCurrentBandRange = propCurrentBand(:,frameRange(1):frameRange(2));
        switch propName
            case {'spDensity',...
                    'f2fDispMag2D','f2fDispMagParaDir',...
                    'f2fDispMarPerpDir','f2fDispSignParaDir',...
                    'f2fDispSignPerpDir','f2fDispMagParaProt',...
                    'f2fDispMagPerpProt','f2fDispSignParaProt',...
                    'f2fDispSignPerpProt','ratioDispMagDir',...
                    'ratioDispSignDir','ratioDispMagProt',...
                    'ratioDispSignProt','asymParam',...
                    'f2fDispMag2DLin','f2fDispMagParaDirLin',...
                    'f2fDispMarPerpDirLin','f2fDispSignParaDirLin',...
                    'f2fDispSignPerpDirLin','f2fDispMagParaProtLin',...
                    'f2fDispMagPerpProtLin','f2fDispSignParaProtLin',...
                    'f2fDispSignPerpProtLin','ratioDispMagDirLin',...
                    'ratioDispSignDirLin','ratioDispMagProtLin',...
                    'ratioDispSignProtLin','asymParamLin',...
                    'diffCoef','confRad'}
                valuesNoNaNs = propCurrentBandRange(:);
                valuesNoNaNs = valuesNoNaNs(~isnan(valuesNoNaNs));
                minValue = prctile(valuesNoNaNs,1);
                maxValue = prctile(valuesNoNaNs,99);
                nanValue = minValue - (maxValue - minValue)/64;
            case {'angleMean','angleStd','angleMeanLin','angleStdLin'}
                nanValue = -180/64;
                maxValue = 180;
            case {'fracUnclass','fracLin','fracIso','fracIsoUnclass',...
                    'fracConf','fracBrown','fracDir'}
                nanValue = -1/64;
                maxValue = 1;
        end
        prop2plot = propCurrentBand;
        prop2plot(isnan(prop2plot)) = nanValue;
        imagesc(prop2plot);
        caxis([nanValue maxValue])
        colormap(cmapExt);
        colorbar
        title([propName ' in Band ' num2str(band2analyze(iBand))])
        xlim([frameRange(1)-0.5 frameRange(2)+0.5])
        xlabel('Time (frames)')
        ylabel('Position along edge (window #)')
        
        %plot cross-correlation of particle property with protrusion vector
        %in each window
        subplot(2,3,5);
        crossCorr2plot = squeeze(crossCorrBandInd(:,1,:,iBand))';
        nanValue = -1 - 2/64;
        crossCorr2plot(isnan(crossCorr2plot)) = nanValue;
        imagesc((-maxLag:maxLag),(1:numWindows),crossCorr2plot);
        caxis([nanValue 1])
        colormap(cmapExt);
        colorbar
        title(['Cross-correlation of ' propName ' with protrusion vector in Band ' ...
            num2str(band2analyze(iBand)) ' - Individual windows'])
        xlabel('Time lag (frames)')
        ylabel('Position along edge (window #)')
        
        %plot cross-correlation combined over all windows
        subplot(2,3,6)
        plot((-maxLag:maxLag),crossCorrBandComb(:,1,iBand));
        myErrorbar((-maxLag:maxLag),crossCorrBandComb(:,1,iBand),crossCorrBandComb(:,2,iBand));
        hold on
        plot([-maxLag-0.1 maxLag+0.1],1.96/sqrt(numPointsCross)*[1 1],'k')
        plot([-maxLag-0.1 maxLag+0.1],-1.96/sqrt(numPointsCross)*[1 1],'k')
        xlim([-maxLag-0.1 maxLag+0.1])
        ylim([-1.2 1.2])
        title(['Cross-correlation of ' propName ' with protrusion vector in Band ' ...
            num2str(band2analyze(iBand)) ' - All windows combined'])
        xlabel('Time lag (frames)')
        ylabel('Cross-correlation')
        
    end %(for iBand = 1 : numBands)
    
    %store output variable
    propertyCrossCorr = struct('timeLag',(-maxLag:maxLag)',...
        'indWindows',crossCorrBandInd,...
        'combined',crossCorrBandComb); %#ok<NASGU>
    eval(['crossCorrSptProp.' propName ' = propertyCrossCorr;']);
    
end %(for iProp = 1 : numProp2analyze)




%% Scatter plots of individual properties vs. protrusion vectors

%NOT TOO USEFUL

% %go over requested properties
% for iProp = 1 : numProp2analyze
%
%     %get property name
%     propName = prop2analyze{iProp};
%
%     %get current property
%     eval(['propCurrent = sptPropInWindow.' propName '.values;']);
%
%     %go over all bands
%     for iBand = 1 : numBands
%         propCurrentBand = squeeze(propCurrent(band2analyze(iBand),:,:));
%         propCurrentBand(isnan(propCurrentBand)) = -1;
%         figure
%         plot(protNormVecMag',propCurrentBand','.')
%         xlabel('Protrusion vector')
%         ylabel(propName)
%     end
%
% end

