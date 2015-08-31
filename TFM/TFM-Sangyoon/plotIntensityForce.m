function []=plotIntensityForce(tracksNA, fileStore,alignEvent,indivColor,varargin)
% plotIntensityForce(tracksNA) plots intensities and forces of tracks
% w.r.t. time frames and store in designated folder
% idGroup1f = find(idGroup1);
% Sangyoon Han, May 2015

ip = inputParser;
ip.addOptional('fileStore',[],@(x) isempty(x) || ischar(x));
ip.addOptional('alignEvent',false,@islogical);
ip.addOptional('indivColor',false,@islogical);
ip.addParamValue('Source','ampTotal',@ischar);

ip.parse(fileStore,alignEvent,indivColor,varargin{:});

alignEvent=ip.Results.alignEvent;
indivColor=ip.Results.indivColor;
source = ip.Results.Source;

h=figure; hold on
nTracks = numel(tracksNA);
% alignment might be necessary:
% events = detectProtrusionEvents(v,dThreshold)
lifeTime = arrayfun(@(x) x.lifeTime,tracksNA);

frameMaxAmp = zeros(nTracks,1);
if alignEvent
    % Find the maxima using some of the core function in detectProtrusionEvents
    for ii=1:nTracks
        d = tracksNA(ii).ampTotal;
        splineParam=0.01;

        %   perform spline filter for ampTotal
        nTime = length(d);
        sd_spline= csaps(1:nTime,d,splineParam);
        sd=ppval(sd_spline,1:nTime);
        % Find the maximum
        [~,curFrameMaxAmp]=max(sd);
        frameMaxAmp(ii)=curFrameMaxAmp;
    end
    % Find the mean time point
    meanFrameDouble = mean(frameMaxAmp);
    meanFrame = floor(meanFrameDouble);
    % Shift each time series w.r.t. the mean time point(meanFrame)
    framesToShift = frameMaxAmp - meanFrame;
    lifeAfterMaxAmp = lifeTime - frameMaxAmp;
    if length(lifeTime)>30
        thresLifeAfterMaxAmp = floor(quantile(lifeAfterMaxAmp,0.8));
    elseif length(lifeTime)>20
        thresLifeAfterMaxAmp = floor(quantile(lifeAfterMaxAmp,0.7));
    else
        thresLifeAfterMaxAmp = floor(quantile(lifeAfterMaxAmp,0.6));
    end
    nSampleFrames = (meanFrame+thresLifeAfterMaxAmp);
    AmpArray = NaN(nTracks,nSampleFrames);
    forceArray = NaN(nTracks,nSampleFrames);
    if strcmp(source,'edgeAdvanceDist')
        edgeDistArray = NaN(nTracks,nSampleFrames);
    end
    for ii=1:nTracks
        curAmp = tracksNA(ii).ampTotal(logical(tracksNA(ii).presence));
        curForce = tracksNA(ii).forceMag(logical(tracksNA(ii).presence));
        
        curFrameRange = tracksNA(ii).iFrame(logical(tracksNA(ii).presence));
        curFrameRangeShifted = curFrameRange - framesToShift(ii);
        curFrameRangeShifted = curFrameRangeShifted(curFrameRangeShifted>0 & curFrameRangeShifted<=nSampleFrames & curFrameRangeShifted<=nSampleFrames);
        AmpArray(ii,curFrameRangeShifted) = curAmp(curFrameRange>framesToShift(ii) & (curFrameRange - framesToShift(ii))<= nSampleFrames);
        forceArray(ii,curFrameRangeShifted) = curForce(curFrameRange>framesToShift(ii) & (curFrameRange - framesToShift(ii))<= nSampleFrames);
        if strcmp(source,'edgeAdvanceDist')
            curEdgeDist = tracksNA(ii).edgeAdvanceDist(logical(tracksNA(ii).presence));
            edgeDistArray(ii,curFrameRangeShifted) = curEdgeDist(curFrameRange>framesToShift(ii) & (curFrameRange - framesToShift(ii))<= nSampleFrames);
        end
    end
    if indivColor
        if strcmp(source,'edgeAdvanceDist')
            subplot(1,3,1), plot(1:nSampleFrames,AmpArray), hold on%,'Color',[0.5 0.5 0.5]), hold on
            subplot(1,3,2), plot(1:nSampleFrames,forceArray), hold on%,'Color',[240/255 128/255 128/255]), hold on
            subplot(1,3,3), plot(1:nSampleFrames,edgeDistArray), hold on%,'Color',[240/255 128/255 128/255]), hold on
        else
            subplot(1,2,1), plot(1:nSampleFrames,AmpArray), hold on%,'Color',[0.5 0.5 0.5]), hold on
            subplot(1,2,2), plot(1:nSampleFrames,forceArray), hold on%,'Color',[240/255 128/255 128/255]), hold on
        end
    else
        if strcmp(source,'edgeAdvanceDist')
            subplot(1,3,1), plot(1:nSampleFrames,AmpArray, 'Color',[0.5 0.5 0.5]), hold on
            subplot(1,3,2), plot(1:nSampleFrames,forceArray, 'Color',[240/255 128/255 128/255]), hold on
            subplot(1,3,3), plot(1:nSampleFrames,edgeDistArray, 'Color',[10/255 220/255 64/255]), hold on
        else
            subplot(1,2,1), plot(1:nSampleFrames,AmpArray, 'Color',[0.5 0.5 0.5]), hold on
            subplot(1,2,2), plot(1:nSampleFrames,forceArray, 'Color',[240/255 128/255 128/255]), hold on
        end
    end        
    maxLifeTime = nSampleFrames;
else    
    for ii=1:nTracks
        d = tracksNA(ii).ampTotal(logical(tracksNA(ii).presence));
        if indivColor
            if strcmp(source,'edgeAdvanceDist')
                subplot(1,3,1), plot(1:tracksNA(ii).lifeTime,d), hold on%,'Color',[0.5 0.5 0.5]), hold on
                subplot(1,3,2), plot(1:tracksNA(ii).lifeTime,tracksNA(ii).forceMag(logical(tracksNA(ii).presence))), hold on%,'Color',[240/255 128/255 128/255]), hold on
                subplot(1,3,3), plot(1:tracksNA(ii).lifeTime,tracksNA(ii).edgeAdvanceDist(logical(tracksNA(ii).presence))), hold on%,'Color',[240/255 128/255 128/255]), hold on
            else
                subplot(1,2,1), plot(1:tracksNA(ii).lifeTime,d), hold on%,'Color',[0.5 0.5 0.5]), hold on
                subplot(1,2,2), plot(1:tracksNA(ii).lifeTime,tracksNA(ii).forceMag(logical(tracksNA(ii).presence))), hold on%,'Color',[240/255 128/255 128/255]), hold on
            end
        else
            if strcmp(source,'edgeAdvanceDist')
                subplot(1,3,1), plot(1:tracksNA(ii).lifeTime,d,'Color',[0.5 0.5 0.5]), hold on
                subplot(1,3,2), plot(1:tracksNA(ii).lifeTime,tracksNA(ii).forceMag(logical(tracksNA(ii).presence)),'Color',[240/255 128/255 128/255]), hold on
                subplot(1,3,3), plot(1:tracksNA(ii).lifeTime,tracksNA(ii).edgeAdvanceDist(logical(tracksNA(ii).presence)),'Color',[10/255 220/255 64/255]), hold on
            else
                subplot(1,2,1), plot(1:tracksNA(ii).lifeTime,d,'Color',[0.5 0.5 0.5]), hold on
                subplot(1,2,2), plot(1:tracksNA(ii).lifeTime,tracksNA(ii).forceMag(logical(tracksNA(ii).presence)),'Color',[240/255 128/255 128/255]), hold on
            end
        end
    end
    % title(['Group 1:' num2str(idGroup1f')])
    % set the life time to be 80 percentile
    if length(lifeTime)>30
        thresLifeTime = quantile(lifeTime,0.8);
    elseif length(lifeTime)>20
        thresLifeTime = quantile(lifeTime,0.7);
    else
        thresLifeTime = quantile(lifeTime,0.6);
    end

    nSampleFrames = floor(thresLifeTime);
    AmpArray = NaN(nTracks,nSampleFrames);
    forceArray = NaN(nTracks,nSampleFrames);
    if strcmp(source,'edgeAdvanceDist')
        edgeDistArray = NaN(nTracks,nSampleFrames);
    end
    p=0;
    for ii=1:nTracks
        p=p+1;
        d = tracksNA(ii).ampTotal(logical(tracksNA(ii).presence));
        curAmp = d;
        fmax = min(nSampleFrames, length(curAmp));
        AmpArray(p,1:fmax) = curAmp(1:fmax);

        curForce = tracksNA(ii).forceMag(logical(tracksNA(ii).presence));
        fmax = min(nSampleFrames, length(curForce));
        forceArray(p,1:fmax) = curForce(1:fmax);

        if strcmp(source,'edgeAdvanceDist')
            curEdgeDist = tracksNA(ii).edgeAdvanceDist(logical(tracksNA(ii).presence));
            edgeDistArray(p,1:fmax) = curEdgeDist(1:fmax);
        end
    end
    maxLifeTime = quantile(lifeTime,0.99);
end    
if nTracks<5
    AmpG1avg = nanmean(AmpArray,1)';
    % AmpG1std = nanstd(AmpArrayG1,1)';
    forceG1avg = nanmean(forceArray,1)';
    % forceG1std = nanstd(forceArrayG1,1)';
    if strcmp(source,'edgeAdvanceDist')
        edgeDistAvg = nanmean(edgeDistArray,1)';
    end
else
    AmpG1avg = nanmedian(AmpArray,1)';
    % AmpG1std = nanstd(AmpArrayG1,1)';
    forceG1avg = nanmedian(forceArray,1)';
    % forceG1std = nanstd(forceArrayG1,1)';
    if strcmp(source,'edgeAdvanceDist')
        edgeDistAvg = nanmedian(edgeDistArray,1)';
    end
end
if strcmp(source,'edgeAdvanceDist')
    set(h,'Position',[200,300,600,200])%,title(['ID:' num2str(ii) ', CC-score:' num2str(tracksNA(ii).CCscore)])
    subplot(1,3,1)
else
    set(h,'Position',[200,400,400,200])%,title(['ID:' num2str(ii) ', CC-score:' num2str(tracksNA(ii).CCscore)])
    subplot(1,2,1) 
end
plot(1:nSampleFrames,AmpG1avg,'k','Linewidth',3)
xlim([0 maxLifeTime])
maxYamp = quantile(nanmax(AmpArray),0.99);
minYamp = quantile(nanmin(AmpArray),0.01);
ylim([minYamp maxYamp])
xlabel('Time (frame)')
ylabel('Fluorescence intensity (a.u.)')
if strcmp(source,'edgeAdvanceDist')
    subplot(1,3,2)
else
    subplot(1,2,2) 
end
plot(1:nSampleFrames,forceG1avg,'r','Linewidth',3)
xlim([0 maxLifeTime])
maxYforce = quantile(nanmax(forceArray),0.95);
minYforce = quantile(nanmin(forceArray),0.01);
ylim([minYforce maxYforce])
xlabel('Time (frame)')
ylabel('Traction (Pa)')

if strcmp(source,'edgeAdvanceDist')
    subplot(1,3,3)
    plot(1:nSampleFrames,edgeDistAvg,'Color',[30/255 160/255 24/255],'Linewidth',3)
    xlim([0 maxLifeTime])
    maxYED = quantile(nanmax(edgeDistAvg),0.95);
    minYED = quantile(nanmin(edgeDistAvg),0.01);
    ylim([minYED maxYED])
    xlabel('Time (frame)')
    ylabel('Edge distance (Pixel)')
end

[pathStore]=fileparts(fileStore);
if ~exist(pathStore,'dir')
    mkdir(pathStore)
end
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-r300', strcat(fileStore));

