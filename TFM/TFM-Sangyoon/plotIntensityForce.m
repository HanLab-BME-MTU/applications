function [longestCohorts,h]=plotIntensityForce(tracksNA, fileStore,alignEvent,indivColor,varargin)
% plotIntensityForce(tracksNA) plots intensities and forces of tracks
% w.r.t. time frames and store in designated folder
% idGroup1f = find(idGroup1);
% Sangyoon Han, May 2015

ip = inputParser;
ip.addOptional('fileStore',[],@(x) isempty(x) || ischar(x));
ip.addOptional('alignEvent',false,@islogical);
ip.addOptional('indivColor',false,@islogical);
ip.addParamValue('Source',{'ampTotal','forceMag'},@iscell);
ip.addParamValue('UseCurrentAxis',false,@islogical);
ip.addParamValue('plotCohorts',false,@islogical);
ip.addParamValue('tInterval',[],@isnumeric);
ip.addParamValue('prePostFrames',[],@isnumeric);
ip.addParamValue('yNormalization',false,@islogical);
ip.addParamValue('onlyFirstMode',false,@islogical);
ip.addParamValue('plotConfInt',false,@islogical);
ip.addParamValue('numCohorts',4,@isnumeric);

ip.parse(fileStore,alignEvent,indivColor,varargin{:});

alignEvent=ip.Results.alignEvent;
indivColor=ip.Results.indivColor;
source = ip.Results.Source;
useCurrentAxis = ip.Results.UseCurrentAxis;
plotCohorts = ip.Results.plotCohorts;
tInterval = ip.Results.tInterval;
prePostFrames = ip.Results.prePostFrames;
yNormalization = ip.Results.yNormalization;
onlyFirstMode = ip.Results.onlyFirstMode;
plotConfInt= ip.Results.plotConfInt;
numCohorts= ip.Results.numCohorts;

if isempty(tInterval)
    tInterval_used = 1;
else
    tInterval_used = tInterval;
end
if isempty(prePostFrames)
    prePostFramesUsed = 0;
else
    prePostFramesUsed = prePostFrames;
end

% make sure the lifetime calculation is okay
tracksNA = recalculateLifeTimeTracks(tracksNA);
nTracks = numel(tracksNA);
if ~useCurrentAxis
    h=figure; 
else
    h=[];
end
hold on
% alignment might be necessary:
% events = detectProtrusionEvents(v,dThreshold)
lifeTime = arrayfun(@(x) x.lifeTime,tracksNA);
nSources=numel(source);

frameMaxAmp = zeros(nTracks,1);
longestCohorts=[];

if nTracks<1
    disp('No data in tracksNA ... Aborting ...')
    return
end

if ~useCurrentAxis
    if strcmp(source,'edgeAdvanceDist')
        set(h,'Position',[200,300,600,200])%,title(['ID:' num2str(ii) ', CC-score:' num2str(tracksNA(ii).CCscore)])
        subplot(1,3,1)
    else
        set(h,'Position',[200,400,400,200])%,title(['ID:' num2str(ii) ', CC-score:' num2str(tracksNA(ii).CCscore)])
    end
end

if alignEvent
    % Find the maxima using some of the core function in detectProtrusionEvents
    for ii=1:nTracks
        d = tracksNA(ii).ampTotal;
        splineParam=0.01;

        %   perform spline filter for ampTotal
        nTime = length(d);
        tRange = 1:nTime;
        numNan = find(isnan(d),1,'last');
        tRange(isnan(d)) = [];
        d(isnan(d)) = [];
        sd_spline= csaps(tRange,d,splineParam);
        sd=ppval(sd_spline,tRange);
%         d = [NaN(1,numNan) d];
%         tRange = [NaN(1,numNan) tRange];
        sd = [NaN(1,numNan) sd];
%         sd(isnan(d)) = NaN;
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
    if length(lifeTime)>1000
        thresLifeAfterMaxAmp = floor(quantile(lifeAfterMaxAmp,0.99));
    elseif length(lifeTime)>100
        thresLifeAfterMaxAmp = floor(quantile(lifeAfterMaxAmp,0.95));
    elseif length(lifeTime)>30
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
        curAmp = (curAmp);%-min(curAmp));%/(max(curAmp)-min(curAmp));
        curForce = tracksNA(ii).forceMag(logical(tracksNA(ii).presence));
        curForce = (curForce);%-min(curForce));%/(max(curForce)-min(curForce));
        
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
        minYamp = nanmin(AmpArray(:));
        maxYamp  = nanmax(AmpArray(:));      
    else
        if strcmp(source,'edgeAdvanceDist')
            subplot(1,3,1), plot(1:nSampleFrames,AmpArray, 'Color',[0.5 0.5 0.5]), hold on
            subplot(1,3,2), plot(1:nSampleFrames,forceArray, 'Color',[240/255 128/255 128/255]), hold on
            subplot(1,3,3), plot(1:nSampleFrames,edgeDistArray, 'Color',[10/255 220/255 64/255]), hold on
        else
            if nSources>1
                subplot(1,2,1)
                if plotConfInt
                    plotTimeSeriesConfInt((1:nSampleFrames)*tInterval,AmpArray, 'Color',[0.5 0.5 0.5])
                else
                    plot((1:nSampleFrames)*tInterval,AmpArray, 'Color',[0.5 0.5 0.5]), hold on
                end
            end
            if strcmp(source{1},'ampTotal') && ~plotConfInt
                plot(1:nSampleFrames,AmpArray, 'Color',[0.5 0.5 0.5]), hold on
                minYamp=nanmin(AmpArray(:));
                maxYamp=nanmax(AmpArray(:));
            end
            if nSources>1
                subplot(1,2,2), 
                if plotConfInt
                    plotTimeSeriesConfInt((1:nSampleFrames)*tInterval,forceArray, 'Color', [240/255 128/255 128/255])
                else
                    plot((1:nSampleFrames)*tInterval,forceArray, 'Color',[240/255 128/255 128/255]), hold on
                end
            end
            if strcmp(source{1},'forceMag') && ~plotConfInt
                plot(1:nSampleFrames,forceArray, 'Color',[240/255 128/255 128/255]), hold on
                minYamp=nanmin(forceArray(:));
                maxYamp=nanmax(forceArray(:));
            end
        end
    end        
    maxLifeTime = nSampleFrames;
else  
    if ~plotCohorts
        for ii=1:nTracks
%             d = getfield(tracksNA(ii),{1},source{1},{find(tracksNA(ii).presence)});
            if indivColor
                d=tracksNA(ii).ampTotal(logical(tracksNA(ii).presence));
                shiftedTime = 1:sum(tracksNA(ii).presence);
                if strcmp(source,'edgeAdvanceDist')
                    subplot(1,3,1), plot(shiftedTime,d), hold on%,'Color',[0.5 0.5 0.5]), hold on
                    subplot(1,3,2), plot(shiftedTime,tracksNA(ii).forceMag(logical(tracksNA(ii).presence))), hold on%,'Color',[240/255 128/255 128/255]), hold on
                    subplot(1,3,3), plot(shiftedTime,tracksNA(ii).edgeAdvanceDist(logical(tracksNA(ii).presence))), hold on%,'Color',[240/255 128/255 128/255]), hold on
                else
                    subplot(1,2,1), plot(shiftedTime,d), hold on%,'Color',[0.5 0.5 0.5]), hold on
                    subplot(1,2,2), plot(shiftedTime,tracksNA(ii).forceMag(logical(tracksNA(ii).presence))), hold on%,'Color',[240/255 128/255 128/255]), hold on
                end
            else
                for kk=1:nSources
                    if ~useCurrentAxis
                        subplot(1,nSources,kk)
                        if ii==1
                            hold on
                        end
                    end
                    plot(1:tracksNA(ii).lifeTime+1,getfield(tracksNA(ii),{1},source{kk},{tracksNA(ii).startingFrameExtra:tracksNA(ii).endingFrameExtra}),'Color',[0.5 0.5 0.5])
                end
    %             if strcmp(source,'edgeAdvanceDist')
    %                 subplot(1,3,1), plot(1:tracksNA(ii).lifeTime,d,'Color',[0.5 0.5 0.5]), hold on
    %                 subplot(1,3,2), plot(1:tracksNA(ii).lifeTime,tracksNA(ii).forceMag(logical(tracksNA(ii).presence)),'Color',[240/255 128/255 128/255]), hold on
    %                 subplot(1,3,3), plot(1:tracksNA(ii).lifeTime,tracksNA(ii).edgeAdvanceDist(logical(tracksNA(ii).presence)),'Color',[10/255 220/255 64/255]), hold on
    %             else
    %                 subplot(1,2,1), plot(1:tracksNA(ii).lifeTime+1,d,'Color',[0.5 0.5 0.5]), hold on
    %                 subplot(1,2,2), plot(1:tracksNA(ii).lifeTime+1,tracksNA(ii).forceMag(logical(tracksNA(ii).presence)),'Color',[240/255 128/255 128/255]), hold on
    %             end
            end
        end
        if strcmp(source,'edgeAdvanceDist')
            subplot(1,3,1), ylabel('F.I. (a.u.)')
            subplot(1,3,2), ylabel('Traction (Pa)')
            subplot(1,3,3), ylabel('Edge Advance (Pixel)')
        else
            subplot(1,2,1), ylabel('F.I. (a.u.)')
            subplot(1,2,2), ylabel('Traction (Pa)')
        end
        
        % title(['Group 1:' num2str(idGroup1f')])
        % set the life time to be 80 percentile
        if length(lifeTime)>1000
            thresLifeTime = quantile(lifeTime,0.99);
        elseif length(lifeTime)>100
            thresLifeTime = quantile(lifeTime,0.95);
        elseif length(lifeTime)>50
            thresLifeTime = quantile(lifeTime,0.9);
        elseif length(lifeTime)>30
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
            d=getfield(tracksNA(ii),{1},source{1},{find(tracksNA(ii).presence)});
%             d = tracksNA(ii).ampTotal(logical(tracksNA(ii).presence));
            curAmp = d;
            fmax = min(nSampleFrames, length(curAmp));
            AmpArray(p,1:fmax) = curAmp(1:fmax);
            
            if numel(source)>1
                curForce = tracksNA(ii).forceMag(logical(tracksNA(ii).presence));
                fmax = min(nSampleFrames, length(curForce));
                forceArray(p,1:fmax) = curForce(1:fmax);
            end
            if strcmp(source,'edgeAdvanceDist')
                curEdgeDist = tracksNA(ii).edgeAdvanceDist(logical(tracksNA(ii).presence));
                edgeDistArray(p,1:fmax) = curEdgeDist(1:fmax);
            end
        end
        maxLifeTime = nSampleFrames;
        maxYamp = quantile(nanmax(AmpArray),0.99);
        minYamp = quantile(nanmin(AmpArray),0.01);
        if isempty(tInterval)
            xlabel('Time (frame)')
        else
            xlabel('Time (sec)')
        end
        maxLifeTime = nSampleFrames;
    else
        % plot cohorts
%         [N,edges,bin]=histcounts(lifeTime);
        % distribution was not normal. Binning depeding on quantile..
%         numCohorts = 2;
        stepPrc=100/(numCohorts); %step percentile increase
        startPrc=0+stepPrc/2;
        endPrc=100-stepPrc/2;
        for ii=1:numel(source)
            subplot(1,2,ii) 
            curSource = source{ii};
            if strcmp(curSource,'ampTotal')
                ciColor = [153/255 255/255 51/255];
                meanColor = [0/255 102/255 0];
            elseif strcmp(curSource,'forceMag')
                ciColor = [255/255 153/255 153/255];
                meanColor = [153/255 0/255 0];
            else
                ciColor = [153/255 255/255 51/255];
                meanColor = [0/255 102/255 0];
            end
    %         prevLT=0;
            for curPrcLT=startPrc:stepPrc:endPrc
                curLT = floor(prctile(lifeTime,curPrcLT));
                upperLT = prctile(lifeTime,curPrcLT+stepPrc/2);
                lowerLT = prctile(lifeTime,curPrcLT-stepPrc/2);
                % get related track profiles
                curTrackIDsWithCurLT = arrayfun(@(x) x.lifeTime>lowerLT & x.lifeTime<=upperLT,tracksNA);
    %             arrayfun(@(x) plot(0:x.endingFrameExtra-x.startingFrameExtra,x.ampTotal(x.startingFrameExtra:x.endingFrameExtra)),tracksNA(curTrackIDsWithCurLT))
                curArray = NaN(sum(curTrackIDsWithCurLT),curLT+1+2*prePostFramesUsed);
                pp=0;
                for kk=find(curTrackIDsWithCurLT)'
                    pp=pp+1;
                    x = tracksNA(kk);
                    try
                        sF = max(x.startingFrameExtra-prePostFramesUsed,x.startingFrameExtraExtra);
                        eF = min(x.endingFrameExtra+prePostFramesUsed,x.endingFrameExtraExtra);
        %                 interpStep = (curLT+2*prePostFramesUsed)/(eF-sF);
                        % We interpolate the series per predetection, detection,
                        % postdetection periods
                        % pre-detection period
                        if x.startingFrameExtra>sF+1
                            interpStepPre = (prePostFramesUsed)/(x.startingFrameExtra-sF);
                            curArray(pp,1:prePostFramesUsed+1)=interp1((1:interpStepPre:prePostFramesUsed+1),...
                                getfield(x,{1},curSource,{sF:x.startingFrameExtra}),1:prePostFramesUsed+1);
                        end
                        % detection period
        %                 interpStep = (curLT+2*prePostFramesUsed)/(eF-sF);
        %                 curArray(pp,1:curLT+2*prePostFramesUsed+1)=interp1((1:interpStep:curLT+2*prePostFramesUsed+1),...
        %                     getfield(x,{1},curSource,{sF:eF}),1:curLT+2*prePostFramesUsed+1);
                        interpStep = (curLT)/(x.endingFrameExtra-x.startingFrameExtra);
                        curArray(pp,1+prePostFramesUsed:curLT+prePostFramesUsed+1)=...
                            interp1((1+prePostFramesUsed:interpStep:curLT+prePostFramesUsed+1),...
                            getfield(x,{1},curSource,{x.startingFrameExtra:x.endingFrameExtra}),1+prePostFramesUsed:curLT+prePostFramesUsed+1);
                        % post-detection period
                        if eF>x.endingFrameExtra+1
                            interpStepPost = (prePostFramesUsed)/(eF-x.endingFrameExtra);
                            curArray(pp,curLT+prePostFramesUsed+1:curLT+2*prePostFramesUsed+1)=...
                                interp1((curLT+prePostFramesUsed+1:interpStepPost:curLT+2*prePostFramesUsed+1),...
                                getfield(x,{1},curSource,{x.endingFrameExtra:eF}),curLT+prePostFramesUsed+1:curLT+2*prePostFramesUsed+1);
                        end
        %                     x.ampTotal(x.startingFrameExtra:x.endingFrameExtra),1:curLT+1);
        %                 curArray(pp,1:(x.endingFrameExtra-x.startingFrameExtra+1))=x.ampTotal(x.startingFrameExtra:x.endingFrameExtra);
                    catch
                        sF = x.startingFrame;
                        eF = x.endingFrame;
                        interpStep = (curLT)/(eF-sF);
                        try
                            curArray(pp,1+prePostFramesUsed:curLT+prePostFramesUsed+1)=...
                                interp1((1+prePostFramesUsed:interpStep:curLT+prePostFramesUsed+1),...
                                getfield(x,{1},curSource,{sF:eF}),1+prePostFramesUsed:curLT+prePostFramesUsed+1);
                        catch
                            disp(['No ' curSource ' found in ' num2str(kk) 'th track. Skipping...'])
                        end
                   end
                end
    %             if strcmp(curSource,'forceMag')
    %                 % In case of force, there can be some effect from large
    %                 % adhesions. Since we are interested in relative change,
    %                 % I'll shift everything to minimum force
    %             [rowsNonNan,colsNonNan]=find(~isnan(curArray));
    %             [~,uniqNonNan] = unique(rowsNonNan);
    %             rowsNonNanUniq = rowsNonNan(uniqNonNan);
    %             colsNonNanUniq = colsNonNan(uniqNonNan);
    %             linearIndNonNan=sub2ind(size(curArray),rowsNonNanUniq, colsNonNanUniq);
    %             medianStartingSig = nanmedian(curArray(linearIndNonNan));
    %             toBeSubtracted = curArray(linearIndNonNan)-medianStartingSig;
                medianStartingSig = prctile(curArray(:,1),5);
                toBeSubtracted = nanmin(curArray,[],2)-medianStartingSig;
                curArray = curArray - repmat(toBeSubtracted,1,curLT+2*prePostFramesUsed+1);
                meanRangeSig = prctile(prctile(curArray,99,2),90)-prctile(prctile(curArray,1,2),10);
                if yNormalization
                    for jj=1:size(curArray,1)
                        curArray(jj,:) = curArray(jj,:)-nanmin(curArray(jj,:))+medianStartingSig;
    %                     curArray(jj,:) = meanRangeSig*(curArray(jj,:)-nanmin(curArray(jj,:)))/(nanmax(curArray(jj,:))-nanmin(curArray(jj,:)))+medianStartingSig;
                    end
                end
    %             end
                curMeanSig = nanmean(curArray,1);
                curSEM = nanstd(curArray,1)/sqrt(size(curArray,1));
                curTScore = tinv([0.025 0.975],size(curArray,1)-1);
                curCI_upper = curMeanSig + curTScore*curSEM;
                curCI_lower = curMeanSig - curTScore*curSEM;
                if onlyFirstMode
                    allColors =  distinguishable_colors(size(curArray,1),'w');
    %                 co = get(gca,'ColorOrder');
                    set(gca, 'ColorOrder', allColors,'NextPlot','replacechildren')
                    idForceTrans=nanmax(curArray(:,1:prePostFrames),[],2)<nanmax(curArray(:,prePostFrames+1:prePostFrames+1+curLT),[],2);
    %                 plot((-prePostFramesUsed:curLT+prePostFramesUsed)*tInterval_used,curArray(~idForceTrans,:),'Linewidth',0.5,'Color',[0.5 0.5 0.5]), hold on
                    plot((-prePostFramesUsed:curLT+prePostFramesUsed)*tInterval_used,curArray(idForceTrans,:),'Linewidth',0.5), hold on
    %                 plot((-prePostFramesUsed:curLT+prePostFramesUsed)*tInterval_used,curArray,'Linewidth',0.5), hold on
                    plot((-prePostFramesUsed:curLT+prePostFramesUsed)*tInterval_used,curMeanSig,'Linewidth',3,'Color','w')     
                    plot((-prePostFramesUsed:curLT+prePostFramesUsed)*tInterval_used,curMeanSig,'Linewidth',1,'Color',meanColor)     
                else
                    fill([-prePostFramesUsed:curLT+prePostFramesUsed curLT+prePostFramesUsed:-1:-prePostFramesUsed]*tInterval_used,[curCI_upper fliplr(curCI_lower)],ciColor,'EdgeColor',ciColor),hold on
                    plot((-prePostFramesUsed:curLT+prePostFramesUsed)*tInterval_used,curMeanSig,'Linewidth',1,'Color',meanColor)
                end
                if ~isempty(prePostFrames)
                    line([curLT*tInterval_used curLT*tInterval_used],[0 curMeanSig(curLT+prePostFramesUsed+1)],'linestyle',':','Color','k')
                end
    %             prevLT = curLT;
                if exist('curPrcLT','var')
                    if onlyFirstMode
                        maxYamp = max(curArray(:));
                        minYamp = min(curArray(:));
                        break
                    else
                        if curPrcLT==startPrc
                            maxYamp = max(curCI_upper);
                            minYamp = min(curCI_lower);
                        else
                            maxYamp = max(maxYamp,max(curCI_upper));
                            minYamp = min(minYamp,min(curCI_lower));
                        end
                    end
                else
                    maxYamp = max([maxYamp curCI_upper]);
                    minYamp = min([minYamp curCI_lower]);
                end
                if ~isempty(curArray)
                    longestCohorts.rawArray = curArray;
                    longestCohorts.curMeanSig = curMeanSig;
                    longestCohorts.curCI_upper = curCI_upper;
                    longestCohorts.curCI_lower = curCI_lower;
                    longestCohorts.curLT = curLT;
                end
            end
            maxLifeTime = (curLT+prePostFramesUsed+1)*tInterval_used;
            if strcmp(curSource,'ampTotal')
                ylabel('F.I. (a.u.)')
            elseif strcmp(curSource,'forceMag')
                ylabel('Traction (Pa)')
            end
            if isempty(tInterval)
                xlabel('Time (frame)')
            else
                xlabel('Time (sec)')
            end
        end    
    end
end
if ~plotCohorts && ~plotConfInt
    if nTracks<5 %|| alignEvent
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
end

if ~plotCohorts && ~plotConfInt
    nEachFrame = sum(AmpArray>0,1);
    startAmpAvgFrame=find(nEachFrame>size(AmpArray,1)*0.02,1);
    if nTracks>1 && strcmp(source{1},'ampTotal') && ~strcmp(source{end},'edgeAdvanceDist')
        subplot(1,2,1), plot(startAmpAvgFrame:nSampleFrames,AmpG1avg(startAmpAvgFrame:nSampleFrames),'k','Linewidth',3)
    end
    if nTracks>1 && strcmp(source{2},'forceMag') && ~strcmp(source{end},'edgeAdvanceDist')
        subplot(1,2,2), plot(startAmpAvgFrame:nSampleFrames,forceG1avg(startAmpAvgFrame:nSampleFrames),'r','Linewidth',3)
    end
end
xlim([-1-prePostFramesUsed maxLifeTime]*tInterval_used)
ylim([minYamp maxYamp])
% if strcmp(source{1},'forceMag')
%     ylim([0 maxYamp])
% end
set(gca,'FontSize',7)

if ~useCurrentAxis
    if strcmp(source,'edgeAdvanceDist')
        subplot(1,3,2)
    else
        subplot(1,2,2) 
    end
end
% if ~plotCohorts && ~plotConfInt && nTracks>1 && (strcmp(source{1},'forceMag') || (length(source)>1 && strcmp(source{2},'forceMag')))
%     plot(1:nSampleFrames,forceG1avg,'r','Linewidth',3)
% end
if ~useCurrentAxis && ~plotCohorts
    xlim([-1-prePostFramesUsed maxLifeTime]*tInterval_used)
%     xlim([0 maxLifeTime])
    maxYforce = quantile(nanmax(forceArray),0.95);
    minYforce = quantile(nanmin(forceArray),0.01);
    ylim([minYforce maxYforce])
    if isempty(tInterval)
        xlabel('Time (frame)')
    else
        xlabel('Time (sec)')
    end
    ylabel('Traction (Pa)')
    set(gca,'FontSize',7)
end
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


if ~isempty(fileStore)
    [pathStore,nameStore]=fileparts(fileStore);
    if ~exist(pathStore,'dir')
        mkdir(pathStore)
    end
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2', '-r300', strcat(fileStore));
    savefig(fullfile(pathStore,[nameStore '.fig']));
end
