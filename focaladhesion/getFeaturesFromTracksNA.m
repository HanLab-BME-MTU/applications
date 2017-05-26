function tracksNA = getFeaturesFromTracksNA(tracksNA,deltaT,getEdgeRelatedFeatures,cropMaskStack)

disp('Post-analysis on adhesion movement...')
tIntervalMin = deltaT/60; % in min
periodMin = 1;
periodFrames = floor(periodMin/tIntervalMin); % early period in frames
frames2min = floor(2/tIntervalMin); % early period in frames
numTracks = numel(tracksNA);
progressText(0,'Post-analysis:');
for k=1:numTracks
    % cross-correlation scores
%     presIdx = logical(tracksNA(k).presence);
    % get the instantaneous velocity
    % get the distance first
%     curTrackVelLength=sum(presIdx)-1;
%     presIdxSeq = find(presIdx);
    presIdxSeq = tracksNA(k).startingFrameExtra:tracksNA(k).endingFrameExtra;
    curTrackVelLength=length(presIdxSeq)-1;
    distTrajec=zeros(curTrackVelLength,1);
    if getEdgeRelatedFeatures
        for kk=1:curTrackVelLength
            real_kk = presIdxSeq(kk);
            distTrajec(kk) = sqrt(sum((tracksNA(k).closestBdPoint(real_kk+1,:)- ...
                tracksNA(k).closestBdPoint(real_kk,:)).^2,2));
            lastPointIntX = round(tracksNA(k).closestBdPoint(real_kk+1,1));
            lastPointIntY = round(tracksNA(k).closestBdPoint(real_kk+1,2));
            if cropMaskStack(lastPointIntY,lastPointIntX,real_kk) %if the last point is in the first mask, it is inward
                distTrajec(kk) = -distTrajec(kk);
            end
        end
        if any(distTrajec~=0)
            [Protrusion,Retraction] = getPersistenceTime(distTrajec,deltaT);%,'plotYes',true)
            if any(isnan(Retraction.persTime)) || sum(Protrusion.persTime) - sum(Retraction.persTime)>0 % this is protrusion for this track
                tracksNA(k).isProtrusion = true;
            else
                tracksNA(k).isProtrusion = false;
            end
            % average velocity (positive for protrusion)
            curProtVel = (Protrusion.Veloc); curProtVel(isnan(curProtVel))=0;
            curProtPersTime = (Protrusion.persTime); curProtPersTime(isnan(curProtPersTime))=0;
            curRetVel = (Retraction.Veloc); curRetVel(isnan(curRetVel))=0;
            curRetPersTime = (Retraction.persTime); curRetPersTime(isnan(curRetPersTime))=0;

            tracksNA(k).edgeVel = (mean(curProtVel.*curProtPersTime)-mean(curRetVel.*curRetPersTime))/mean([curProtPersTime;curRetPersTime]);
        else
            tracksNA(k).edgeVel = 0;
        end
    end
    % lifetime information
    try
        sFextend=tracksNA(k).startingFrameExtraExtra;
        eFextend=tracksNA(k).endingFrameExtraExtra;
        sF=tracksNA(k).startingFrameExtra;
        eF=tracksNA(k).endingFrameExtra;
        if isempty(sF)
            sF=tracksNA(k).startingFrame;
        end
        if isempty(eF)
            eF=tracksNA(k).endingFrame;
        end
    catch
        sF=tracksNA(k).startingFrame;
        eF=tracksNA(k).endingFrame;
        tracksNA(k).lifeTime = eF-sF;    
    end
    tracksNA(k).lifeTime = eF-sF;    
    % Inital intensity slope for one min
    timeInterval = deltaT/60; % in min
    earlyPeriod = floor(1/timeInterval); % frames per minute
    lastFrame = min(sum(~isnan(tracksNA(k).amp)),sF+earlyPeriod-1);
    lastFrameFromOne = lastFrame - sF+1;
%     lastFrameFromOne = sF;
%     lastFrame = min(sum(~isnan(tracksNA(k).amp)),sF+earlyPeriod-1);
    [curR,curM] = regression(timeInterval*(1:lastFrameFromOne),tracksNA(k).amp(sF:lastFrame));
    tracksNA(k).ampSlope = curM; % in a.u./min
    tracksNA(k).ampSlopeR = curR; % Pearson's correlation coefficient
    
    curEndFrame = min(tracksNA(k).startingFrameExtra+periodFrames-1,tracksNA(k).endingFrame);
    curEarlyPeriod = curEndFrame - tracksNA(k).startingFrameExtra+1;
    [~,curM] = regression(tIntervalMin*(1:curEarlyPeriod),tracksNA(k).ampTotal(tracksNA(k).startingFrameExtra:curEndFrame));
    tracksNA(k).earlyAmpSlope = curM; % in a.u./min

    % Assembly rate: Slope from emergence to maximum - added 10/27/2016 for
    % Michelle's analysis % This output might contain an error when the
    % value has noisy maximum. So it should be filtered with
    % earlyAmpSlope..
    splineParam=0.01; 
    tRange = tracksNA(k).startingFrameExtra:tracksNA(k).endingFrameExtra;
    curAmpTotal =  tracksNA(k).ampTotal(tRange);
    sd_spline= csaps(tRange,curAmpTotal,splineParam);
    sd=ppval(sd_spline,tRange);
    % Find the maximum
    [~,maxSdInd] = max(sd);
    maxAmpFrame = tRange(maxSdInd);
    
%     [~,assemRate] = regression(tIntervalMin*(tRange(1:maxSdInd)),tracksNA(k).ampTotal(tracksNA(k).startingFrameExtra:maxAmpFrame));
    nSampleStart=min(9,floor((maxSdInd)/2));
    if nSampleStart>4 && ttest2(curAmpTotal(1:nSampleStart),curAmpTotal(maxSdInd-nSampleStart+1:maxSdInd)) && ...
            mean(curAmpTotal(1:nSampleStart))<mean(curAmpTotal(maxSdInd-nSampleStart+1:maxSdInd))
        [~,assemRate] = regression(tIntervalMin*(tRange(1:maxSdInd)),...
            log(tracksNA(k).ampTotal(tracksNA(k).startingFrameExtra:maxAmpFrame)/...
            tracksNA(k).ampTotal(tracksNA(k).startingFrameExtra)));
    else
        assemRate = NaN;
    end
    tracksNA(k).assemRate = assemRate; % in 1/min
    
    % Disassembly rate: Slope from maximum to end
%     [~,disassemRate] = regression(tIntervalMin*(tRange(maxSdInd:end)),tracksNA(k).ampTotal(maxAmpFrame:tracksNA(k).endingFrameExtra));
    % I decided to exclude tracks whose end point amplitude is still
    % hanging, i.e., ending amplitude is still high enough compared to
    % starting point, or the last 10 points are not different from 10
    % poihnts near the maximum
    nSampleEnd=min(9,floor((length(tRange)-maxSdInd)*2/3));
    if nSampleEnd>4 && ttest2(curAmpTotal(end-nSampleEnd:end),curAmpTotal(maxSdInd:maxSdInd+nSampleEnd)) && ...
            mean(curAmpTotal(end-nSampleEnd:end))<mean(curAmpTotal(maxSdInd:maxSdInd+nSampleEnd))
        [~,disassemRate] = regression(tIntervalMin*(tRange(maxSdInd:end)),...
            log(curAmpTotal(maxSdInd) ./curAmpTotal(maxSdInd:end)));
    else
        disassemRate = NaN;
    end
    tracksNA(k).disassemRate = disassemRate; % in 1/min

    nSampleEndLate=min(9,floor((tracksNA(k).endingFrameExtraExtra-maxAmpFrame)*2/3));
    curStartFrame = max(tracksNA(k).startingFrame,tracksNA(k).endingFrameExtraExtra-periodFrames+1);
    curLatePeriod = tracksNA(k).endingFrameExtraExtra - curStartFrame+1;
    if nSampleEndLate>4 && ttest2(tracksNA(k).ampTotal(tracksNA(k).endingFrameExtraExtra-nSampleEndLate:tracksNA(k).endingFrameExtraExtra),...
            tracksNA(k).ampTotal(maxAmpFrame:maxAmpFrame+nSampleEndLate)) && ...
            mean(tracksNA(k).ampTotal(tracksNA(k).endingFrameExtraExtra-nSampleEndLate:tracksNA(k).endingFrameExtraExtra))...
            <mean(tracksNA(k).ampTotal(maxAmpFrame:maxAmpFrame+nSampleEndLate))
        [~,curMlate] = regression(tIntervalMin*(1:curLatePeriod),tracksNA(k).ampTotal(curStartFrame:tracksNA(k).endingFrameExtraExtra));
    else
        curMlate = NaN;
    end
    tracksNA(k).lateAmpSlope = curMlate; % in a.u./min
    

    curEndFrame = min(sF+periodFrames-1,eF);
    curEarlyPeriod = curEndFrame - sF+1;
    if getEdgeRelatedFeatures
        [~,curMdist] = regression(tIntervalMin*(1:curEarlyPeriod),tracksNA(k).distToEdge(sF:curEndFrame));
        tracksNA(k).distToEdgeSlope = curMdist; % in a.u./min
        tracksNA(k).distToEdgeChange = (tracksNA(k).distToEdge(end)-tracksNA(k).distToEdge(tracksNA(k).startingFrame)); % in pixel
        % Determining protrusion/retraction based on closestBdPoint and [xCoord
        % yCoord]. If the edge at one time point does not cross the adhesion
        % tracks over entire life time, there is no problem. But that's not
        % always the case: adhesion tracks crosses the cell boundary at first
        % or last time point. And we don't know if the adhesion tracks are
        % moving in the direction of protrusion or retraction. Thus, what I
        % will do is to calculate the inner product of vectors from the
        % adhesion tracks from the closestBdPoint at the first frame. If the
        % product is positive, it means both adhesion points are in the same
        % side. And if distance from the boundary to the last point is larger
        % than the distance to the first point, it means the adhesion track is
        % retracting. In the opposite case, the track is protruding (but it
        % will happen less likely because usually the track would cross the
        % first frame boundary if it is in the protrusion direction). 

        % We need to find out boundary points projected on the line of adhesion track 
        try
            fitobj = fit(tracksNA(k).xCoord(sF:eF)',tracksNA(k).yCoord(sF:eF)','poly1'); % this is an average linear line fit of the adhesion track
        catch
            tracksNA(k)=readIntensityFromTracks(tracksNA(k),imgStack,1,'extraLength',30,'movieData',MD,'retrack',reTrack);
        end
        x0=nanmedian(tracksNA(k).xCoord);
        y0=fitobj(x0);
        dx = 1;
        dy = fitobj.p1;
        trackLine = createLineGeom2d(x0,y0,dx,dy); % this is a geometric average linear line of the adhesion track

    %     trackLine = edgeToLine(edge);
        firstBdPoint = [tracksNA(k).closestBdPoint(sF,1) tracksNA(k).closestBdPoint(sF,2)];
        firstBdPointProjected = projPointOnLine(firstBdPoint, trackLine); % this is an edge boundary point at the first time point projected on the average line of track.
        % try to record advanceDist and edgeAdvanceDist for every single time
        % point ...
        for ii=sFextend:eFextend
            curBdPoint = [tracksNA(k).closestBdPoint(ii,1) tracksNA(k).closestBdPoint(ii,2)];
            curBdPointProjected = projPointOnLine(curBdPoint, trackLine); % this is an edge boundary point at the last time point projected on the average line of track.

            fromFirstBdPointToFirstAdh = [tracksNA(k).xCoord(sF)-firstBdPointProjected(1), tracksNA(k).yCoord(sF)-firstBdPointProjected(2)]; % a vector from the first edge point to the first track point
            fromFirstBdPointToLastAdh = [tracksNA(k).xCoord(ii)-firstBdPointProjected(1), tracksNA(k).yCoord(ii)-firstBdPointProjected(2)]; % a vector from the first edge point to the last track point
            fromCurBdPointToFirstAdh = [tracksNA(k).xCoord(sF)-curBdPointProjected(1), tracksNA(k).yCoord(sF)-curBdPointProjected(2)]; % a vector from the last edge point to the first track point
            fromCurBdPointToLastAdh = [tracksNA(k).xCoord(ii)-curBdPointProjected(1), tracksNA(k).yCoord(ii)-curBdPointProjected(2)]; % a vector from the last edge point to the last track point
            firstBDproduct=fromFirstBdPointToFirstAdh*fromFirstBdPointToLastAdh';
            curBDproduct=fromCurBdPointToFirstAdh*fromCurBdPointToLastAdh';
            if firstBDproduct>0 && firstBDproduct>curBDproduct% both adhesion points are in the same side
                tracksNA(k).advanceDist(ii) = (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
                tracksNA(k).edgeAdvanceDist(ii) = (fromCurBdPointToLastAdh(1)^2 + fromCurBdPointToLastAdh(2)^2)^0.5 - ...
                                                                (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
            else
                if curBDproduct>0 % both adhesion points are in the same side w.r.t. last boundary point
                    tracksNA(k).advanceDist(ii) = (fromCurBdPointToFirstAdh(1)^2 + fromCurBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                    (fromCurBdPointToLastAdh(1)^2 + fromCurBdPointToLastAdh(2)^2)^0.5; % in pixel
                    tracksNA(k).edgeAdvanceDist(ii) = (fromCurBdPointToFirstAdh(1)^2 + fromCurBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                    (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5; % in pixel
                    % this code is nice to check:
        %             figure, imshow(paxImgStack(:,:,tracksNA(k).endingFrame),[]), hold on, plot(tracksNA(k).xCoord,tracksNA(k).yCoord,'w'), plot(tracksNA(k).closestBdPoint(:,1),tracksNA(k).closestBdPoint(:,2),'r')
        %             plot(firstBdPointProjected(1),firstBdPointProjected(2),'co'),plot(lastBdPointProjected(1),lastBdPointProjected(2),'bo')
        %             plot(tracksNA(k).xCoord(tracksNA(k).startingFrame),tracksNA(k).yCoord(tracksNA(k).startingFrame),'yo'),plot(tracksNA(k).xCoord(tracksNA(k).endingFrame),tracksNA(k).yCoord(tracksNA(k).endingFrame),'mo')
                else % Neither products are positive. This means the track crossed both the first and last boundaries. These would show shear movement. Relative comparison is performed.
        %             disp(['Adhesion track ' num2str(k) ' crosses both the first and last boundaries. These would show shear movement. Relative comparison is performed...'])
                    % Using actual BD points instead of projected ones because
                    % somehow the track might be tilted...
                    fromFirstBdPointToFirstAdh = [tracksNA(k).xCoord(sF)-tracksNA(k).closestBdPoint(sF,1), tracksNA(k).yCoord(sF)-tracksNA(k).closestBdPoint(sF,2)];
                    fromFirstBdPointToLastAdh = [tracksNA(k).xCoord(ii)-tracksNA(k).closestBdPoint(sF,1), tracksNA(k).yCoord(ii)-tracksNA(k).closestBdPoint(sF,2)];
                    fromCurBdPointToFirstAdh = [tracksNA(k).xCoord(sF)-tracksNA(k).closestBdPoint(ii,1), tracksNA(k).yCoord(sF)-tracksNA(k).closestBdPoint(ii,2)];
                    fromCurBdPointToLastAdh = [tracksNA(k).xCoord(ii)-tracksNA(k).closestBdPoint(ii,1), tracksNA(k).yCoord(ii)-tracksNA(k).closestBdPoint(ii,2)];
                    firstBDproduct=fromFirstBdPointToFirstAdh*fromFirstBdPointToLastAdh';
                    curBDproduct=fromCurBdPointToFirstAdh*fromCurBdPointToLastAdh';
                    if firstBDproduct>curBDproduct % First BD point is in more distant position from the two adhesion points than the current BD point is.
                        tracksNA(k).advanceDist(ii) = (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                        (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
                        tracksNA(k).edgeAdvanceDist(ii) = (fromCurBdPointToLastAdh(1)^2 + fromCurBdPointToLastAdh(2)^2)^0.5 - ...
                                                                        (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
                    else
                        tracksNA(k).advanceDist(ii) = (fromCurBdPointToFirstAdh(1)^2 + fromCurBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                        (fromCurBdPointToLastAdh(1)^2 + fromCurBdPointToLastAdh(2)^2)^0.5; % in pixel
                        tracksNA(k).edgeAdvanceDist(ii) = (fromCurBdPointToFirstAdh(1)^2 + fromCurBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                        (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5; % in pixel
                    end
                end
            end
        end
        % Record average advanceDist and edgeAdvanceDist for entire lifetime,
        % last 2 min, and every 2 minutes, and maximum of those
        % protrusion/retraction distance to figure out if the adhesion
        % experienced any previous protrusion ...
        for ii=sF:eF
            i2minBefore = max(sF,ii-frames2min);
            tracksNA(k).advanceDistChange2min(ii) = tracksNA(k).advanceDist(ii)-tracksNA(k).advanceDist(i2minBefore);
            tracksNA(k).edgeAdvanceDistChange2min(ii) = tracksNA(k).edgeAdvanceDist(ii)-tracksNA(k).edgeAdvanceDist(i2minBefore);
        end
        % Get the maximum of them. 
        tracksNA(k).maxAdvanceDistChange = max(tracksNA(k).advanceDistChange2min(sF:eF-1));
        tracksNA(k).maxEdgeAdvanceDistChange = max(tracksNA(k).edgeAdvanceDistChange2min(sF:eF-1));
    %     lastBdPoint = [tracksNA(k).closestBdPoint(eF,1) tracksNA(k).closestBdPoint(eF,2)];
    %     lastBdPointProjected = projPointOnLine(lastBdPoint, trackLine); % this is an edge boundary point at the last time point projected on the average line of track.
    %     
    %     fromFirstBdPointToFirstAdh = [tracksNA(k).xCoord(sF)-firstBdPointProjected(1), tracksNA(k).yCoord(sF)-firstBdPointProjected(2)]; % a vector from the first edge point to the first track point
    %     fromFirstBdPointToLastAdh = [tracksNA(k).xCoord(eF)-firstBdPointProjected(1), tracksNA(k).yCoord(eF)-firstBdPointProjected(2)]; % a vector from the first edge point to the last track point
    %     fromLastBdPointToFirstAdh = [tracksNA(k).xCoord(sF)-lastBdPointProjected(1), tracksNA(k).yCoord(sF)-lastBdPointProjected(2)]; % a vector from the last edge point to the first track point
    %     fromLastBdPointToLastAdh = [tracksNA(k).xCoord(eF)-lastBdPointProjected(1), tracksNA(k).yCoord(eF)-lastBdPointProjected(2)]; % a vector from the last edge point to the last track point
    %     firstBDproduct=fromFirstBdPointToFirstAdh*fromFirstBdPointToLastAdh';
    %     lastBDproduct=fromLastBdPointToFirstAdh*fromLastBdPointToLastAdh';
    %     if firstBDproduct>0 && firstBDproduct>lastBDproduct% both adhesion points are in the same side
    %         tracksNA(k).advanceDist = (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5 - ...
    %                                                         (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
    %         tracksNA(k).edgeAdvanceDist = (fromLastBdPointToLastAdh(1)^2 + fromLastBdPointToLastAdh(2)^2)^0.5 - ...
    %                                                         (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
    %     else
    %         if lastBDproduct>0 % both adhesion points are in the same side w.r.t. last boundary point
    %             tracksNA(k).advanceDist = (fromLastBdPointToFirstAdh(1)^2 + fromLastBdPointToFirstAdh(2)^2)^0.5 - ...
    %                                                             (fromLastBdPointToLastAdh(1)^2 + fromLastBdPointToLastAdh(2)^2)^0.5; % in pixel
    %             tracksNA(k).edgeAdvanceDist = (fromLastBdPointToFirstAdh(1)^2 + fromLastBdPointToFirstAdh(2)^2)^0.5 - ...
    %                                                             (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5; % in pixel
    %             % this code is nice to check:
    % %             figure, imshow(paxImgStack(:,:,tracksNA(k).endingFrame),[]), hold on, plot(tracksNA(k).xCoord,tracksNA(k).yCoord,'w'), plot(tracksNA(k).closestBdPoint(:,1),tracksNA(k).closestBdPoint(:,2),'r')
    % %             plot(firstBdPointProjected(1),firstBdPointProjected(2),'co'),plot(lastBdPointProjected(1),lastBdPointProjected(2),'bo')
    % %             plot(tracksNA(k).xCoord(tracksNA(k).startingFrame),tracksNA(k).yCoord(tracksNA(k).startingFrame),'yo'),plot(tracksNA(k).xCoord(tracksNA(k).endingFrame),tracksNA(k).yCoord(tracksNA(k).endingFrame),'mo')
    %         else
    % %             disp(['Adhesion track ' num2str(k) ' crosses both the first and last boundaries. These would show shear movement. Relative comparison is performed...'])
    %             fromFirstBdPointToFirstAdh = [tracksNA(k).xCoord(sF)-tracksNA(k).closestBdPoint(sF,1), tracksNA(k).yCoord(sF)-tracksNA(k).closestBdPoint(sF,2)];
    %             fromFirstBdPointToLastAdh = [tracksNA(k).xCoord(eF)-tracksNA(k).closestBdPoint(sF,1), tracksNA(k).yCoord(eF)-tracksNA(k).closestBdPoint(sF,2)];
    %             fromLastBdPointToFirstAdh = [tracksNA(k).xCoord(sF)-tracksNA(k).closestBdPoint(eF,1), tracksNA(k).yCoord(sF)-tracksNA(k).closestBdPoint(eF,2)];
    %             fromLastBdPointToLastAdh = [tracksNA(k).xCoord(eF)-tracksNA(k).closestBdPoint(eF,1), tracksNA(k).yCoord(eF)-tracksNA(k).closestBdPoint(eF,2)];
    %             firstBDproduct=fromFirstBdPointToFirstAdh*fromFirstBdPointToLastAdh';
    %             lastBDproduct=fromLastBdPointToFirstAdh*fromLastBdPointToLastAdh';
    %             if firstBDproduct>lastBDproduct
    %                 tracksNA(k).advanceDist = (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5 - ...
    %                                                                 (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
    %                 tracksNA(k).edgeAdvanceDist = (fromLastBdPointToLastAdh(1)^2 + fromLastBdPointToLastAdh(2)^2)^0.5 - ...
    %                                                                 (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
    %             else
    %                 tracksNA(k).advanceDist = (fromLastBdPointToFirstAdh(1)^2 + fromLastBdPointToFirstAdh(2)^2)^0.5 - ...
    %                                                                 (fromLastBdPointToLastAdh(1)^2 + fromLastBdPointToLastAdh(2)^2)^0.5; % in pixel
    %                 tracksNA(k).edgeAdvanceDist = (fromLastBdPointToFirstAdh(1)^2 + fromLastBdPointToFirstAdh(2)^2)^0.5 - ...
    %                                                                 (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5; % in pixel
    %             end
    %         end
    %     end
    end
    %Mean Squared Displacement
    meanX = mean(tracksNA(k).xCoord(logical(tracksNA(k).presence)));
    meanY = mean(tracksNA(k).yCoord(logical(tracksNA(k).presence)));
    tracksNA(k).MSD=sum((tracksNA(k).xCoord(logical(tracksNA(k).presence))'-meanX).^2+...
        (tracksNA(k).yCoord(logical(tracksNA(k).presence))'-meanY).^2);
    tracksNA(k).MSDrate = tracksNA(k).MSD/tracksNA(k).lifeTime;
    progressText(k/(numTracks-1));
end
