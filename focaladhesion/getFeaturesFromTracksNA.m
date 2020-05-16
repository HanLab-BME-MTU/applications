function tracksNA = getFeaturesFromTracksNA(tracksNA,deltaT,getEdgeRelatedFeatures,onlyNewAnalysis)

%% input set up
if nargin<4
    onlyNewAnalysis=false;
end

%%
disp('Post-analysis on adhesion movement...')
tIntervalMin = deltaT/60; % in min
periodMin = 1;
periodFrames = floor(periodMin/tIntervalMin); % early period in frames
frames2min = floor(2/tIntervalMin); % early period in frames
numTracks = numel(tracksNA);
progressText(0,'Post-analysis:');
%% run again
if ~isfield(tracksNA,'lifeTime')
    tracksNA(end).lifeTime=[];
end
if ~isfield(tracksNA,'ampSlope')
    tracksNA(end).ampSlope=[];
end
if ~isfield(tracksNA,'earlyAmpSlope')
    tracksNA(end).earlyAmpSlope=[];
end
if ~isfield(tracksNA,'assemRate')
    tracksNA(end).assemRate=[];
end
if ~isfield(tracksNA,'disassemRate')
    tracksNA(end).disassemRate=[];
end
if ~isfield(tracksNA,'lateAmpSlope')
    tracksNA(end).lateAmpSlope=[];
end
if getEdgeRelatedFeatures
    if ~isfield(tracksNA,'distToEdgeSlope')
        tracksNA(end).distToEdgeSlope=[];
    end
    if ~isfield(tracksNA,'distToEdgeChange')
        tracksNA(end).distToEdgeChange=[];
    end
    if ~isfield(tracksNA,'advanceDist')
        tracksNA(end).advanceDist=[];
    end
    if ~isfield(tracksNA,'edgeAdvanceDist')
        tracksNA(end).edgeAdvanceDist=[];
    end
    if ~isfield(tracksNA,'advanceDistNaive')
        tracksNA(end).advanceDistNaive=[];
    end
    if ~isfield(tracksNA,'edgeAdvanceDistNaive')
        tracksNA(end).edgeAdvanceDistNaive=[];
    end
    if ~isfield(tracksNA,'advanceDistChange2min')
        tracksNA(end).advanceDistChange2min=[];
    end
    if ~isfield(tracksNA,'edgeAdvanceDistChange2min')
        tracksNA(end).edgeAdvanceDistChange2min=[];
    end
    if ~isfield(tracksNA,'maxAdvanceDistChange')
        tracksNA(end).maxAdvanceDistChange=[];
    end
    if ~isfield(tracksNA,'maxEdgeAdvanceDistChange')
        tracksNA(end).maxEdgeAdvanceDistChange=[];
    end
    if ~isfield(tracksNA,'MSD')
        tracksNA(end).MSD=[];
    end
    if ~isfield(tracksNA,'MSDrate')
        tracksNA(end).MSDrate=[];
    end
end
parfor k=1:numTracks
    % cross-correlation scores
%     presIdx = logical(curTrack.presence);
    % get the instantaneous velocity
    % get the distance first
%     curTrackVelLength=sum(presIdx)-1;
%     presIdxSeq = find(presIdx);
%     % edgeVel and distTrajec look to be erroneous. deleting it...
%     presIdxSeq = curTrack.startingFrameExtra:curTrack.endingFrameExtra;
%     curTrackVelLength=length(presIdxSeq)-1;
%     distTrajec=zeros(curTrackVelLength,1);
%     if getEdgeRelatedFeatures
%         for kk=1:curTrackVelLength
%             real_kk = presIdxSeq(kk);
%             distTrajec(kk) = sqrt(sum((curTrack.closestBdPoint(real_kk+1,:)- ...
%                 curTrack.closestBdPoint(real_kk,:)).^2,2));
%             lastPointIntX = round(curTrack.closestBdPoint(real_kk+1,1));
%             lastPointIntY = round(curTrack.closestBdPoint(real_kk+1,2));
%             if ~isnan(lastPointIntX) 
%                 if cropMaskStack(lastPointIntY,lastPointIntX,real_kk) %if the last point is in the first mask, it is inward
%                     distTrajec(kk) = -distTrajec(kk);
%                 end
%             end
%         end
%         if any(distTrajec~=0)
%             try
%                 [Protrusion,Retraction] = getPersistenceTime(distTrajec,deltaT);%,'plotYes',true)
%                 if any(isnan(Retraction.persTime)) || sum(Protrusion.persTime) - sum(Retraction.persTime)>0 % this is protrusion for this track
%                     curTrack.isProtrusion = true;
%                 else
%                     curTrack.isProtrusion = false;
%                 end
%                 % average velocity (positive for protrusion)
%                 curProtVel = (Protrusion.Veloc); curProtVel(isnan(curProtVel))=0;
%                 curProtPersTime = (Protrusion.persTime); curProtPersTime(isnan(curProtPersTime))=0;
%                 curRetVel = (Retraction.Veloc); curRetVel(isnan(curRetVel))=0;
%                 curRetPersTime = (Retraction.persTime); curRetPersTime(isnan(curRetPersTime))=0;
% 
%                 curTrack.edgeVel = (mean(curProtVel.*curProtPersTime)-mean(curRetVel.*curRetPersTime))/mean([curProtPersTime;curRetPersTime]);
%             catch
%                 curTrack.edgeVel = 0;
%             end
%         else
%             curTrack.edgeVel = 0;
%         end
%     end
    curTrack=tracksNA(k);
    % lifetime information
    try
        sFextend=curTrack.startingFrameExtraExtra;
        eFextend=curTrack.endingFrameExtraExtra;
        sF=curTrack.startingFrameExtra;
        eF=curTrack.endingFrameExtra;
        if isempty(sF)
            sF=tracksNA(k).startingFrame;
        end
        if isempty(eF)
            eF=tracksNA(k).endingFrame;
        end
    catch
        sF=tracksNA(k).startingFrame;
        eF=tracksNA(k).endingFrame;
    end
    curTrack.lifeTime = eF-sF;   
    tRange = sF:eF;
    % Inital intensity slope for one min
    timeInterval = deltaT/60; % in min
    entirePeriod = eF-sF; % floor(1/timeInterval); % frames per minute
    lastFrame = min(sum(~isnan(curTrack.amp)),sF+entirePeriod-1);
    lastFrameFromOne = lastFrame - sF+1;
    if ~onlyNewAnalysis
    %     lastFrameFromOne = sF;
    %     lastFrame = min(sum(~isnan(curTrack.amp)),sF+earlyPeriod-1);
        [~,curM] = regression(timeInterval*(1:lastFrameFromOne),curTrack.amp(sF:lastFrame));
        curTrack.ampSlope = curM; % in a.u./min
    %     curTrack.ampSlopeR = curR; % Pearson's correlation coefficient

        curEndFrame = min(curTrack.startingFrameExtra+periodFrames-1,curTrack.endingFrame);
        curEarlyPeriod = curEndFrame - curTrack.startingFrameExtra+1;
        [~,curM] = regression(tIntervalMin*(1:curEarlyPeriod),curTrack.ampTotal(curTrack.startingFrameExtra:curEndFrame));
        curTrack.earlyAmpSlope = curM; % in a.u./min

        % Assembly rate: Slope from emergence to maximum - added 10/27/2016 for
        % Michelle's analysis % This output might contain an error when the
        % value has noisy maximum. So it should be filtered with
        % earlyAmpSlope..
        
        % Need to update this calculate this part more wisely
        assemRate = getAssemRateFromTracks(curTrack,tIntervalMin);
        curTrack.assemRate = assemRate; % in 1/min

        % Disassembly rate: Slope from maximum to end
    %     [~,disassemRate] = regression(tIntervalMin*(tRange(maxSdInd:end)),curTrack.ampTotal(maxAmpFrame:curTrack.endingFrameExtra));
        % I decided to exclude tracks whose end point amplitude is still
        % hanging, i.e., ending amplitude is still high enough compared to
        % starting point, or the last 10 points are not different from 10
        % poihnts near the maximum
        disassemRate = getDisassemRateFromTracks(curTrack,tIntervalMin);
        curTrack.disassemRate = disassemRate; % in 1/min
        
        % MaxAmpFrame estimation
        splineParam=0.01; 
        try
            sd_spline= csaps(tRange,curTrack.amp(tRange),splineParam);
        catch
            tRange = find(~isnan(curTrack.amp));
            sd_spline= csaps(tRange,curTrack.amp(tRange),splineParam);
        end
        sd=ppval(sd_spline,tRange);
        [~,maxSdInd] = max(sd);
        maxAmpFrame = tRange(maxSdInd);

        nSampleEndLate=min(9,floor((curTrack.endingFrameExtraExtra-maxAmpFrame)*2/3));
        curStartFrame = max(curTrack.startingFrame,curTrack.endingFrameExtraExtra-periodFrames+1);
        curLatePeriod = curTrack.endingFrameExtraExtra - curStartFrame+1;
        sigTtest=ttest2(curTrack.ampTotal(curTrack.endingFrameExtraExtra-nSampleEndLate:curTrack.endingFrameExtraExtra),...
                curTrack.ampTotal(maxAmpFrame:maxAmpFrame+nSampleEndLate));
        if ~isnan(sigTtest)
            if nSampleEndLate>4 && sigTtest && ...
                    mean(curTrack.ampTotal(curTrack.endingFrameExtraExtra-nSampleEndLate:curTrack.endingFrameExtraExtra))...
                    <mean(curTrack.ampTotal(maxAmpFrame:maxAmpFrame+nSampleEndLate))
                [~,curMlate] = regression(tIntervalMin*(1:curLatePeriod),curTrack.ampTotal(curStartFrame:curTrack.endingFrameExtraExtra));
            else
                curMlate = NaN;
            end
        else
            curMlate = NaN;
        end
        curTrack.lateAmpSlope = curMlate; % in a.u./min
    end

    curEndFrame = min(sF+periodFrames-1,eF);
    curEarlyPeriod = curEndFrame - sF+1;
    if getEdgeRelatedFeatures
        [~,curMdist] = regression(tIntervalMin*(1:curEarlyPeriod),curTrack.distToEdge(sF:curEndFrame));
        curTrack.distToEdgeSlope = curMdist; % in a.u./min
        curTrack.distToEdgeChange = (curTrack.distToEdge(end)-curTrack.distToEdge(curTrack.startingFrame)); % in pixel
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
            fitobj = fit(curTrack.xCoord(sF:eF)',curTrack.yCoord(sF:eF)','poly1'); % this is an average linear line fit of the adhesion track
        catch
            % This means fitting a line with the track has failed, highly
            % likely due to too short track lifetime or too variable
            % lacations or having NaNs
            curX = curTrack.xCoord(sF:eF); curY=curTrack.yCoord(sF:eF);
            indNaNX = isnan(curX);
            t = 1:length(curX);
            t_nn = t(~indNaNX);
            curX2 = interp1(t_nn,curX(~indNaNX),t,'linear','extrap');
            curY2 = interp1(t_nn,curY(~indNaNX),t,'linear','extrap');

            fitobj = fit(curX2(~isnan(curX2))',curY2(~isnan(curY2))','poly1'); % this is an average linear line fit of the adhesion track
        end
        x0=nanmedian(curTrack.xCoord);
        y0=fitobj(x0);
        dx = 1;
        dy = fitobj.p1;
        trackLine = createLineGeom2d(x0,y0,dx,dy); % this is a geometric average linear line of the adhesion track

    %     trackLine = edgeToLine(edge);
        firstBdPoint = [curTrack.closestBdPoint(sF,1) curTrack.closestBdPoint(sF,2)];
        firstBdPointProjected = projPointOnLine(firstBdPoint, trackLine); % this is an edge boundary point at the first time point projected on the average line of track.
        % try to record advanceDist and edgeAdvanceDist for every single time
        % point ...
        for ii=sFextend:eFextend
            if ~onlyNewAnalysis
                curBdPoint = [curTrack.closestBdPoint(ii,1) curTrack.closestBdPoint(ii,2)];
                curBdPointProjected = projPointOnLine(curBdPoint, trackLine); % this is an edge boundary point at the last time point projected on the average line of track.

                fromFirstBdPointToFirstAdh = [curTrack.xCoord(sF)-firstBdPointProjected(1), curTrack.yCoord(sF)-firstBdPointProjected(2)]; % a vector from the first edge point to the first track point
                fromFirstBdPointToLastAdh = [curTrack.xCoord(ii)-firstBdPointProjected(1), curTrack.yCoord(ii)-firstBdPointProjected(2)]; % a vector from the first edge point to the last track point
                fromCurBdPointToFirstAdh = [curTrack.xCoord(sF)-curBdPointProjected(1), curTrack.yCoord(sF)-curBdPointProjected(2)]; % a vector from the last edge point to the first track point
                fromCurBdPointToLastAdh = [curTrack.xCoord(ii)-curBdPointProjected(1), curTrack.yCoord(ii)-curBdPointProjected(2)]; % a vector from the last edge point to the last track point
                firstBDproduct=fromFirstBdPointToFirstAdh*fromFirstBdPointToLastAdh';
                curBDproduct=fromCurBdPointToFirstAdh*fromCurBdPointToLastAdh';

                if firstBDproduct>0 && firstBDproduct>curBDproduct% both adhesion points are in the same side
                    curTrack.advanceDist(ii) = (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                    (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
                    curTrack.edgeAdvanceDist(ii) = (fromCurBdPointToLastAdh(1)^2 + fromCurBdPointToLastAdh(2)^2)^0.5 - ...
                                                                    (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
                else
                    if curBDproduct>0 % both adhesion points are in the same side w.r.t. last boundary point
                        curTrack.advanceDist(ii) = (fromCurBdPointToFirstAdh(1)^2 + fromCurBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                        (fromCurBdPointToLastAdh(1)^2 + fromCurBdPointToLastAdh(2)^2)^0.5; % in pixel
                        curTrack.edgeAdvanceDist(ii) = (fromCurBdPointToFirstAdh(1)^2 + fromCurBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                        (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5; % in pixel
                        % this code is nice to check:
            %             figure, imshow(paxImgStack(:,:,curTrack.endingFrame),[]), hold on, plot(curTrack.xCoord,curTrack.yCoord,'w'), plot(curTrack.closestBdPoint(:,1),curTrack.closestBdPoint(:,2),'r')
            %             plot(firstBdPointProjected(1),firstBdPointProjected(2),'co'),plot(lastBdPointProjected(1),lastBdPointProjected(2),'bo')
            %             plot(curTrack.xCoord(curTrack.startingFrame),curTrack.yCoord(curTrack.startingFrame),'yo'),plot(curTrack.xCoord(curTrack.endingFrame),curTrack.yCoord(curTrack.endingFrame),'mo')
                    else % Neither products are positive. This means the track crossed both the first and last boundaries. These would show shear movement. Relative comparison is performed.
            %             disp(['Adhesion track ' num2str(k) ' crosses both the first and last boundaries. These would show shear movement. Relative comparison is performed...'])
                        % Using actual BD points instead of projected ones because
                        % somehow the track might be tilted...
                        fromFirstBdPointToFirstAdh = [curTrack.xCoord(sF)-curTrack.closestBdPoint(sF,1), curTrack.yCoord(sF)-curTrack.closestBdPoint(sF,2)];
                        fromFirstBdPointToLastAdh = [curTrack.xCoord(ii)-curTrack.closestBdPoint(sF,1), curTrack.yCoord(ii)-curTrack.closestBdPoint(sF,2)];
                        fromCurBdPointToFirstAdh = [curTrack.xCoord(sF)-curTrack.closestBdPoint(ii,1), curTrack.yCoord(sF)-curTrack.closestBdPoint(ii,2)];
                        fromCurBdPointToLastAdh = [curTrack.xCoord(ii)-curTrack.closestBdPoint(ii,1), curTrack.yCoord(ii)-curTrack.closestBdPoint(ii,2)];
                        firstBDproduct=fromFirstBdPointToFirstAdh*fromFirstBdPointToLastAdh';
                        curBDproduct=fromCurBdPointToFirstAdh*fromCurBdPointToLastAdh';
                        if firstBDproduct>curBDproduct % First BD point is in more distant position from the two adhesion points than the current BD point is.
                            curTrack.advanceDist(ii) = (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                            (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
                            curTrack.edgeAdvanceDist(ii) = (fromCurBdPointToLastAdh(1)^2 + fromCurBdPointToLastAdh(2)^2)^0.5 - ...
                                                                            (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
                        else
                            curTrack.advanceDist(ii) = (fromCurBdPointToFirstAdh(1)^2 + fromCurBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                            (fromCurBdPointToLastAdh(1)^2 + fromCurBdPointToLastAdh(2)^2)^0.5; % in pixel
                            curTrack.edgeAdvanceDist(ii) = (fromCurBdPointToFirstAdh(1)^2 + fromCurBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                            (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5; % in pixel
                        end
                    end
                end
            end
            % Now I'll copy everything to advanceDistNaive and
            % edgeAdvanceDistNaive -Nov 2018 Sangyoon
            curBdPointNaive = [curTrack.closestBdPointNaive(ii,1) curTrack.closestBdPoint(ii,2)];
            curBdPointProjected = projPointOnLine(curBdPointNaive, trackLine); % this is an edge boundary point at the last time point projected on the average line of track.

            fromFirstBdPointToFirstAdh = [curTrack.xCoord(sF)-firstBdPointProjected(1), curTrack.yCoord(sF)-firstBdPointProjected(2)]; % a vector from the first edge point to the first track point
            fromFirstBdPointToLastAdh = [curTrack.xCoord(ii)-firstBdPointProjected(1), curTrack.yCoord(ii)-firstBdPointProjected(2)]; % a vector from the first edge point to the last track point
            fromCurBdPointToFirstAdh = [curTrack.xCoord(sF)-curBdPointProjected(1), curTrack.yCoord(sF)-curBdPointProjected(2)]; % a vector from the last edge point to the first track point
            fromCurBdPointToLastAdh = [curTrack.xCoord(ii)-curBdPointProjected(1), curTrack.yCoord(ii)-curBdPointProjected(2)]; % a vector from the last edge point to the last track point
            firstBDproduct=fromFirstBdPointToFirstAdh*fromFirstBdPointToLastAdh';
            curBDproduct=fromCurBdPointToFirstAdh*fromCurBdPointToLastAdh';
            
            if firstBDproduct>0 && firstBDproduct>curBDproduct% both adhesion points are in the same side
                curTrack.advanceDistNaive(ii) = (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
                curTrack.edgeAdvanceDistNaive(ii) = (fromCurBdPointToLastAdh(1)^2 + fromCurBdPointToLastAdh(2)^2)^0.5 - ...
                                                                (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
            else
                if curBDproduct>0 % both adhesion points are in the same side w.r.t. last boundary point
                    curTrack.advanceDistNaive(ii) = (fromCurBdPointToFirstAdh(1)^2 + fromCurBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                    (fromCurBdPointToLastAdh(1)^2 + fromCurBdPointToLastAdh(2)^2)^0.5; % in pixel
                    curTrack.edgeAdvanceDistNaive(ii) = (fromCurBdPointToFirstAdh(1)^2 + fromCurBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                    (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5; % in pixel
                    % this code is nice to check:
        %             figure, imshow(paxImgStack(:,:,curTrack.endingFrame),[]), hold on, plot(curTrack.xCoord,curTrack.yCoord,'w'), plot(curTrack.closestBdPoint(:,1),curTrack.closestBdPoint(:,2),'r')
        %             plot(firstBdPointProjected(1),firstBdPointProjected(2),'co'),plot(lastBdPointProjected(1),lastBdPointProjected(2),'bo')
        %             plot(curTrack.xCoord(curTrack.startingFrame),curTrack.yCoord(curTrack.startingFrame),'yo'),plot(curTrack.xCoord(curTrack.endingFrame),curTrack.yCoord(curTrack.endingFrame),'mo')
                else % Neither products are positive. This means the track crossed both the first and last boundaries. These would show shear movement. Relative comparison is performed.
        %             disp(['Adhesion track ' num2str(k) ' crosses both the first and last boundaries. These would show shear movement. Relative comparison is performed...'])
                    % Using actual BD points instead of projected ones because
                    % somehow the track might be tilted...
                    fromFirstBdPointToFirstAdh = [curTrack.xCoord(sF)-curTrack.closestBdPointNaive(sF,1), curTrack.yCoord(sF)-curTrack.closestBdPointNaive(sF,2)];
                    fromFirstBdPointToLastAdh = [curTrack.xCoord(ii)-curTrack.closestBdPointNaive(sF,1), curTrack.yCoord(ii)-curTrack.closestBdPointNaive(sF,2)];
                    fromCurBdPointToFirstAdh = [curTrack.xCoord(sF)-curTrack.closestBdPointNaive(ii,1), curTrack.yCoord(sF)-curTrack.closestBdPointNaive(ii,2)];
                    fromCurBdPointToLastAdh = [curTrack.xCoord(ii)-curTrack.closestBdPointNaive(ii,1), curTrack.yCoord(ii)-curTrack.closestBdPointNaive(ii,2)];
                    firstBDproduct=fromFirstBdPointToFirstAdh*fromFirstBdPointToLastAdh';
                    curBDproduct=fromCurBdPointToFirstAdh*fromCurBdPointToLastAdh';
                    if firstBDproduct>curBDproduct % First BD point is in more distant position from the two adhesion points than the current BD point is.
                        curTrack.advanceDistNaive(ii) = (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                        (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
                        curTrack.edgeAdvanceDistNaive(ii) = (fromCurBdPointToLastAdh(1)^2 + fromCurBdPointToLastAdh(2)^2)^0.5 - ...
                                                                        (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
                    else
                        curTrack.advanceDistNaive(ii) = (fromCurBdPointToFirstAdh(1)^2 + fromCurBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                        (fromCurBdPointToLastAdh(1)^2 + fromCurBdPointToLastAdh(2)^2)^0.5; % in pixel
                        curTrack.edgeAdvanceDistNaive(ii) = (fromCurBdPointToFirstAdh(1)^2 + fromCurBdPointToFirstAdh(2)^2)^0.5 - ...
                                                                        (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5; % in pixel
                    end
                end
            end
        end
        if ~onlyNewAnalysis
            % Record average advanceDist and edgeAdvanceDist for entire lifetime,
            % last 2 min, and every 2 minutes, and maximum of those
            % protrusion/retraction distance to figure out if the adhesion
            % experienced any previous protrusion ...
            for ii=sF:eF
                i2minBefore = max(sF,ii-frames2min);
                curTrack.advanceDistChange2min(ii) = curTrack.advanceDist(ii)-curTrack.advanceDist(i2minBefore);
                curTrack.edgeAdvanceDistChange2min(ii) = curTrack.edgeAdvanceDist(ii)-curTrack.edgeAdvanceDist(i2minBefore);
            end
            % Get the maximum of them. 
            curTrack.maxAdvanceDistChange = max(curTrack.advanceDistChange2min(sF:eF-1));
            curTrack.maxEdgeAdvanceDistChange = max(curTrack.edgeAdvanceDistChange2min(sF:eF-1));
        end
    %     lastBdPoint = [curTrack.closestBdPoint(eF,1) curTrack.closestBdPoint(eF,2)];
    %     lastBdPointProjected = projPointOnLine(lastBdPoint, trackLine); % this is an edge boundary point at the last time point projected on the average line of track.
    %     
    %     fromFirstBdPointToFirstAdh = [curTrack.xCoord(sF)-firstBdPointProjected(1), curTrack.yCoord(sF)-firstBdPointProjected(2)]; % a vector from the first edge point to the first track point
    %     fromFirstBdPointToLastAdh = [curTrack.xCoord(eF)-firstBdPointProjected(1), curTrack.yCoord(eF)-firstBdPointProjected(2)]; % a vector from the first edge point to the last track point
    %     fromLastBdPointToFirstAdh = [curTrack.xCoord(sF)-lastBdPointProjected(1), curTrack.yCoord(sF)-lastBdPointProjected(2)]; % a vector from the last edge point to the first track point
    %     fromLastBdPointToLastAdh = [curTrack.xCoord(eF)-lastBdPointProjected(1), curTrack.yCoord(eF)-lastBdPointProjected(2)]; % a vector from the last edge point to the last track point
    %     firstBDproduct=fromFirstBdPointToFirstAdh*fromFirstBdPointToLastAdh';
    %     lastBDproduct=fromLastBdPointToFirstAdh*fromLastBdPointToLastAdh';
    %     if firstBDproduct>0 && firstBDproduct>lastBDproduct% both adhesion points are in the same side
    %         curTrack.advanceDist = (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5 - ...
    %                                                         (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
    %         curTrack.edgeAdvanceDist = (fromLastBdPointToLastAdh(1)^2 + fromLastBdPointToLastAdh(2)^2)^0.5 - ...
    %                                                         (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
    %     else
    %         if lastBDproduct>0 % both adhesion points are in the same side w.r.t. last boundary point
    %             curTrack.advanceDist = (fromLastBdPointToFirstAdh(1)^2 + fromLastBdPointToFirstAdh(2)^2)^0.5 - ...
    %                                                             (fromLastBdPointToLastAdh(1)^2 + fromLastBdPointToLastAdh(2)^2)^0.5; % in pixel
    %             curTrack.edgeAdvanceDist = (fromLastBdPointToFirstAdh(1)^2 + fromLastBdPointToFirstAdh(2)^2)^0.5 - ...
    %                                                             (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5; % in pixel
    %             % this code is nice to check:
    % %             figure, imshow(paxImgStack(:,:,curTrack.endingFrame),[]), hold on, plot(curTrack.xCoord,curTrack.yCoord,'w'), plot(curTrack.closestBdPoint(:,1),curTrack.closestBdPoint(:,2),'r')
    % %             plot(firstBdPointProjected(1),firstBdPointProjected(2),'co'),plot(lastBdPointProjected(1),lastBdPointProjected(2),'bo')
    % %             plot(curTrack.xCoord(curTrack.startingFrame),curTrack.yCoord(curTrack.startingFrame),'yo'),plot(curTrack.xCoord(curTrack.endingFrame),curTrack.yCoord(curTrack.endingFrame),'mo')
    %         else
    % %             disp(['Adhesion track ' num2str(k) ' crosses both the first and last boundaries. These would show shear movement. Relative comparison is performed...'])
    %             fromFirstBdPointToFirstAdh = [curTrack.xCoord(sF)-curTrack.closestBdPoint(sF,1), curTrack.yCoord(sF)-curTrack.closestBdPoint(sF,2)];
    %             fromFirstBdPointToLastAdh = [curTrack.xCoord(eF)-curTrack.closestBdPoint(sF,1), curTrack.yCoord(eF)-curTrack.closestBdPoint(sF,2)];
    %             fromLastBdPointToFirstAdh = [curTrack.xCoord(sF)-curTrack.closestBdPoint(eF,1), curTrack.yCoord(sF)-curTrack.closestBdPoint(eF,2)];
    %             fromLastBdPointToLastAdh = [curTrack.xCoord(eF)-curTrack.closestBdPoint(eF,1), curTrack.yCoord(eF)-curTrack.closestBdPoint(eF,2)];
    %             firstBDproduct=fromFirstBdPointToFirstAdh*fromFirstBdPointToLastAdh';
    %             lastBDproduct=fromLastBdPointToFirstAdh*fromLastBdPointToLastAdh';
    %             if firstBDproduct>lastBDproduct
    %                 curTrack.advanceDist = (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5 - ...
    %                                                                 (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
    %                 curTrack.edgeAdvanceDist = (fromLastBdPointToLastAdh(1)^2 + fromLastBdPointToLastAdh(2)^2)^0.5 - ...
    %                                                                 (fromFirstBdPointToLastAdh(1)^2 + fromFirstBdPointToLastAdh(2)^2)^0.5; % in pixel
    %             else
    %                 curTrack.advanceDist = (fromLastBdPointToFirstAdh(1)^2 + fromLastBdPointToFirstAdh(2)^2)^0.5 - ...
    %                                                                 (fromLastBdPointToLastAdh(1)^2 + fromLastBdPointToLastAdh(2)^2)^0.5; % in pixel
    %                 curTrack.edgeAdvanceDist = (fromLastBdPointToFirstAdh(1)^2 + fromLastBdPointToFirstAdh(2)^2)^0.5 - ...
    %                                                                 (fromFirstBdPointToFirstAdh(1)^2 + fromFirstBdPointToFirstAdh(2)^2)^0.5; % in pixel
    %             end
    %         end
    %     end
    end
    if ~onlyNewAnalysis    
        %Mean Squared Displacement
        meanX = mean(curTrack.xCoord(logical(curTrack.presence)));
        meanY = mean(curTrack.yCoord(logical(curTrack.presence)));
        curTrack.MSD=sum((curTrack.xCoord(logical(curTrack.presence))'-meanX).^2+...
            (curTrack.yCoord(logical(curTrack.presence))'-meanY).^2);
        curTrack.MSDrate = curTrack.MSD/curTrack.lifeTime;
    end
    progressText(k/(numTracks-1));
    tracksNA(k) = curTrack;
end
