function [idlist,maxTagIdx] = linkerSecondLAP(idlist, goodTimes, goodTimesM, nSpots, intList, tagIndices, maxTagIdx, constants, verbose, intAx)
%linkerSecondLAP links the source-frames (with the "good" number of spots) to the spots in the target frames
% idlist is sorted at the end of this function, and linkup/linkdown are
% added, and zero-rows are filled in for the future estimated spots

% build "trackingStrategy" to link tags with "wrong" n.
% trackingStrategy will be a numSources+1-by-nTargets matrix.

% find future "target" frames
goodTimesT = setdiff(goodTimes,goodTimesM);


% check whether there is anything to link still, or whether all is well
% already. If all timepoints are in the list of timepoints with the
% correct number of spots, we're perfectly happy - however, we still
% have to sort the linklists and to add linkup/linkdown
if ~isempty(goodTimesT)

    nTargets = length(goodTimesT);
    nSources = min(constants.numSources,length(goodTimesM));
    linkingStrategy = zeros(nSources+1,nTargets);

    % strategy: closest tags. If ambiguous, take the earlier one - its SNR
    % should be lower
    for iTarget = 1:nTargets
        deltaT = abs(goodTimesM - goodTimesT(iTarget));
        [deltaSorted, deltaIdx] = sort(deltaT);
        % take top n sources
        linkingStrategy(:,iTarget) = [iTarget; goodTimesM(deltaIdx(1:nSources))];
    end

    % allocate arrays for delta amplitude and delta position
    [deltaAmp, deltaXyz] = deal(zeros(sum(nSources * nSpots(goodTimesT) * maxTagIdx),1));
    deltaCounter = 0;

    % link according to linkingStrategy. Sum the costMatrices from the
    % different sources, but be careful about potential births/deaths
    for t = linkingStrategy

        % build costMatrix. We cannot correct by a COM, unfortunately, because
        % we do not know by what we should correct. However, we can correct the
        % amplitudes by dividing the amplitudes of the sources. Then we need to
        % go through the same routine as for the source-linking to find the
        % proper sigmas so that we can weigh distance and amplitudes.
        tTarget = goodTimesT(t(1));
        targetAmp = idlist(tTarget).linklist(:,8);
        targetXyz = idlist(tTarget).linklist(:,9:11) - ...
            repmat(idlist(tTarget).centroid  .* constants.useCOM, nSpots(tTarget),1);
        
        % admittedly, this part of the function is a mess and would really
        % benefit from a bit of clean-up.
        % Basically, what the function needs to do is to fill distMats
        % where the rows correspond to all the tags that appear in the
        % source frames, and the columns to the spot indices of the target
        % frame. Where there is no data (which happens if there is a
        % source-tag that doesn't appear in all sources), there should be
        % NaN.
        
        

        % get number of different indices we have in the source frames
        sourceIdxList = unique(tagIndices(t(2:end),:));
        sourceIdxList(~sourceIdxList) = []; % remove zero
        minSourceIdxMinusOne = sourceIdxList(1)-1; % unique sorts

        %     if verbose
        %         figure
        %         plot3(targetXyz(1),targetXyz(2),targetXyz(3),'.r');
        %         hold on
        %         grid on
        %         posAx = gca;
        %     end

        % at the end we will have distMats again that will be stored in the
        % target frames
        nSourceSpots = max(sourceIdxList)-minSourceIdxMinusOne;
        distMatAmp = repmat(NaN,...
            [nSourceSpots,...
            nSpots(goodTimesT(t(1))),...
            nSources]);
        distMatXyz = repmat(NaN,...
            [nSourceSpots,...
            nSpots(goodTimesT(t(1))),...
            nSources]);

        for sourceIdx = 1:nSources
            sourceCoord = [];
            sourceAmp = zeros(nSourceSpots,1);
            sourceXyz = zeros(nSourceSpots,3);

            tSource = t(sourceIdx + 1);
            ampFact = ...
                exp(idlist(1).stats.intFit.xFit(end) * ...
                diff([goodTimesT(t(1)),tSource]));
            % only have entries in sourceCoords wherever we have a tagIdx
            sourceAmp(idlist(tSource).linklist(:,4) - minSourceIdxMinusOne,1) = ...
                [idlist(tSource).linklist(:, 8)./ampFact];
            sourceXyz(idlist(tSource).linklist(:,4) - minSourceIdxMinusOne,1:3) = ...
                idlist(tSource).linklist(:, 9:11) - ...
                repmat(idlist(tSource).centroid, nSpots(tSource),1) * constants.useCOM;
            %         if verbose
            %             markers = '*dph';
            %             for i=1:nSpots(tSource)
            %             plot3(posAx,[sourceXyz(i,1);targetXyz(1)],...
            %                 [sourceXyz(i,2);targetXyz(2)],...
            %                 [sourceXyz(i,3);targetXyz(3)],...
            %                 'Marker',markers(i),...
            %                 'Color',extendedColors(sourceIdx+1),...
            %                 'DisplayName',num2str(tSource));
            %             end
            %         end

            % it is possible that there are more indices than maxNumSpots!
            % also, there could be less. Therefore, use tagIndicesIdx
            % furthermore, sort tagIndices. 
            tagIndicesIdx = find(tagIndices(tSource,:));
            distMatAmp(sort(tagIndices(tSource,tagIndicesIdx)),:,sourceIdx) = ...
                distMat2(sourceAmp(find(sourceAmp)),targetAmp);
            distMatXyz(sort(tagIndices(tSource,tagIndicesIdx)),:,sourceIdx) = ...
                distMat2(sourceXyz(find(sourceAmp),:),targetXyz);

            % instead of settingNaNs, we can just initialize the matrices
            % with NaN (I hope)
%             % set NaN wherever sourceAmp had a zero-row - there should
%             % NEVER be a spot with no amplitude.
%             badIdx = ~distMatAmp(:,:,sourceIdx);
%             distMatAmp(badIdx,:,sourceIdx) = NaN;
%             distMatXyz(badIdx,:,sourceIdx) = NaN;
        end

        % sstore values in lists so that we can
        % build the new sigmas
        deltaCounterNew = deltaCounter + numel(distMatAmp);
        deltaAmp(deltaCounter + 1: deltaCounterNew) = distMatAmp(:);
        deltaXyz(deltaCounter + 1: deltaCounterNew) = distMatXyz(:);
        deltaCounter = deltaCounterNew;

        % the final costMatrices are the nanMean along the third dimension.
        % This takes into account the NaN-rows.
        distMatAmp = nanmean(distMatAmp,3);
        distMatXyz = nanmean(distMatXyz,3);

        % there could be rows resulting from an index that doesn't appear in
        % any of the sources. Remove these NaN-cols!
        distMatAmp(isnan(distMatAmp(:,1)),:) = [];
        distMatXyz(isnan(distMatXyz(:,1)),:) = [];

        % store values in idlist
        idlist(goodTimesT(t(1))).distMatAmp = distMatAmp;
        idlist(goodTimesT(t(1))).distMatXyz = distMatXyz;

        % store sourceIdxList so that we'll know what linked to what ;)
        idlist(goodTimesT(t(1))).sourceIdxList = sourceIdxList;

    end % for linkingStrategy: build costMatrices

    % deltaLists could be a bit too long - we initialized them assuming that
    % we would always get the maximum number of possible tag indices
    deltaAmp = deltaAmp(1:deltaCounter);
    deltaXyz = deltaXyz(1:deltaCounter);

    % if there are more tags than spots, there can be many zeros in deltaAmp,
    % deltaXYZ
    deltaAmp(isnan(deltaAmp)) = [];
    deltaXyz(isnan(deltaXyz)) = [];
    deltaAmp(~deltaAmp) = [];
    deltaXyz(~deltaXyz) = [];

    % find sigmaAmpT, sigmaXyzT
    % find cutoff
    [dummy, maxAmp] = cutFirstHistMode(deltaAmp,verbose>1);
    if verbose > 1
        set(gcf,'Name','Amplitude Histogram S->T')
    end
    [dummy, maxXyz] = cutFirstHistMode(deltaXyz,verbose>1);
    if verbose > 1
        set(gcf,'Name','Displacement Histogram S->T')
    end

    % the average deltaAmplitude is the mean of all deltaAmps that are below
    % the cutoff
    sigmaAmpT = mean(deltaAmp(deltaAmp<=maxAmp));
    % same for position
    sigmaXyzT = mean(deltaXyz(deltaXyz<=maxXyz));

    % store sigmas
    idlist(1).stats.sigmaAmpT = sigmaAmpT;
    idlist(1).stats.sigmaXyzT = sigmaXyzT;

    % fill distMat. Set -1 if distMat > maxDistance. Make sure that there's
    % something to link at all - otherwise LAP dies
    for t = goodTimesT'
        idlist(t).distMat = sqrt(...
            idlist(t).distMatAmp.^2/sigmaAmpT^2 * constants.relAmpWeight +...
            idlist(t).distMatXyz.^2/sigmaXyzT^2);
        % only check for maxDistance if not -1
        if constants.relativeMaxDistance > 0
            tooLargeDistance = ...
                idlist(t).distMatXyz/sigmaXyzT > constants.relativeMaxDistance;
            if any(tooLargeDistance(:))
                idlist(t).distMat(tooLargeDistance) = -1;
            end
        end
        if constants.absoluteMaxDistance > 0
            tooLargeDistance = ...
                idlist(t).distMatXyz > constants.absoluteMaxDistance;
            if any(tooLargeDistance(:))
                idlist(t).distMat(tooLargeDistance) = -1;
            end
        end

    end
    
    % DEBUG
    if 0
        disp([reshape([idlist(goodTimesT).distMat;goodTimesT'],[],1),...
            reshape([idlist(goodTimesT).distMatAmp;goodTimesT'*sigmaAmpT],[],1)/sigmaAmpT,...
            reshape([idlist(goodTimesT).distMatXyz;goodTimesT'*sigmaXyzT],[],1)/sigmaXyzT])
    end
    
    % loop through the goodTimesT and LAP. Assign only tagIndices so far, no
    % linkup, linkdown. We'll do this later when we know all the links.
    % Also, remember tagIndices.
    for t = goodTimesT'

        % read sourceIdx
        sourceIdxList = idlist(t).sourceIdxList;

        % make sure there are non nans in the distMat, and do LAP
        dm = idlist(t).distMat;
        dm(~isfinite(dm)) = -1;
        [up,down] = lap(dm,-1,0,1);

        % down is where to we link from target.
        down = down(1:nSpots(t));
        down(down > length(sourceIdxList)) = -1;
        oldDown = down > 0;
        newDown = down < 0;

        % write the old tagIndices
        idlist(t).linklist(oldDown,4) = ...
            sourceIdxList(down(oldDown));
        newMaxTagIdx = maxTagIdx + sum(newDown);
        % write the new tagIndices
        idlist(t).linklist(newDown,4) = [maxTagIdx + 1 : newMaxTagIdx]';
        % These are not estimated positions!
        idlist(t).linklist(newDown,3) = 0;
        % store maxIdx, tagIndices
        maxTagIdx = newMaxTagIdx;
        tagIndices(t,1:nSpots(t)) = idlist(t).linklist(:,4)';
        intList(t,idlist(t).linklist(:,4)') = idlist(t).linklist(:,8)';
    end % for goodTimesT: Link

    % plot new links in intensity plot
    if verbose
        plot(intAx,intList,'-o');
    end

end % if isempty(goodTimesT)


% Before going back, sort and augment linklist, however, and fill in linkup and
% linkdown
for t=[goodTimes';NaN,goodTimes(1:end-1)']
    ll = idlist(t(1)).linklist;
    idlist(t(1)).linklist = ...
        [repmat(t(1),maxTagIdx,1),zeros(maxTagIdx,size(ll,2)-1)];
    idlist(t(1)).linklist(ll(:,4),2:end) = ll(:,2:end);
    % write indices for all tags
    idlist(t(1)).linklist(:,4) = (1:maxTagIdx)';
    % keep spot number of all not-found tags at 0

    % linkup, linkdown are nothing else than the spot numbers in the
    % previous frame for a tag-sorted linklist. A link to 0 means that we
    % have effectively a birth or death
    if ~isnan(t(2))
        % we can do linkup to t-1, linkdown for t-1
        idlist(t(1)).linklist(:,6) = idlist(t(2)).linklist(:,2);
        idlist(t(2)).linklist(:,7) = idlist(t(1)).linklist(:,2);
    end
end