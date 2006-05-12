function [idlist] = linkerEstimatePositions(idlist, maxTagIdx, goodTimes, constants, verbose, intAx)
%likerEstimatePositions estimates positions, amplitudes and serchradii for
%lost tags

% loop through tagIndices. Find where we have interruptions. Loop through
% interruptions. Find start and end (if exist). Estimate XYZ-position and
% set search radius for every frame and fill in.

for i = 1:maxTagIdx

    % collect timepoint-#, amplitude and positions
    timeAmpPos = catStruct(1,sprintf('idlist.linklist(%i,[1,8:11])',i));

    % This list only contains entries for good timepoints.
    % Find indices where we found the tag
    goodTagIdx = find(timeAmpPos(:,2));
    badTagIdx  = ~timeAmpPos(:,2); % badTagIdx is logical

    % use goodTagIdxPos to estimate positions (take fusions into account!)
    goodTagIdxPos = find(timeAmpPos(:,3));
    badTagIdxPos =  ~timeAmpPos(:,3); % badTagIdx is logical

    % Fit the amplitudes. We already know the time constant, and only need
    % to find the multiplier. Therefore, we can use robustMean.
    % We will replace the amplitudes of outliers with their fitted amps.
    % of course, only fit if there is more than one point...

    expValue = ...
        exp(idlist(1).stats.intFit.xFit(end) * timeAmpPos(goodTagIdx,1));
    if size(goodTagIdx,1) > 1
        [factor, dummy, inlierIdx] = ...
            robustMean(timeAmpPos(goodTagIdx,2)./expValue);
        % plot if verbose
        if verbose
            plotMatrix = zeros(nnz(badTagIdx),maxTagIdx);
            plotMatrix(:,i) = ...
                exp(idlist(1).stats.intFit.xFit(end) * ...
                timeAmpPos(badTagIdx,1)) * factor;
            plot(intAx,timeAmpPos(badTagIdx,1),plotMatrix,'d');%,...
            %'MarkerFaceColor',...
            %colorOrder(wraparound(i,[1;size(colorOrder,1)]),:))
        end
    else
        factor = timeAmpPos(goodTagIdx,2)./expValue;
        inlierIdx = 1;
        % if a tag only appears once, it probably signifies a superfluous
        % spot. Therefore, flag with 2/3 - we don't want to have the index
        % appear everywhere, and we probably want to delete it.
        for t = goodTimes'
            idlist(t).linklist(i,5) = 3;
        end
        % of course, we want to show the tag where it appears
        idlist(timeAmpPos(goodTagIdx,1)).linklist(i,5) = 2;

    end

    % store int-fit in idlist.stats
    idlist(1).stats.intFit.tagFactor(i) = factor;

    if constants.replaceAmplitudes
        outlierIdx = missingIndices(inlierIdx,length(goodTagIdx));

        % replace amplitudes and set flag 1
        if ~isempty(outlierIdx)
            for o=outlierIdx
                t = timeAmpPos(goodTagIdx(o),1);
                idlist(t).linklist(i,8) = expValue(o) * factor;
                idlist(t).linklist(i,5) = 1;
                if verbose
                    plotMatrix = zeros(1,maxTagIdx);
                    plotMatrix(i) = idlist(t).linklist(i,8);
                    plot(intAx,t,plotMatrix,'s');
                end
            end
        end
    end
    % write estimated amplitudes
    timeAmpPos(badTagIdx,2) = ...
        factor * ...
        exp(timeAmpPos(badTagIdx,1) * idlist(1).stats.intFit.xFit(end));

    % because of fusions, we write amplitudes separately
    for t = 1:length(badTagIdx)
        idlist(timeAmpPos(t,1)).linklist(i,8) = timeAmpPos(t,2);
    end


    % find improved sigmaPos, if there are enough tag positions
    % -- later --
    % subtract COM
    timeAmpPos(goodTagIdxPos,3:5) = timeAmpPos(goodTagIdxPos,3:5) - ...
        cat(1,idlist(timeAmpPos(goodTagIdxPos,1)).centroid) * constants.useCOM;

    % bwlabel badTagIdx, so that we identify interruptions
    badGroups = bwlabel(badTagIdxPos);

    % loop through interruptions. Find first and last frame (and good frame
    % before and/or after). Then loop through frames and fill in estimated
    % position. Also write idlist.trackInit as [tagIdx, xDelta, yDelta,
    % zDelta]. xyz are in pixels.
    for j=1:max(badGroups)
        currentGroupIdx = find(badGroups==j);
        startEnd = [0,0]; % startEnd is an index into timeAmpPos
        startEnd(1) = currentGroupIdx(1) - 1; % remains 0 if cgi(1)==1
        if currentGroupIdx(end) == size(timeAmpPos,1)
            % startEnd(2) remains 0
        else
            startEnd(2) = currentGroupIdx(end) + 1;
        end

        switch (startEnd(1) == 0) + 2 * (startEnd(2) == 0)
            case 0 % we have a start and an end
                pos1 = timeAmpPos(startEnd(1),3:5);
                pos2 = timeAmpPos(startEnd(2),3:5);

                startT = timeAmpPos(startEnd(1),1);
                deltaT = timeAmpPos(startEnd(2),1) - startT;
                % loop through interruptions. Calculate stuff
                for tapIdx = startEnd(1)+1:startEnd(2)-1
                    tIdlist = timeAmpPos(tapIdx,1);
                    tau1 = tIdlist - startT;



                    % calculate position and fill linklist. Remember to add
                    % COM back again. Also, write lostChi2
                    newPos = (pos1*(deltaT-tau1) + pos2*tau1)/deltaT + ...
                        idlist(tIdlist).centroid .* constants.useCOM;
                    idlist(tIdlist).linklist(i,9:12) = ...
                        [newPos, constants.lostChi2];
                    % flag spot as estimated position
                    if idlist(tIdlist).linklist(i,3) == 0
                        idlist(tIdlist).linklist(i,3) = 1;
                    end
                    % calculate delta for searchRadius. Convert to pixels
                    searchRadius = constants.diffusionFactor * ...
                        idlist(1).stats.sigmaXyzS * ...
                        sqrt(tau1-tau1^2/deltaT) ./ constants.pix2mu;
                    idlist(tIdlist).trackInit = ...
                        [idlist(tIdlist).trackInit;...
                        i, searchRadius];

                end
            case 1 % only end (at beginning of movie)
                pos2 = timeAmpPos(startEnd(2),3:5);

                startT = 0;
                deltaT = timeAmpPos(startEnd(2),1) - startT;
                % loop through interruptions. Calculate stuff
                for tapIdx = startEnd(1)+1:startEnd(2)-1
                    tIdlist = timeAmpPos(tapIdx,1);
                    tau2 = deltaT - tIdlist;



                    % We keep pos2. Remember to add
                    % COM back again
                    newPos = pos2 + ...
                        idlist(tIdlist).centroid .* constants.useCOM;
                    idlist(tIdlist).linklist(i,9:12) = ...
                        [newPos, constants.lostChi2];
                    % flag spot as estimated position
                    if idlist(tIdlist).linklist(i,3) == 0
                        idlist(tIdlist).linklist(i,3) = 1;
                    end
                    % calculate delta for searchRadius. Convert to pixels
                    searchRadius = constants.diffusionFactor * ...
                        idlist(1).stats.sigmaXyzS * ...
                        sqrt(tau2) ./ constants.pix2mu;
                    idlist(tIdlist).trackInit = ...
                        [idlist(tIdlist).trackInit;...
                        i, searchRadius];

                end
            case 2 % only start (at end of movie)
                pos1 = timeAmpPos(startEnd(1),3:5);

                startT = timeAmpPos(startEnd(1),1);
                %deltaT = max(goodTimes) - startT;

                % loop through interruptions. Calculate stuff
                for tapIdx = startEnd(1)+1:size(timeAmpPos,1)
                    tIdlist = timeAmpPos(tapIdx,1);
                    tau1 = tIdlist - startT;



                    % Only pos1. Remember to add
                    % COM back again
                    newPos = pos1 + ...
                        idlist(tIdlist).centroid .* constants.useCOM;
                    idlist(tIdlist).linklist(i,9:12) = ...
                        [newPos, constants.lostChi2];
                    % flag spot as estimated position
                    if idlist(tIdlist).linklist(i,3) == 0
                        idlist(tIdlist).linklist(i,3) = 1;
                    end
                    % calculate delta for searchRadius. Convert to pixels
                    searchRadius = constants.diffusionFactor * ...
                        idlist(1).stats.sigmaXyzS * ...
                        sqrt(tau1) ./ constants.pix2mu;
                    idlist(tIdlist).trackInit = ...
                        [idlist(tIdlist).trackInit;...
                        i, searchRadius];

                end
            otherwise
                error('One-frame-movie or major bug in linker')
        end
    end
end