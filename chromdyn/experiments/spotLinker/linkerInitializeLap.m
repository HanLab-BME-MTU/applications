function [idlist, goodTimesM, t1t2, intAx] = linkerInitializeLap(idlist, nSpots, ampList, goodIdx, nTimepoints, constants, verbose)
%linkerInitializeLap calculates distanceMatrices and cutoffs for the subsequent LAP
%
% INPUT  verbose - 0 none, 1 intFigure, 2 histograms
%
% OUTPUT idlist - idlist with stored intensity and distance stats and the
%                 distance matrices
%        goodTimesM - timepoints where the correct number has been found
%        t1t2 - pairs of consecutive good timepoints between which the
%               tags will be linked
%        intAx - axes of the intensity plot
%           

% here, we need to
% (1) find bleaching (we know drift already)
% (2) build distMat of corrected intensities and corrected distances
% (3) calculate sigmaInt, sigmaPos from corrected data

% since we have to calculate the distance matrices here already, we can
% store them for re-use


% --- (1) ---
% the intensity fit has to return the time-constant. Do robust weighted
% least squares. Allow seperate I0 for different numbers of tags, but
% require the same time constant for all the fits.

% log(y) = a + b*x, weight: y
% build design matrices so that we minimize (Ax-B)'*W*(Ax-B)
% A: [1, t]. Since there are several numbers of spots, we need a
% "one-column" for every n.
% Unique makes sure that there is no zero-column and we also don't want to
% have entries for n=0
A = [nSpots(goodIdx) * (1./unique(nSpots(goodIdx))') == 1, find(goodIdx)];
[xFit, stdX, goodRows, intAx] = robustExponentialFit2(ampList(goodIdx),A,verbose);

tBleach = xFit(end);
badRows = setdiff(1:nTimepoints,goodRows);

if verbose
    set(get(intAx,'Parent'),...
        'Name',sprintf('Amplitude fit for %s',constants.name));
end


% store stuff so that we can recreate the int-Figure at a later time
idlist(1).stats.intFit.xFit = xFit; % xFit(2) is the time constant (in 1/s)
idlist(1).stats.intFit.badRows = badRows;
idlist(1).stats.ampList = ampList;

% --- (2) ---
% calculate corrected distanceMatrices. Store at t1.
% store all distances and deltaAmps in two additional vectors for
% histograms
% Since we do a two-step linking, we need two sets of distanceMatrices: One
% linking the "sources", one linking "sources" to "targets". However, to
% link the "targets", we need to have linked the sources first.

% Sources
maxNumSpots = min(max(nSpots),constants.MAXSPOTS);
goodTimesM = 0;
delta = 0;
if nnz(goodIdx) > 2
    while length(goodTimesM) < 2 && maxNumSpots-delta > 0
        goodTimesM = find(nSpots == maxNumSpots-delta);
        delta = delta + 1;
    end
    
    % test that we are not at zero spots now
    if maxNumSpots-delta == 0
        % take the lowest nspots with the highest count
        [uniqueEntries,numberOfOccurences] = countEntries(nSpots(nSpots>0));
        [dummy,idx] = max(numberOfOccurences);
        maxNumSpots = uniqueEntries(idx);
        goodTimesM = find(nSpots == maxNumSpots);
    end
else
    [dummy,goodTimesM] = max(nSpots);
end
t1t2 = [goodTimesM(1:end-1)';goodTimesM(2:end)'];

% there are t1*t2 deltas per comparison
[deltaAmp, deltaXyz] = deal(zeros(sum(prod(nSpots(t1t2),1)),1));

deltaCounter = 0;
for t = t1t2
    % t contains [t1;t2]

    % correct amplitude by dividing the second time by exp(tBleach*dt)
    ampFact = exp(tBleach * diff(t));
    idlist(t(1)).distMatAmp = distMat2(...
        idlist(t(1)).linklist(:,8),idlist(t(2)).linklist(:,8) ./ ampFact);

    % correct distMat by subtracting COM from both sets of coordinates.
    % Whenever the number of spots changes, this can be fatal. Therefore,
    % don't correct if change of nSpots - not applicably anymore, b/c only
    % sources are being linked

    %if nSpots(t(1)) == nSpots(t(2))
    idlist(t(1)).distMatXyz = distMat2(...
        idlist(t(1)).linklist(:,9:11) - ...
        repmat(idlist(t(1)).centroid .* constants.useCOM, nSpots(t(1)),1),...
        idlist(t(2)).linklist(:,9:11) - ...
        repmat(idlist(t(2)).centroid .* constants.useCOM, nSpots(t(2)),1));
    %     else
    %          idlist(t(1)).distMatXyz = distMat2(...
    %         idlist(t(1)).linklist(:,9:11),...
    %         idlist(t(2)).linklist(:,9:11));
    %     end

    % store deltas in lists
    deltaCounterEnd = deltaCounter + prod(nSpots(t));
    deltaAmp(deltaCounter+1:deltaCounterEnd) = idlist(t(1)).distMatAmp(:);
    deltaXyz(deltaCounter+1:deltaCounterEnd) = idlist(t(1)).distMatXyz(:);
    deltaCounter = deltaCounterEnd;
end


% --- (3) ---
% find cutoff
[dummy, maxAmp] = cutFirstHistMode(deltaAmp,verbose>1);
if verbose > 1
    set(gcf,'Name','Amplitude Histogram')
end
[dummy, maxXyz] = cutFirstHistMode(deltaXyz,verbose>1);
if verbose > 1
    set(gcf,'Name','Displacement Histogram')
end

% the average deltaAmplitude is the mean of all deltaAmps that are below
% the cutoff
sigmaAmpS = mean(deltaAmp(deltaAmp<=maxAmp));
% same for position
sigmaXyzS = mean(deltaXyz(deltaXyz<=maxXyz));

% fill distMat. Set -1 if distMat > maxDistance
for t = goodTimesM'
    idlist(t).distMat = sqrt(...
        idlist(t).distMatAmp.^2/sigmaAmpS^2 * constants.relAmpWeight +...
        idlist(t).distMatXyz.^2/sigmaXyzS^2);
    % only check for maxDistance if not -1
    if constants.relativeMaxDistance > 0
        tooLargeDistance = ...
            idlist(t).distMatXyz/sigmaXyzS > constants.relativeMaxDistance;
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

% store sigmas
idlist(1).stats.sigmaAmpS = sigmaAmpS;
idlist(1).stats.sigmaXyzS = sigmaXyzS;