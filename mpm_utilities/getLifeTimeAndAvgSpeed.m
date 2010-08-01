function [lifeTime avgSpeed ADD_EVERY_ELEMENTS_ON_THE_BOARD] = getLifeTimeAndAvgSpeed(MPM,distToEdge,dLims)

[nrows ncols] = size(MPM);
nFrames = ncols / 2;
nBands = numel(dLims)-1;

%% Compute the distance of each track point away from the cell edge.

% track foot print
trackMask = mpm2trackMask(MPM,true);
ind = find(trackMask);
[iRow, t] = ind2sub(size(trackMask),ind);

% y-coordinate of each track point
trackY = zeros(size(trackMask));
trackY(ind) = MPM(sub2ind(size(MPM), iRow, 2 * t - 1));

% x-coordinate of each track point
trackX = zeros(size(trackMask));
trackX(ind) = MPM(sub2ind(size(MPM), iRow, 2 * t));

% Distance of each track point away from the cell edge
ind2 = sub2ind(size(distToEdge),trackY(ind),trackX(ind),t);
trackDistPoints = zeros(size(trackMask));
trackDistPoints(ind) = distToEdge(ind2);

%% Compute the pairwise distance between consecutive track points
% A single point track has a distance equals to 0. A two-point track has a
% distance equals to the distance between the 2 points, etc.
trackPWDPoints = zeros(size(trackMask));
trackX2 = [trackX(:,2:end), zeros(nrows,1)];
trackY2 = [trackY(:,2:end), zeros(nrows,1)];
trackPWDPoints(ind) = sqrt((trackX(ind) - trackX2(ind)).^2 + (trackY(ind) - trackY2(ind)).^2);
% clear the value at the last frame of each track
trackPWDPoints(trackMask & ~([trackMask(:,2:end) false(nrows,1)])) = 0;

%% Compute lifetime and avg speed
trackID = (1:nrows)';

accu = zeros(size(trackID));
accuPWD = zeros(size(trackID));

inBand = false(size(trackID,1), nBands);
    
lifeTime = cell(nFrames,nBands);
avgSpeed = cell(nFrames,nBands);

for iFrame = 1:nFrames
    idxLive = trackID(trackMask(:,iFrame));
    idxDead = trackID(~trackMask(:,iFrame));
    
    % accumulate lifetime
    accu(idxLive) = accu(idxLive) + 1;
    
    % accumulate pair-wise distance
    accuPWD(idxLive) = accuPWD(idxLive) + trackPWDPoints(idxLive,iFrame);
    
    for iBand = 1:nBands
        % accumulate iBand
        inBand(idxLive,iBand) = inBand(idxLive,iBand) | ...
            trackDistPoints(idxLive,iFrame) >= dLims(iBand) & ...
            trackDistPoints(idxLive,iFrame) <= dLims(iBand+1);
        
        % Indices to Keep
        ind2keep = ~trackMask(:,iFrame) & accu > 1 & inBand(:,iBand);
        
        lifeTime{iFrame,iBand} = accu(ind2keep);
        
        avgSpeed{iFrame,iBand} = (1 ./ accu(ind2keep)) .* accuPWD(ind2keep);
        
        % Reset
        inBand(idxDead,iBand) = false;
    end
    
    % Reset
    accu(idxDead) = 0;
    accuPWD(idxDead) = 0;
end

lifeTime = arrayfun(@(iBand) vertcat(lifeTime{:,iBand}), 1:nBands, 'UniformOutput', false);
avgSpeed = arrayfun(@(iBand) vertcat(avgSpeed{:,iBand}), 1:nBands, 'UniformOutput', false);
