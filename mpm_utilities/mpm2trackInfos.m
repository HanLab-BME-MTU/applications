function [trackInfos xMap yMap] = mpm2trackInfos(MPM,distToEdge,dLims,minLifetime,keepCompleteTracksOnly)
% Input values:
%
% MPM:                    the Magic Position Matrix (see qfsm for details)
%
% distToEdge:             the distance transform of the cell's footprint.
%                         It is a matrix of size NxMxT where T is the
%                         number of frames.
%
% dLims:                  a vector of distances away from the cell edge.
%                         These distances have to be ordered.
%
% minLifetime:            minimum number of frames a track needs to live
%                         for. Default is 1.
%
% keepCompleteTracksOnly: a logical value that allows to discard track that
%                         are either beginning at the very first frame of
%                         the movie or ending at the very last frame.
%
% Output values:
%
% tracks:                 is a cell array of size = numel(dLims)-1. Each
%                         cell element i contains a array Nx3 where N is
%                         the number of tracks falling into a band away
%                         from the cell edge defined by dLims(i)
%                         dLims(i+1). Each track contains:
%
%        tracks(i,1) = row index in the MPM
%        tracks(i,2) = first frame
%        tracks(i,3) = last frame
%
% xMap:                   the x-coordinate of track point. xMap is a matrix
%                         which number of rows equals MPM's number of rows
%                         and its number of columns equals half MPM's
%                         number of columns (= number of frames).
%
% yMap:                   the y-coordinate of track point.
%
% Sylvain Berlemont, 2010

if nargin < 4 || isempty(minLifetime)
    minLifetime = 1;
end

if nargin < 5 || isempty(keepCompleteTracksOnly)
    keepCompleteTracksOnly = true;
end

[nrows ncols] = size(MPM);
nFrames = ncols / 2;
nBands = numel(dLims)-1;

%% Compute the distance of each track point away from the cell edge.

% track foot print
trackMask = mpm2trackMask(MPM,keepCompleteTracksOnly);
ind = find(trackMask);
[iRow, t] = ind2sub(size(trackMask),ind);

% y-coordinate of each track point
yMap = zeros(size(trackMask));
yMap(ind) = MPM(sub2ind(size(MPM), iRow, 2 * t - 1));

% x-coordinate of each track point
xMap = zeros(size(trackMask));
xMap(ind) = MPM(sub2ind(size(MPM), iRow, 2 * t));

% Distance of each track point away from the cell edge
ind2 = sub2ind(size(distToEdge),yMap(ind),xMap(ind),t);
trackDistPoints = zeros(size(trackMask));
trackDistPoints(ind) = distToEdge(ind2);

%% Get first and last frames
trackID = (1:nrows)';

accu = zeros(size(trackID));
%accuPWD = zeros(size(trackID));

inBand = false(size(trackID,1), nBands);
   
firstFrame = zeros(size(trackID));

trackInfos = cell(nFrames,nBands);

for iFrame = 1:nFrames
    idxLive = trackID(trackMask(:,iFrame));
    idxDead = trackID(~trackMask(:,iFrame));
    
    % accumulate lifetime
    accu(idxLive) = accu(idxLive) + 1;
    
    % keep first frame
    firstFrame(trackMask(:,iFrame) & firstFrame == 0) = iFrame;

    for iBand = 1:nBands
        % accumulate iBand
        inBand(idxLive,iBand) = inBand(idxLive,iBand) | ...
            trackDistPoints(idxLive,iFrame) >= dLims(iBand) & ...
            trackDistPoints(idxLive,iFrame) <= dLims(iBand+1);
        
        % Indices to Keep
        ind2keep = ~trackMask(:,iFrame) & accu >= minLifetime & inBand(:,iBand);
        
        trackInfos{iFrame,iBand} = horzcat(...
            trackID(ind2keep),...
            firstFrame(ind2keep),...
            repmat(iFrame-1,nnz(ind2keep),1));
            
        % Reset
        inBand(idxDead,iBand) = false;
    end
    
    % Reset
    accu(idxDead) = 0;
    firstFrame(idxDead) = 0;
end

trackInfos = arrayfun(@(iBand) vertcat(trackInfos{:,iBand}), 1:nBands,...
    'UniformOutput', false);
