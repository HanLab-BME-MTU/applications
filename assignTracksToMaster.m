
% Note: it is assumed that the first channel was used for tracking in each track set
% Both structures must contain single-segment tracks only


function [idxMap unassignedSlave unassignedMaster] = assignTracksToMaster(slaveTracks, masterTracks, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('slaveTracks', @(t) isstruct(t) && all([t.nSeg]==1));
ip.addRequired('masterTracks', @(t) isstruct(t) && all([t.nSeg]==1));
ip.addParamValue('MinOverlap', 1, @isscalar);
ip.addParamValue('MaxDistance', 10, @isscalar);
ip.parse(slaveTracks, masterTracks, varargin{:});
minOverlap = ip.Results.MinOverlap;
R = ip.Results.MaxDistance;

% mean positions
mMeans = arrayfun(@(t) [mean(t.x(1,:)) mean(t.y(1,:))], masterTracks, 'UniformOutput', false);
mMeans = vertcat(mMeans{:});
sMeans = arrayfun(@(t) [mean(t.x(1,:)) mean(t.y(1,:))], slaveTracks, 'UniformOutput', false);
sMeans = vertcat(sMeans{:});

% slave tracks in vicinity to master tracks
idxMap = KDTreeBallQuery(sMeans, mMeans, R);

masterLengths = arrayfun(@(i) numel(i.t), masterTracks);
slaveLengths = arrayfun(@(i) numel(i.t), slaveTracks);

% parse each set of assignments and check for overlap
nm = numel(masterTracks);
for k = 1:nm
    if ~isempty(idxMap{k})
        overlapStart = max([slaveTracks(idxMap{k}).start], masterTracks(k).start);
        overlapEnd = min([slaveTracks(idxMap{k}).end], masterTracks(k).end);
        overlap = overlapEnd - overlapStart + 1;
        % for positive overlaps, compute average distance
        sel = find(overlap>=minOverlap);
        idxMap{k} = idxMap{k}(sel);
        if ~isempty(sel)
            %
            %overlapStart = overlapStart(sel);
            %overlapEnd = overlapEnd(sel);
            dist = zeros(1,numel(sel));
            for i = 1:numel(sel)
                masterIdx = masterTracks(k).f>=overlapStart(sel(i)) & masterTracks(k).f<=overlapEnd(sel(i));
                slaveIdx = slaveTracks(idxMap{k}(i)).f>=overlapStart(sel(i)) & slaveTracks(idxMap{k}(i)).f<=overlapEnd(sel(i));
                dist(i) = min(sqrt( (masterTracks(k).x(1,masterIdx)-slaveTracks(idxMap{k}(i)).x(1,slaveIdx)).^2 +...
                                (masterTracks(k).y(1,masterIdx)-slaveTracks(idxMap{k}(i)).y(1,slaveIdx)).^2));
            end
            idxMap{k} = idxMap{k}(dist<=5);
        end
        %idxMap{k} = idxMap{k}((masterLengths(k) >= slaveLengths(idxMap{k})));
    end
end

unassignedSlave = setdiff(1:numel(slaveTracks), vertcat(idxMap{:}));
unassignedMaster = find(cellfun(@(i) isempty(i), idxMap))';
