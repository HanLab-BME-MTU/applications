
% Note: it is assumed that the first channel was used for tracking in each track set
% Both structure must contain single segment tracks only


function idx = assignTracksToMaster(slaveTracks, masterTracks, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('slaveTracks', @(t) isstruct(t) && all([t.nSeg]==1));
ip.addRequired('masterTracks', @(t) isstruct(t) && all([t.nSeg]==1));
ip.addParamValue('MinOverlap', 1, @isscalar);
ip.addParamValue('MaxDistance', 5, @isscalar);
ip.parse(slaveTracks, masterTracks, varargin{:});
minOverlap = ip.Results.MinOverlap;
R = ip.Results.MaxDistance;

% mean positions
mMeans = arrayfun(@(t) [mean(t.x{1}(1,:)) mean(t.y{1}(1,:))], masterTracks, 'UniformOutput', false);
mMeans = vertcat(mMeans{:});
sMeans = arrayfun(@(t) [mean(t.x{1}(1,:)) mean(t.y{1}(1,:))], slaveTracks, 'UniformOutput', false);
sMeans = vertcat(sMeans{:});

% tracks overlapping in space
idx = KDTreeBallQuery(sMeans, mMeans, R);

masterLengths = arrayfun(@(t) numel(t.t{1}), masterTracks);
slaveLengths = arrayfun(@(t) numel(t.t{1}), slaveTracks);

% parse each set of assignments and check for overlap
nm = numel(masterTracks);
for k = 1:nm
    if ~isempty(idx{k})
        overlap = min([slaveTracks(idx{k}).end], masterTracks(k).end) - ...
            max([slaveTracks(idx{k}).start], masterTracks(k).start) + 1;
        idx{k} = idx{k}(overlap>=minOverlap & (masterLengths(k) >= slaveLengths(idx{k})));
    end
end