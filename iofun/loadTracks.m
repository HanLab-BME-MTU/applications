% Francois Aguet, July 18 2011

function tracks = loadTracks(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) isstruct(x) & numel(x)==1);
ip.addParamValue('FileName', 'trackAnalysis.mat', @ischar); 
ip.addParamValue('Cutoff', 4, @isscalar);
ip.addParamValue('Sort', 'on', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.addParamValue('Type', 'valid', @(x) strcmpi(x, 'valid') | strcmpi(x, 'all'));
ip.parse(data, varargin{:});

cutoff = ip.Results.Cutoff * data.framerate;

ta = load([data.source 'Tracking' filesep ip.Results.FileName]);
tracks = ta.tracks;

if strcmpi(ip.Results.Type, 'valid')
    tracks = tracks([tracks.valid]==1 & [tracks.lifetime_s] >= cutoff & arrayfun(@(t) ~iscell(t.x), tracks));
else % if all
    tracks = tracks([tracks.lifetime_s] >= cutoff);
end

if strcmpi(ip.Results.Sort, 'on')
    rTrackIdx = find(arrayfun(@(t) ~iscell(t.x), tracks));
    mTrackIdx = setdiff(1:length(tracks), rTrackIdx);
    [~, rSortIdx] = sort([tracks(rTrackIdx).lifetime_s], 'descend');
    [~, mSortIdx] = sort([tracks(mTrackIdx).lifetime_s], 'descend');
    tracks = tracks([rSortIdx mSortIdx+rTrackIdx(end)]);
end