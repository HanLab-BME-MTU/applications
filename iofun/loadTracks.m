% Francois Aguet, July 18 2011

function tracks = loadTracks(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) isstruct(x) & numel(x)==1);
ip.addParamValue('FileName', 'trackAnalysis.mat', @ischar); 
ip.addParamValue('Cutoff', 4, @isscalar);
ip.addParamValue('Sort', 'on', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.addParamValue('PostProc', [], @isscalar);
ip.addParamValue('Type', 'valid', @(x) strcmpi(x, 'valid') | strcmpi(x, 'all'));
ip.parse(data, varargin{:});

cutoff = ip.Results.Cutoff * data.framerate;

ta = load([data.source 'Tracking' filesep ip.Results.FileName]);
tracks = ta.tracks;

if strcmpi(ip.Results.Type, 'valid')
    tracks = tracks([tracks.valid]==1 & [tracks.lifetime_s] >= cutoff & arrayfun(@(t) ~iscell(t.x), tracks));
    if ~isempty(ip.Results.PostProc)
        kLevel = norminv(1-0.05/2.0, 0, 1); % ~2 std above background
        sb = arrayfun(@(t) sum(sum(t.startBuffer.A(1,:) > t.startBuffer.sigma_r(1,:)*kLevel))>1, tracks);
        eb = arrayfun(@(t) sum(sum(t.endBuffer.A(1,:) > t.endBuffer.sigma_r(1,:)*kLevel))>1, tracks);
        tracks(sb | eb) = [];
        % threshold max intensity
        maxRatio = arrayfun(@(t) max(t.A(1,:) ./ (t.sigma_r(1,:)*kLevel)), tracks);
        tracks(maxRatio < ip.Results.PostProc) = [];
        
    end
else % if all
    tracks = tracks([tracks.lifetime_s] >= cutoff);
end

if strcmpi(ip.Results.Sort, 'on')
    rTrackIdx = find([tracks.type]==1);%find(arrayfun(@(t) ~iscell(t.x), tracks));
    mTrackIdx = setdiff(1:length(tracks), rTrackIdx);
    [~, rSortIdx] = sort([tracks(rTrackIdx).lifetime_s], 'descend');
    [~, mSortIdx] = sort([tracks(mTrackIdx).lifetime_s], 'descend');
    tracks = tracks([rSortIdx mSortIdx+rTrackIdx(end)]);
end