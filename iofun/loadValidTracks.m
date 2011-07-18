% Francois Aguet, July 18 2011

function tracks = loadValidTracks(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) isstruct(x) & numel(x)==1);
ip.addParamValue('FileName', 'trackAnalysis.mat', @ischar); 
ip.addParamValue('Cutoff', 4, @isscalar);
ip.addParamValue('Sort', 'on', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));

ip.parse(data, varargin{:});

cutoff = ip.Results.Cutoff * data.framerate;

ta = load([data.source 'Tracking' filesep ip.Results.FileName]);
tracks = ta.tracks;
tracks = tracks([tracks.valid]==1 & [tracks.lifetime_s] >= cutoff);

if strcmpi(ip.Results.Sort, 'on')
    [~, sortIdx] = sort([tracks.lifetime_s], 'descend');
    tracks = tracks(sortIdx);
end