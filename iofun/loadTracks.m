%tracks = loadTracks(data, varargin)
%
% Inputs:
%
%  data      : data structure
% 
% Options:
%
% 'Category' : 'Ia'  Single tracks with valid gaps
%              'Ib'  Single tracks with invalid gaps
%              'Ic'  Single tracks cut at beginning or end
%              'Id'  Single tracks, persistent
%              'IIa' Compound tracks with valid gaps
%              'IIb' Compound tracks with invalid gaps
%              'IIc' Compound tracks cut at beginning or end
%              'IId' Compound tracks, persistent

% Francois Aguet (last modified: 02/06/12)

function tracks = loadTracks(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) isstruct(x) & numel(x)==1);
ip.addParamValue('FileName', 'trackAnalysis.mat', @ischar); 
ip.addParamValue('Cutoff', 4, @isscalar);
ip.addParamValue('Sort', true, @islogical);
% ip.addParamValue('PostProc', [], @isscalar);
ip.addParamValue('Category', 'all', @(x) any(strcmpi(x, {'all', 'Ia', 'Ib', 'Ic', 'Id', 'IIa', 'IIb', 'IIc', 'IId'})));
ip.parse(data, varargin{:});

cutoff_s = ip.Results.Cutoff * data.framerate;

load([data.source 'Tracking' filesep ip.Results.FileName]);

if ip.Results.Sort
    [~, sortIdx] = sort([tracks.lifetime_s], 'descend'); %#ok<NODEF>
    tracks = tracks(sortIdx);
end

singleIdx = [tracks.nSeg]==1;
validGaps = arrayfun(@(t) max([t.gapStatus 4]), tracks)==4;
vis = [tracks.visibility];
switch ip.Results.Type
    case 'Ia'
        idx = singleIdx & validGaps & vis==1;
    case 'Ib'
        idx = singleIdx & ~validGaps & vis==1;
    case 'Ic'
        idx = singleIdx & vis==2;
    case 'Id'
        idx = singleIdx & vis==3;
    case 'IIa'
        idx = ~singleIdx & validGaps & vis==1;
    case 'IIb'
        idx = ~singleIdx & ~validGaps & vis==1;
    case 'IIc'
        idx = ~singleIdx & vis==2;
    case 'IId'
        idx = ~singleIdx & vis==3;
    case 'all'
        idx = 1:numel(tracks);
end
idx = idx & [tracks.lifetime_s] >= cutoff_s;
tracks = tracks(idx);
    
%     if ~isempty(ip.Results.PostProc)
%         kLevel = norminv(1-0.05/2.0, 0, 1); % ~2 std above background
%         sb = arrayfun(@(t) sum(sum(t.startBuffer.A(1,:) > t.startBuffer.sigma_r(1,:)*kLevel))>1, tracks);
%         eb = arrayfun(@(t) sum(sum(t.endBuffer.A(1,:) > t.endBuffer.sigma_r(1,:)*kLevel))>1, tracks);
%         tracks(sb | eb) = [];
%         % threshold max intensity
%         maxRatio = arrayfun(@(t) max(t.A(1,:) ./ (t.sigma_r(1,:)*kLevel)), tracks);
%         tracks(maxRatio < ip.Results.PostProc) = [];
%     end
