% Francois Aguet, July 18 2011

function [tracks trackInfo] = loadTracks(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) isstruct(x) & numel(x)==1);
ip.addParamValue('FileName', 'trackAnalysis.mat', @ischar); 
ip.addParamValue('Cutoff', 4, @isscalar);
ip.addParamValue('Sort', true, @islogical);
% ip.addParamValue('PostProc', [], @isscalar);
ip.addParamValue('Type', 'valid', @(x) any(strcmpi(x, {'valid', 'all'})));
ip.parse(data, varargin{:});

cutoff = ip.Results.Cutoff * data.framerate;

load([data.source 'Tracking' filesep ip.Results.FileName]);

if ip.Results.Sort
    [~, sortIdx] = sort([tracks.lifetime_s], 'descend');
    tracks = tracks(sortIdx);
    trackInfo.track2segIndex = trackInfo.track2segIndex(sortIdx);
end


if strcmpi(ip.Results.Type, 'valid')
    idx = [tracks.valid]==1 & [tracks.lifetime_s] >= cutoff & [tracks.nSeg]==1;
    
%     if ~isempty(ip.Results.PostProc)
%         kLevel = norminv(1-0.05/2.0, 0, 1); % ~2 std above background
%         sb = arrayfun(@(t) sum(sum(t.startBuffer.A(1,:) > t.startBuffer.sigma_r(1,:)*kLevel))>1, tracks);
%         eb = arrayfun(@(t) sum(sum(t.endBuffer.A(1,:) > t.endBuffer.sigma_r(1,:)*kLevel))>1, tracks);
%         tracks(sb | eb) = [];
%         % threshold max intensity
%         maxRatio = arrayfun(@(t) max(t.A(1,:) ./ (t.sigma_r(1,:)*kLevel)), tracks);
%         tracks(maxRatio < ip.Results.PostProc) = [];
%         
%     end
else % if all
    idx = [tracks.lifetime_s] >= cutoff;
end
tracks = tracks(idx);






if nargout>1
    
    % segments to retain
    sidx = trackInfo.track2segIndex(idx);
    sidx = [sidx{:}];
    
    [gapMap, segStarts, segEnds, seg2trackIndex, track2segIndex] = catTrackFields(tracks, data.movieLength, 'gapVect', 1);
   
    trackInfo.x = trackInfo.x(sidx,:);
    trackInfo.y = trackInfo.y(sidx,:);
    trackInfo.gapMap = gapMap;
    trackInfo.segStarts = segStarts;
    trackInfo.segEnds = segEnds;
    trackInfo.seg2trackIndex = seg2trackIndex;
    trackInfo.track2segIndex = track2segIndex;
    
    trackInfo.nSeg = [tracks.nSeg];
    trackInfo.status = [tracks.status];
    trackInfo.valid = [tracks.valid];
    
    trackInfo.lifetimes_f = [tracks.end]-[tracks.start]+1;
end

