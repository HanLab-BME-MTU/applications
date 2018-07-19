%M = catTrackFields(data, tracks, fieldname, ch) concatenates a selected field from all tracks into a matrix

% Francois Aguet, 06/24/11 (Last modified: 09/20/11)

function [M segStarts segEnds seg2trackIndex track2segIndex] = catTrackFields(tracks, movieLength, fieldname, ch)

if nargin<4
    ch = 1;
end

% number of segments in each track
ns = [tracks.nSeg];
M = NaN(sum(ns), movieLength);

segStarts = [tracks(:).segmentStarts];
segEnds = [tracks(:).segmentEnds];

i = 1;
nt = length(tracks);
seg2trackIndex = cell(1,nt);
track2segIndex = cell(1,nt);
for t = 1:nt
    seg2trackIndex{t} = t*ones(1,ns(t));
    track2segIndex{t} = (1:ns(t))+i-1;
    for s = 1:ns(t)
        M(i, segStarts(i):segEnds(i)) = tracks(t).(fieldname){s}(ch,:);
        i = i+1;
    end
end
seg2trackIndex = [seg2trackIndex{:}];
