%M = catTrackFields(data, tracks, fieldname, ch) concatenates a selected field from all tracks into a matrix

% Francois Aguet, 06/24/11 (Last modified: 08/12/11)

function [M segStarts segEnds trackIndex] = catTrackFields(tracks, movieLength, fieldname, ch)

if nargin<4
    ch = 1;
end

nt = length(tracks);

% tracks must be sorted as regular followed by merging/splitting
idx = find(arrayfun(@(t) ~iscell(t.x), tracks)==1, 1, 'last');
%idx = find([tracks.type]==1, 1, 'last');
rIdx = 1:idx; % regular tracks
cIdx = (idx+1):nt; % compound tracks

M1 = NaN(idx, movieLength);
starts1 = [tracks(rIdx).start]';
ends1 = [tracks(rIdx).end]';
for t = 1:idx
    M1(t, starts1(t):ends1(t)) = tracks(t).(fieldname)(ch,:);
end

% number of segments in each compound track
ns = arrayfun(@(t) length(t.x), tracks(cIdx));
M2 = NaN(sum(ns), movieLength);

i = 1;
starts2 = NaN(sum(ns),1);
ends2 = NaN(sum(ns),1);
nc = length(cIdx);
trackIndex = cell(1,nc);
for t = 1:nc
    trackIndex{t} = t*ones(1,ns(t))+idx;
    for s = 1:ns(t)
        ti = tracks(cIdx(t));
        starts2(i) = ti.segmentStarts{s};
        ends2(i) = ti.segmentEnds{s};
        M2(i, starts2(i):ends2(i)) = ti.(fieldname){s}(ch,:);
        i = i+1;
    end
end
trackIndex = [1:idx trackIndex{:}];

M = [M1; M2];

if nargout > 1
    segStarts = [starts1; starts2];
    segEnds = [ends1; ends2];
end